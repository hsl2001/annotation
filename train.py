#!/usr/bin/env python3
"""Train and evaluate anno on sampled CRF windows with 7-fold cross-validation.

Usage:
    python3 train.py -n 2
    python3 train.py -n 3 -p 8
    python3 train.py -n 2 --fasta truth.nip.chr5.fasta truth.slm.chr5.fasta --gff truth.nip.chr5.gff3 truth.slm.chr5.gff3
"""

import argparse
from collections import defaultdict
from dataclasses import dataclass
import io
import math
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

from accuracy import (
    filtered_gff,
    gene_level,
    read_fasta,
    reverse_complement,
    run_gffcompare,
    scores,
)
from make_dataset import (
    download,
    get_selection,
    log,
    pair_inputs,
    parse_explicit_pairs,
    prefix_for,
    rewrite_attributes,
    source_files,
    write_annotation,
    FASTA_SUFFIXES,
    GFF_SUFFIXES,
)

CRF_WINDOW_SIZE = 65536
DEFAULT_WINDOWS_PER_GENOME = 1000
DEFAULT_CWT_PAD = 1024
DEFAULT_FOLDS = 7


@dataclass(frozen=True)
class SampledWindow:
    genome_index: int
    genome_name: str
    seqid: str
    start: int          # 0-based active window start in contig
    end: int            # 0-based active window end (exclusive) in contig
    ext_start: int      # 0-based padded extract start in contig
    ext_end: int        # 0-based padded extract end (exclusive) in contig
    gene_count: int


def valid_cds_sequence(seq: str, intervals: list[tuple[int, int]], strand: str) -> bool:
    total = sum(end - start + 1 for start, end in intervals)
    if total < 6 or total % 3 != 0:
        return False
    fragments = [seq[start - 1:end] for start, end in intervals]
    coding = "".join(fragments).upper()
    if strand == "-":
        coding = reverse_complement(coding)
    if any(base not in "ACGT" for base in coding):
        return False
    if coding[:3] != "ATG" or coding[-3:] not in {"TAA", "TAG", "TGA"}:
        return False
    for index in range(3, len(coding) - 3, 3):
        if coding[index:index + 3] in {"TAA", "TAG", "TGA"}:
            return False
    return True


def compute_contig_gene_spans(
    sequences: dict[str, str], gff_path: Path
) -> tuple[dict[str, list[tuple[int, int]]], dict[str, list[str]]]:
    """Replicate anno.c label_cds and store normalized GFF lines per contig."""
    lengths = {name: len(seq) for name, seq in sequences.items()}
    normalized_buf = io.StringIO()
    write_annotation(gff_path, normalized_buf, "", lengths)
    normalized_lines = normalized_buf.getvalue().splitlines()

    groups: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    gff_by_seqid: dict[str, list[str]] = defaultdict(list)
    for line in normalized_lines:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        seqid = fields[0]
        if seqid not in lengths:
            continue
        gff_by_seqid[seqid].append(line)
        if fields[2] != "CDS":
            continue
        start, end = int(fields[3]), int(fields[4])
        parent = None
        for attr in fields[8].split(";"):
            if attr.startswith("Parent="):
                parent = attr[7:].split(",")[0]
                break
        if parent:
            groups[parent].append((seqid, fields[6], start, end))

    candidates: list[tuple[int, str, str, str, list[tuple[int, int]]]] = []
    for parent, items in groups.items():
        items.sort(key=lambda item: item[2])
        total = sum(end - start + 1 for _, _, start, end in items)
        seqid, strand = items[0][0], items[0][1]
        intervals = [(start, end) for _, _, start, end in items]
        candidates.append((total, parent, seqid, strand, intervals))
    candidates.sort(key=lambda item: (-item[0], item[1]))

    occupied: dict[str, bytearray] = {name: bytearray(length) for name, length in lengths.items()}
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for _, _, seqid, strand, intervals in candidates:
        if strand not in {"+", "-"}:
            continue
        if any(curr[0] <= prev[1] + 2 for prev, curr in zip(intervals, intervals[1:])):
            continue
        if not valid_cds_sequence(sequences[seqid], intervals, strand):
            continue
        span_start = intervals[0][0] - 1
        span_end = intervals[-1][1]
        label = occupied[seqid]
        if any(label[pos] for pos in range(span_start, span_end)):
            continue
        for pos in range(span_start, span_end):
            label[pos] = 1
        spans[seqid].append((span_start, span_end))

    for seqid in spans:
        spans[seqid].sort()
    return spans, gff_by_seqid


def partition_crf_windows(
    genome_index: int,
    genome_name: str,
    sequences: dict[str, str],
    spans: dict[str, list[tuple[int, int]]],
    cwt_pad: int,
) -> list[SampledWindow]:
    """Partition contigs into CRF windows matching crf.h make_windows."""
    windows: list[SampledWindow] = []
    for seqid, seq in sequences.items():
        length = len(seq)
        contig_spans = spans.get(seqid, [])
        start = 0
        span_idx = 0
        num_spans = len(contig_spans)
        while start < length:
            end = start + CRF_WINDOW_SIZE if length - start > CRF_WINDOW_SIZE else length
            if end < length:
                for gs, ge in contig_spans:
                    if gs < end < ge:
                        end = gs if (gs > start and end - gs <= ge - end) else ge
                        break
                    if gs >= end:
                        break
            if end <= start:
                end = min(length, start + CRF_WINDOW_SIZE)

            while span_idx < num_spans and contig_spans[span_idx][1] <= start:
                span_idx += 1
            k = span_idx
            gene_count = 0
            while k < num_spans and contig_spans[k][0] < end:
                gs, ge = contig_spans[k]
                if gs >= start and ge <= end:
                    gene_count += 1
                k += 1

            ext_start = max(0, start - cwt_pad)
            ext_end = min(length, end + cwt_pad)
            windows.append(
                SampledWindow(
                    genome_index=genome_index,
                    genome_name=genome_name,
                    seqid=seqid,
                    start=start,
                    end=end,
                    ext_start=ext_start,
                    ext_end=ext_end,
                    gene_count=gene_count,
                )
            )
            start = end
    return windows


def select_windows_for_genome(
    all_windows: list[SampledWindow],
    target_count: int,
    rng: random.Random,
) -> list[SampledWindow]:
    """Select up to target_count windows per genome, prioritizing coding windows."""
    coding = [w for w in all_windows if w.gene_count > 0]
    noncoding = [w for w in all_windows if w.gene_count == 0]
    if len(coding) >= target_count:
        selected = rng.sample(coding, target_count)
    else:
        need = min(target_count - len(coding), len(noncoding))
        selected = coding + (rng.sample(noncoding, need) if need > 0 else [])
    selected.sort(key=lambda w: (w.seqid, w.start))
    return selected


def window_contig_id(window: SampledWindow, index: int) -> str:
    safe_seqid = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in window.seqid)
    return f"g{window.genome_index}_w{index:05d}_{safe_seqid}_{window.start + 1}_{window.end}"


def write_window_dataset(
    windows: list[SampledWindow],
    genome_sequences: dict[int, dict[str, str]],
    genome_gff_lines: dict[int, dict[str, list[str]]],
    fasta_path: Path,
    gff_path: Path,
) -> None:
    """Write padded FASTA (with crop=offset:length) and shifted GFF3 for sampled windows."""
    with fasta_path.open("w", encoding="utf-8") as fasta_out, gff_path.open("w", encoding="utf-8") as gff_out:
        gff_out.write("##gff-version 3\n")
        for idx, win in enumerate(windows):
            contig_id = window_contig_id(win, idx)
            full_seq = genome_sequences[win.genome_index][win.seqid]
            subseq = full_seq[win.ext_start:win.ext_end]
            crop_offset = win.start - win.ext_start
            crop_length = win.end - win.start
            fasta_out.write(
                f">{contig_id} crop={crop_offset}:{crop_length} "
                f"genome={win.genome_name} source={win.seqid}:{win.start + 1}-{win.end}\n"
            )
            for pos in range(0, len(subseq), 80):
                fasta_out.write(subseq[pos:pos + 80] + "\n")

            active_start_1 = win.start + 1
            active_end_1 = win.end
            for line in genome_gff_lines[win.genome_index].get(win.seqid, []):
                fields = line.split("\t")
                if fields[2] == "region":
                    continue
                f_start, f_end = int(fields[3]), int(fields[4])
                if f_start < active_start_1 or f_end > active_end_1:
                    continue
                fields[0] = contig_id
                fields[3] = str(f_start - win.ext_start)
                fields[4] = str(f_end - win.ext_start)
                fields[8] = rewrite_attributes(fields[8], f"{contig_id}__")
                gff_out.write("\t".join(fields) + "\n")


def evaluate_predictions(
    reference_gff: Path,
    query_gff: Path,
    genome_fasta_path: Path,
    workdir: Path,
    threads: int,
) -> dict[str, tuple[float, float, float]]:
    """Compute protein gene-level and gffcompare CDS exon-level (precision, recall, f1)."""
    gene_metrics, _ = gene_level(
        reference_gff,
        query_gff,
        genome_fasta_path,
        workdir,
        bitscore=50.0,
        evalue=1e-5,
        coverage=0.70,
        threads=threads,
    )
    genome = {name: seq.upper() for name, seq in read_fasta(genome_fasta_path)}
    ref_longest = workdir / "reference.longest.cds.gff3"
    qry_cds = workdir / "query.cds.gff3"
    ref_count = filtered_gff(reference_gff, genome, ref_longest, longest_only=True)
    qry_count = filtered_gff(query_gff, genome, qry_cds)
    exon_metrics = scores(0, qry_count, ref_count)
    if ref_count and qry_count:
        exon_metrics = run_gffcompare(ref_longest, qry_cds, workdir, "gffcompare_exon")["exon"]
    return {"gene": gene_metrics, "exon": exon_metrics}


def resolve_genome_inputs(args: argparse.Namespace, parser: argparse.ArgumentParser) -> list[tuple[str, Path, Path]]:
    has_pair = bool(args.pair)
    has_fasta_gff = bool(args.fasta or args.gff)
    has_dir = args.directory is not None
    has_ncbi_opts = (
        args.cache_dir is not None
        or args.min_busco is not None
        or args.workers is not None
    )

    local_sources = sum((has_pair, has_fasta_gff, has_dir))
    if local_sources > 1:
        parser.error("options -d/--directory, --pair, and --fasta/--gff are mutually exclusive")
    if local_sources > 0 and has_ncbi_opts:
        parser.error("--cache-dir, --min-busco, and --workers cannot be used with local inputs (-d, --pair, --fasta/--gff)")

    if has_pair:
        raw_pairs = parse_explicit_pairs(args.pair)
        if len(raw_pairs) < args.num_genomes:
            parser.error(f"requested -n {args.num_genomes}, but only {len(raw_pairs)} --pair given")
        selected_pairs = raw_pairs[:args.num_genomes]
        return [(prefix_for(f).rstrip("_"), f, g) for f, g in selected_pairs]

    if has_fasta_gff:
        if not args.fasta or not args.gff or len(args.fasta) != len(args.gff):
            parser.error("--fasta and --gff must be supplied together with the same number of files")
        if len(args.fasta) < args.num_genomes:
            parser.error(f"requested -n {args.num_genomes}, but only {len(args.fasta)} files given")
        selected_pairs = list(zip(args.fasta[:args.num_genomes], args.gff[:args.num_genomes]))
        return [(prefix_for(f).rstrip("_"), f, g) for f, g in selected_pairs]

    if has_dir:
        fastas = source_files(args.directory, FASTA_SUFFIXES)
        gffs = source_files(args.directory, GFF_SUFFIXES)
        raw_pairs = pair_inputs(fastas, gffs)
        if len(raw_pairs) < args.num_genomes:
            parser.error(f"requested -n {args.num_genomes}, but only {len(raw_pairs)} pairs found in {args.directory}")
        selected_pairs = raw_pairs[:args.num_genomes]
        return [(prefix_for(f).rstrip("_"), f, g) for f, g in selected_pairs]

    min_busco = 0.95 if args.min_busco is None else args.min_busco
    workers = 10 if args.workers is None else args.workers
    if not 0.0 <= min_busco <= 1.0:
        parser.error("--min-busco must be between 0 and 1")
    if not 1 <= workers <= 30:
        parser.error("--workers must be between 1 and 30")

    cache = (args.cache_dir or Path(".ncbi-datasets")).resolve()
    cache.mkdir(parents=True, exist_ok=True)
    records = get_selection(cache, min_busco, args.num_genomes)
    pairs = download(cache, records, workers)
    return [
        (f"{rec['accession']} ({rec['organism']['organism_name']})", f, g)
        for rec, (f, g) in zip(records, pairs)
    ]


def run_anno(
    anno_bin: Path,
    fasta_path: Path,
    model_path: Path,
    threads: int,
    train_gff: Path | None = None,
    out_gff: Path | None = None,
    iterations: int | None = None,
) -> None:
    cmd = [str(anno_bin), "-m", str(model_path), "-p", str(threads)]
    if iterations is not None and iterations > 0:
        cmd.extend(["-e", str(iterations)])
    cmd.append(str(fasta_path))
    if train_gff is not None:
        cmd.append(str(train_gff))
    log("+ " + " ".join(cmd))
    if out_gff is not None:
        with out_gff.open("w", encoding="utf-8") as stdout_file:
            subprocess.run(cmd, check=True, stdout=stdout_file)
    else:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def assign_stratified_folds(
    windows_by_genome: dict[int, list[SampledWindow]],
    folds: int,
    rng: random.Random,
) -> list[list[SampledWindow]]:
    fold_buckets: list[list[SampledWindow]] = [[] for _ in range(folds)]
    for genome_idx in sorted(windows_by_genome):
        pool = list(windows_by_genome[genome_idx])
        rng.shuffle(pool)
        for idx, win in enumerate(pool):
            fold_buckets[idx % folds].append(win)
    for bucket in fold_buckets:
        bucket.sort(key=lambda w: (w.genome_index, w.seqid, w.start))
    return fold_buckets


def mean_std(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    m = sum(values) / len(values)
    if len(values) < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in values) / (len(values) - 1)
    return m, math.sqrt(var)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-n", "--num-genomes", type=int, required=True,
                        help="number of genomes to sample and train on (e.g. -n 2, -n 3)")
    parser.add_argument("-w", "--windows-per-genome", type=int, default=DEFAULT_WINDOWS_PER_GENOME,
                        help=f"number of CRF windows sampled per genome (default: {DEFAULT_WINDOWS_PER_GENOME})")
    parser.add_argument("-k", "--folds", type=int, default=DEFAULT_FOLDS,
                        help=f"number of cross-validation folds (default: {DEFAULT_FOLDS})")
    parser.add_argument("--cwt-pad", type=int, default=DEFAULT_CWT_PAD,
                        help=f"flanking sequence padding in bp for CWT boundary protection (default: {DEFAULT_CWT_PAD})")
    parser.add_argument("-p", "--threads", type=int, default=4,
                        help="number of worker threads for anno and blastp (default: 4)")
    parser.add_argument("-e", "--iterations", type=int, default=0,
                        help="maximum L-BFGS iterations for anno (default: anno default)")
    parser.add_argument("-m", "--model", type=Path,
                        help="output path for the final trained model (default: anno.model)")
    parser.add_argument("--anno", type=Path, default=Path("./anno"),
                        help="path to anno executable (default: ./anno)")
    parser.add_argument("--seed", type=int, default=42,
                        help="random seed for window sampling and fold assignment (default: 42)")
    parser.add_argument("--no-final-train", action="store_true",
                        help="skip training the final model on all sampled windows after CV")
    parser.add_argument("--workdir", type=Path,
                        help="optional directory to keep CV fold files and intermediate outputs")
    parser.add_argument("-d", "--directory", type=Path,
                        help="optional local directory containing FASTA and GFF3 files")
    parser.add_argument("--pair", action="append", default=[], metavar="FASTA:GFF3",
                        help="explicit local FASTA:GFF3 pair (repeatable)")
    parser.add_argument("--fasta", nargs="+", type=Path,
                        help="explicit local FASTA files, paired with --gff")
    parser.add_argument("--gff", nargs="+", type=Path,
                        help="explicit local GFF3 files, paired with --fasta")
    parser.add_argument("--cache-dir", type=Path,
                        help="NCBI FTP download cache directory (default: .ncbi-datasets)")
    parser.add_argument("--min-busco", type=float,
                        help="minimum annotation BUSCO complete fraction (default: 0.95)")
    parser.add_argument("--workers", type=int,
                        help="FTP download workers (1-30, default: 10)")
    args = parser.parse_args()

    if args.no_final_train and args.model is not None:
        parser.error("-m / --model cannot be used with --no-final-train")
    if args.model is None:
        args.model = Path("anno.model")

    if args.num_genomes < 1:
        parser.error("-n / --num-genomes must be >= 1")
    if args.windows_per_genome < 1:
        parser.error("-w / --windows-per-genome must be >= 1")
    if args.folds < 2:
        parser.error("-k / --folds must be >= 2")
    if args.cwt_pad < 0:
        parser.error("--cwt-pad must be >= 0")
    if args.threads < 1:
        parser.error("-p / --threads must be >= 1")

    anno_bin = args.anno.resolve()
    if not anno_bin.is_file():
        parser.error(f"anno executable not found: {anno_bin} (run 'make' first)")

    rng = random.Random(args.seed)
    genomes = resolve_genome_inputs(args, parser)
    log(f"Selected {len(genomes)} genome(s) for training:")

    genome_sequences: dict[int, dict[str, str]] = {}
    genome_gff_lines: dict[int, dict[str, list[str]]] = {}
    windows_by_genome: dict[int, list[SampledWindow]] = {}

    for idx, (label, fasta_path, gff_path) in enumerate(genomes, 1):
        log(f"  [{idx}/{len(genomes)}] Loading {label}: {fasta_path.name} + {gff_path.name}")
        seqs = {name: seq for name, seq in read_fasta(fasta_path)}
        spans, gff_by_seqid = compute_contig_gene_spans(seqs, gff_path)
        all_windows = partition_crf_windows(idx, label, seqs, spans, args.cwt_pad)
        sampled = select_windows_for_genome(all_windows, args.windows_per_genome, rng)
        if len(sampled) < args.folds:
            raise ValueError(
                f"Genome {label} yielded only {len(sampled)} CRF window(s), fewer than {args.folds} folds"
            )
        coding_win = sum(1 for w in sampled if w.gene_count > 0)
        total_genes = sum(w.gene_count for w in sampled)
        log(
            f"    Total CRF windows: {len(all_windows):,} -> sampled {len(sampled):,} "
            f"({coding_win:,} coding windows, {total_genes:,} non-overlapping valid genes, "
            f"CWT pad={args.cwt_pad} bp)"
        )
        genome_sequences[idx] = seqs
        genome_gff_lines[idx] = gff_by_seqid
        windows_by_genome[idx] = sampled

    fold_buckets = assign_stratified_folds(windows_by_genome, args.folds, rng)
    total_sampled = sum(len(b) for b in fold_buckets)
    log(f"Starting {args.folds}-fold cross-validation across {total_sampled:,} sampled CRF windows...")

    keep_workdir = args.workdir is not None
    if keep_workdir:
        workdir = args.workdir.resolve()
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        workdir = Path(tempfile.mkdtemp(prefix="anno-cv-", dir="."))

    can_evaluate = (
        shutil.which("gffread") is not None
        and shutil.which("makeblastdb") is not None
        and shutil.which("blastp") is not None
        and shutil.which("gffcompare") is not None
    )
    if not can_evaluate:
        log("WARNING: gffread/blastp/gffcompare not found in PATH; only prediction gene counts will be reported.")

    fold_metrics: list[dict[str, tuple[float, float, float]]] = []
    iterations = args.iterations if args.iterations > 0 else None

    try:
        for fold_idx in range(args.folds):
            fold_num = fold_idx + 1
            fold_dir = workdir / f"fold-{fold_num}"
            fold_dir.mkdir(parents=True, exist_ok=True)

            val_windows = fold_buckets[fold_idx]
            train_windows = [
                w for i, bucket in enumerate(fold_buckets) if i != fold_idx for w in bucket
            ]
            train_windows.sort(key=lambda w: (w.genome_index, w.seqid, w.start))

            train_fa = fold_dir / "train.fasta"
            train_gff = fold_dir / "train.gff3"
            val_fa = fold_dir / "val.fasta"
            val_gff = fold_dir / "val.gff3"
            pred_gff = fold_dir / "predict.gff3"
            fold_model = fold_dir / "fold.model"

            write_window_dataset(train_windows, genome_sequences, genome_gff_lines, train_fa, train_gff)
            write_window_dataset(val_windows, genome_sequences, genome_gff_lines, val_fa, val_gff)

            log(
                f"=== Fold {fold_num}/{args.folds}: "
                f"train={len(train_windows):,} windows, val={len(val_windows):,} windows ==="
            )
            run_anno(anno_bin, train_fa, fold_model, args.threads, train_gff=train_gff, iterations=iterations)
            run_anno(anno_bin, val_fa, fold_model, args.threads, out_gff=pred_gff)

            if can_evaluate:
                eval_dir = fold_dir / "eval"
                eval_dir.mkdir(exist_ok=True)
                metrics = evaluate_predictions(val_gff, pred_gff, val_fa, eval_dir, args.threads)
                fold_metrics.append(metrics)
                g_p, g_r, g_f1 = metrics["gene"]
                e_p, e_r, e_f1 = metrics["exon"]
                print(
                    f"Fold {fold_num}/{args.folds}: "
                    f"gene F1={g_f1:.4%} (P={g_p:.4%}, R={g_r:.4%}) | "
                    f"exon F1={e_f1:.4%} (P={e_p:.4%}, R={e_r:.4%})",
                    flush=True,
                )

        if fold_metrics:
            g_f1_m, g_f1_s = mean_std([m["gene"][2] for m in fold_metrics])
            g_p_m, g_p_s = mean_std([m["gene"][0] for m in fold_metrics])
            g_r_m, g_r_s = mean_std([m["gene"][1] for m in fold_metrics])
            e_f1_m, e_f1_s = mean_std([m["exon"][2] for m in fold_metrics])
            e_p_m, e_p_s = mean_std([m["exon"][0] for m in fold_metrics])
            e_r_m, e_r_s = mean_std([m["exon"][1] for m in fold_metrics])
            print(
                f"\n{args.folds}-Fold Cross-Validation Summary ({len(genomes)} genome(s), "
                f"{total_sampled:,} windows):"
            )
            print(
                f"  gene-level (protein):       F1={g_f1_m:.4%} ± {g_f1_s:.4%} "
                f"(precision={g_p_m:.4%} ± {g_p_s:.4%}, recall={g_r_m:.4%} ± {g_r_s:.4%})"
            )
            print(
                f"  gffcompare CDS exon-level:  F1={e_f1_m:.4%} ± {e_f1_s:.4%} "
                f"(precision={e_p_m:.4%} ± {e_p_s:.4%}, recall={e_r_m:.4%} ± {e_r_s:.4%})"
            )

        if not args.no_final_train:
            all_sampled = [
                w for idx in sorted(windows_by_genome) for w in windows_by_genome[idx]
            ]
            final_dir = workdir / "final"
            final_dir.mkdir(parents=True, exist_ok=True)
            final_fa = final_dir / "train.fasta"
            final_gff = final_dir / "train.gff3"
            write_window_dataset(all_sampled, genome_sequences, genome_gff_lines, final_fa, final_gff)
            log(f"=== Final training on all {len(all_sampled):,} sampled windows -> {args.model} ===")
            run_anno(
                anno_bin,
                final_fa,
                args.model.resolve(),
                args.threads,
                train_gff=final_gff,
                iterations=iterations,
            )
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())

