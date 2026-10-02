#!/usr/bin/env python3
"""Accuracy metrics for reference and query GFF3 files.
    micromamba run -n anno python3 accuracy.py \
        -r reference.gff3 -q anno.gff3 -g genome.fasta.gz
"""

import argparse
import gzip
import math
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def attributes(text):
    result = {}
    for field in text.split(";"):
        field = field.strip()
        if not field:
            continue
        key, separator, value = field.partition("=")
        if not separator:
            parts = field.split(None, 1)
            if len(parts) == 2:
                key, value = parts
            else:
                continue
        result[key] = value.strip('"')
    return result


def gff_records(path):
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.rstrip() == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected GFF3/GTF with 9 columns")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 1 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            if fields[2].lower() in ("exon", "cds") and fields[6] not in ("+", "-"):
                raise ValueError(f"{path}:{line_number}: invalid strand")
            yield fields, attributes(fields[8]), start, end


def feature_parents(info):
    return [parent for parent in (info.get("Parent") or info.get("transcript_id") or info.get("gene_id") or "").split(",")
            if parent]


def transcript_intervals(path, feature):
    metadata, intervals, plain = {}, defaultdict(list), set()
    for fields, info, start, end in gff_records(path):
        identifier = info.get("ID") or info.get("transcript_id")
        if fields[2].lower() in ("mrna", "transcript") and identifier:
            gene = (info.get("Parent") or info.get("gene_id") or identifier).split(",")[0]
            metadata[identifier] = (gene or identifier, fields[0], fields[6])
        elif fields[2].lower() == feature:
            interval = (fields[0], fields[6], start, end)
            parents = feature_parents(info)
            if not parents:
                plain.add(interval)
            for parent in parents:
                intervals[parent].append(interval)
                if info.get("gene_id"):
                    metadata.setdefault(parent, (info["gene_id"], fields[0], fields[6]))
    models = {}
    for transcript, segments in intervals.items():
        seqid, strand = segments[0][:2]
        gene, seqid, strand = metadata.get(transcript, (transcript, seqid, strand))
        if all(segment[:2] == (seqid, strand) for segment in segments):
            models[transcript] = gene, seqid, strand, segments
    return models, plain


def read_exons(path):
    transcripts, plain = transcript_intervals(path, "exon")
    selected = {}
    for gene, seqid, strand, intervals in transcripts.values():
        length = sum(end - start + 1 for _, _, start, end in intervals)
        key = seqid, strand, gene
        if key not in selected or length > selected[key][0]:
            selected[key] = length, intervals
    return plain | {interval for _, intervals in selected.values() for interval in intervals}


def read_genes(path):
    """Map each coding gene to its longest-CDS transcript and CDS-spanning region.

    Returns gene_id -> (seqid, start, end, strand, transcript_id).
    """
    transcripts, _ = transcript_intervals(path, "cds")
    genes = {}
    for transcript, (gene, seqid, strand, intervals) in transcripts.items():
        length = sum(end - start + 1 for _, _, start, end in intervals)
        start = min(interval[2] for interval in intervals)
        end = max(interval[3] for interval in intervals)
        previous = genes.get(gene)
        if previous is None or length > previous[0]:
            genes[gene] = (length, seqid, start, end, strand, transcript)
    return {gene: value[1:] for gene, value in genes.items()}


def read_genome(path):
    sequences = {name: sequence.upper() for name, sequence in read_fasta(path)}
    if not sequences:
        raise ValueError(f"{path}: no FASTA sequences found")
    return sequences


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def coding_transcripts(path, genome):
    """Read coding transcripts and apply the paper's invalid-model filters."""
    transcripts, _ = transcript_intervals(path, "cds")
    valid = {}
    for transcript, (gene, seqid, strand, intervals) in transcripts.items():
        if seqid not in genome:
            raise ValueError(f"{path}: sequence {seqid} is missing from the genome FASTA")
        ordered = sorted(intervals, key=lambda item: item[2])
        if any(end > len(genome[seqid]) for _, _, _, end in ordered):
            raise ValueError(f"{path}: CDS extends beyond sequence {seqid}")
        if any(next_start - end - 1 <= 1
               for (_, _, _, end), (_, _, next_start, _) in zip(ordered, ordered[1:])):
            continue
        sequence = "".join(
            genome[seqid][start - 1:end] for _, _, start, end in ordered
        )
        if strand == "-":
            sequence = reverse_complement(sequence)
        if (len(sequence) < 6 or len(sequence) % 3 or sequence[:3] != "ATG" or
                sequence[-3:] not in {"TAA", "TAG", "TGA"} or
                any(sequence[index:index + 3] in {"TAA", "TAG", "TGA"}
                    for index in range(3, len(sequence) - 3, 3))):
            continue
        coding_length = sum(end - start + 1 for _, _, start, end in intervals)
        valid[transcript] = {
            "gene": gene,
            "seqid": seqid,
            "strand": strand,
            "intervals": [(start, end) for _, _, start, end in intervals],
            "length": coding_length,
        }
    return valid


def filtered_gff(path, genome, output, longest_only=False):
    transcripts = coding_transcripts(path, genome)
    if longest_only:
        selected = {}
        for transcript, record in transcripts.items():
            key = (record["seqid"], record["strand"], record["gene"])
            previous = selected.get(key)
            if previous is None or (record["length"], transcript) > (previous[1]["length"], previous[0]):
                selected[key] = (transcript, record)
        transcripts = {transcript: record for transcript, record in selected.values()}

    genes = {}
    for record in transcripts.values():
        key = record["seqid"], record["strand"], record["gene"]
        start = min(interval[0] for interval in record["intervals"])
        end = max(interval[1] for interval in record["intervals"])
        previous = genes.get(key, (start, end))
        genes[key] = min(start, previous[0]), max(end, previous[1])
    written = set()
    with open(output, "w") as stream:
        stream.write("##gff-version 3\n")
        for transcript, record in sorted(transcripts.items()):
            start = min(interval[0] for interval in record["intervals"])
            end = max(interval[1] for interval in record["intervals"])
            gene = record["gene"]
            gene_key = (record["seqid"], record["strand"], gene)
            if gene_key not in written:
                written.add(gene_key)
                gene_start, gene_end = genes[gene_key]
                stream.write(
                    f"{record['seqid']}\taccuracy\tgene\t{gene_start}\t{gene_end}\t.\t"
                    f"{record['strand']}\t.\tID={gene}\n"
                )
            stream.write(
                f"{record['seqid']}\taccuracy\ttranscript\t{start}\t{end}\t.\t"
                f"{record['strand']}\t.\tID={transcript};Parent={gene}\n"
            )
            for number, (exon_start, exon_end) in enumerate(sorted(record["intervals"]), 1):
                stream.write(
                    f"{record['seqid']}\taccuracy\texon\t{exon_start}\t{exon_end}\t.\t"
                    f"{record['strand']}\t.\tID={transcript}.exon{number};Parent={transcript}\n"
                )
    return len(transcripts)


def parse_gffcompare_stats(path):
    metrics = {}
    pattern = re.compile(
        r"^\s*(Base|Exon|Intron|Intron chain|Transcript|Locus) level:\s*"
        r"([0-9.]+)\s*\|\s*([0-9.]+)"
    )
    with open(path) as stream:
        for line in stream:
            match = pattern.match(line)
            if not match:
                continue
            name = match.group(1).lower().replace(" ", "_")
            recall, precision = float(match.group(2)) / 100, float(match.group(3)) / 100
            tp = precision * recall
            metrics[name] = scores(
                tp,
                recall - tp,
                precision - tp,
            )
    if not metrics:
        raise ValueError(f"could not parse gffcompare statistics: {path}")
    return metrics


def run_gffcompare(reference, query, workdir, label):
    prefix = workdir / label
    subprocess.run(
        ["gffcompare", "--no-exon-merge", "--strict-match", "-r", str(reference),
         "-o", str(prefix), str(query)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    stats_path = pathlib.Path(f"{prefix}.stats")
    if not stats_path.exists() and prefix.is_file():
        stats_path = prefix
    return parse_gffcompare_stats(stats_path)


def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def overlap(left, right):
    total = i = j = 0
    while i < len(left) and j < len(right):
        start = max(left[i][0], right[j][0])
        end = min(left[i][1], right[j][1])
        total += max(0, end - start + 1)
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def scores(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def evaluate(reference, query):
    exact_tp = len(reference & query)
    exact_fp = len(query - reference)
    exact_fn = len(reference - query)
    ref_groups = defaultdict(list)
    query_groups = defaultdict(list)
    for seqid, strand, start, end in reference:
        ref_groups[(seqid, strand)].append((start, end))
    for seqid, strand, start, end in query:
        query_groups[(seqid, strand)].append((start, end))
    bp_tp = bp_fp = bp_fn = 0
    for group in set(ref_groups) | set(query_groups):
        ref = merge(ref_groups[group])
        pred = merge(query_groups[group])
        shared = overlap(ref, pred)
        ref_bp = sum(end - start + 1 for start, end in ref)
        query_bp = sum(end - start + 1 for start, end in pred)
        bp_tp += shared
        bp_fp += query_bp - shared
        bp_fn += ref_bp - shared
    return scores(bp_tp, bp_fp, bp_fn), scores(exact_tp, exact_fp, exact_fn), (bp_tp, bp_fp, bp_fn), (exact_tp, exact_fp, exact_fn)


def genome_fasta(genome, workdir):
    """Return a plain, gffread-indexable FASTA path (decompressing .gz if needed)."""
    if str(genome).endswith(".gz"):
        with gzip.GzipFile(filename=genome, mode="rb") as source, \
                tempfile.NamedTemporaryFile(mode="wb", prefix="genome-", suffix=".fa", dir=workdir, delete=False) as target:
            shutil.copyfileobj(source, target)
            return pathlib.Path(target.name)
    link = workdir / "genome.fa"
    source = pathlib.Path(genome).resolve()
    if pathlib.Path(genome).absolute() == link.absolute():
        return link
    if link.is_symlink():
        (workdir / "genome.fa.fai").unlink(missing_ok=True)
    if link.resolve() == source:
        return link
    if link.is_symlink():
        link.unlink()
    if link.exists():
        raise ValueError(f"Refusing to overwrite existing FASTA: {link}")
    link.symlink_to(source)
    return link


def write_protein_gff(path, genes, output):
    selected_genes = set(genes)
    selected_transcripts = {value[4] for value in genes.values()}
    with open(output, "w") as target:
        target.write("##gff-version 3\n")
        for fields, info, start, end in gff_records(path):
            feature = fields[2].lower()
            identifier = info.get("ID") or info.get("transcript_id")
            parents = [parent for parent in feature_parents(info) if parent in selected_transcripts]
            keep = ((feature == "gene" and identifier in selected_genes) or
                    (feature in ("mrna", "transcript") and identifier in selected_transcripts) or
                    (feature == "cds" and bool(parents)))
            if not keep:
                continue
            if fields[6] not in ("+", "-"):
                raise ValueError(f"{path}: selected coding feature has invalid strand")
            if feature == "cds" and "Parent" in info:
                info["Parent"] = ",".join(parents)
                fields[8] = ";".join(f"{key}={value}" for key, value in info.items())
            target.write("\t".join(fields) + "\n")


def read_fasta(path):
    name, chunks = None, []
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].split()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name, chunks = header[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"{path}:{line_number}: FASTA sequence before header")
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def extract_proteins(gff, genome, genes, out_fasta, workdir, tag):
    """Translate the longest-CDS transcript of each gene, header renamed to gene id."""
    if not genes:
        out_fasta.write_text("")
        return 0
    raw = workdir / f"{tag}.raw.faa"
    coding_gff = workdir / f"{tag}.coding.gff3"
    write_protein_gff(gff, genes, coding_gff)
    try:
        subprocess.run(["gffread", str(coding_gff), "-g", str(genome), "-y", str(raw)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as error:
        detail = error.stderr.strip().splitlines()[-1] if error.stderr else "unknown error"
        raise ValueError(f"gffread failed for {gff}: {detail}") from error
    wanted = {transcript: gene for gene, (_, _, _, _, transcript) in genes.items()}
    written = 0
    with open(out_fasta, "w") as target:
        for name, sequence in read_fasta(raw):
            gene = wanted.get(name)
            if gene is None:
                continue
            sequence = sequence.replace(".", "").replace("*", "")
            if not sequence:
                continue
            target.write(f">{gene}\n{sequence}\n")
            written += 1
    return written


def blast_hits(query_fasta, ref_fasta, workdir, evalue, threads):
    database = workdir / "refdb"
    subprocess.run(["makeblastdb", "-in", str(ref_fasta), "-dbtype", "prot", "-out", str(database)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    columns = "6 qseqid sseqid bitscore evalue qstart qend qlen sstart send slen"
    result = subprocess.run(
        ["blastp", "-query", str(query_fasta), "-db", str(database), "-evalue", str(evalue),
         "-outfmt", columns, "-max_target_seqs", "5", "-num_threads", str(threads)],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    hits = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 10:
            continue
        query, subject = parts[0], parts[1]
        hit_score, hit_evalue = float(parts[2]), float(parts[3])
        qstart, qend, qlen = int(parts[4]), int(parts[5]), int(parts[6])
        sstart, send, slen = int(parts[7]), int(parts[8]), int(parts[9])
        qcov = (qend - qstart + 1) / qlen if qlen else 0.0
        scov = (send - sstart + 1) / slen if slen else 0.0
        hits.append((query, subject, hit_score, hit_evalue, qcov, scov))
    return hits


def region_overlap(a, b):
    return a[0] == b[0] and a[1] <= b[2] and b[1] <= a[2]


def gene_level(reference_gff, query_gff, genome, workdir, *, bitscore, evalue, coverage, threads):
    ref_genes = read_genes(reference_gff)
    query_genes = read_genes(query_gff)
    indexed = genome_fasta(genome, workdir)
    ref_fasta, query_fasta = workdir / "reference.faa", workdir / "query.faa"
    positives = extract_proteins(reference_gff, indexed, ref_genes, ref_fasta, workdir, "reference")
    predicted = extract_proteins(query_gff, indexed, query_genes, query_fasta, workdir, "query")
    if positives == 0 or predicted == 0:
        return scores(0, predicted, positives), (0, positives, predicted)

    ref_region = {gene: (seqid, start, end) for gene, (seqid, start, end, *_) in ref_genes.items()}
    query_region = {gene: (seqid, start, end) for gene, (seqid, start, end, *_) in query_genes.items()}

    candidates = []
    for query, subject, hit_score, hit_evalue, qcov, scov in blast_hits(query_fasta, ref_fasta, workdir, evalue, threads):
        if hit_score < bitscore or hit_evalue > evalue:
            continue
        if qcov <= coverage or scov <= coverage:
            continue
        region_q, region_r = query_region.get(query), ref_region.get(subject)
        if (region_q is None or region_r is None or not region_overlap(region_q, region_r)
            or query_genes[query][3] != ref_genes[subject][3]):
            continue
        candidates.append((hit_score, query, subject))

    candidates.sort(reverse=True)
    used_query, used_ref, tp = set(), set(), 0
    for _, query, subject in candidates:
        if query in used_query or subject in used_ref:
            continue
        used_query.add(query)
        used_ref.add(subject)
        tp += 1
    return scores(tp, predicted - tp, positives - tp), (tp, positives, predicted)


def main():
    parser = argparse.ArgumentParser(description="Overall exon, bp and protein-level gene accuracy")
    parser.add_argument("-r", "--reference", required=True, help="reference GFF3/GTF(.gz)")
    parser.add_argument("-q", "--query", required=True, help="query GFF3/GTF(.gz)")
    parser.add_argument("-g", "--genome", help="genome FASTA(.gz); enables protein-level gene F1")
    parser.add_argument("--bitscore", type=float, default=50.0, help="minimum blastp bit score")
    parser.add_argument("--evalue", type=float, default=1e-5, help="maximum blastp e-value")
    parser.add_argument("--coverage", type=float, default=0.70,
                        help="minimum best-HSP coverage of both proteins (fraction)")
    parser.add_argument("--threads", type=int, default=4, help="blastp threads")
    parser.add_argument("--workdir", type=pathlib.Path, help="keep protein/blast intermediates here")
    parser.add_argument("--gffcompare", action="store_true",
                        help="run CDS-only gffcompare exon/locus evaluation")
    args = parser.parse_args()
    if args.gffcompare and not args.genome:
        parser.error("--gffcompare requires --genome")
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not math.isfinite(args.bitscore) or args.bitscore < 0:
        parser.error("--bitscore must be finite and nonnegative")
    if not math.isfinite(args.evalue) or args.evalue <= 0:
        parser.error("--evalue must be finite and positive")
    if not 0 <= args.coverage < 1:
        parser.error("--coverage must be in [0, 1)")

    try:
        reference = read_exons(args.reference)
        query = read_exons(args.query)
        bp, exon, bp_counts, exon_counts = evaluate(reference, query)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"bp-level: F1={bp[2]:.4%} precision={bp[0]:.4%} recall={bp[1]:.4%}")
    print(f"  bases TP={bp_counts[0]:,} FP={bp_counts[1]:,} FN={bp_counts[2]:,}")
    print(f"exon-level: F1={exon[2]:.4%} precision={exon[0]:.4%} recall={exon[1]:.4%}")
    print(f"  exons TP={exon_counts[0]:,} FP={exon_counts[1]:,} FN={exon_counts[2]:,}")

    if not args.genome:
        return 0

    keep = args.workdir is not None
    if keep:
        workdir = args.workdir
    else:
        base = "tmp" if pathlib.Path("tmp").is_dir() else "."
        workdir = pathlib.Path(tempfile.mkdtemp(prefix="anno-gene-", dir=base))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        gene, counts = gene_level(args.reference, args.query, args.genome, workdir,
                                  bitscore=args.bitscore, evalue=args.evalue,
                                  coverage=args.coverage, threads=args.threads)
        if args.gffcompare:
            genome = read_genome(args.genome)
            reference_all = workdir / "reference.cds.gff3"
            reference_longest = workdir / "reference.longest.cds.gff3"
            query_cds = workdir / "query.cds.gff3"
            reference_count = filtered_gff(args.reference, genome, reference_all)
            filtered_gff(args.reference, genome, reference_longest, longest_only=True)
            query_count = filtered_gff(args.query, genome, query_cds)
            exon_gff = locus_gff = scores(0, query_count, reference_count)
            if reference_count and query_count:
                exon_gff = run_gffcompare(reference_longest, query_cds, workdir, "gffcompare_exon")["exon"]
                locus_gff = run_gffcompare(reference_all, query_cds, workdir, "gffcompare_locus")["locus"]
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)
    print(f"gene-level (protein): F1={gene[2]:.4%} precision={gene[0]:.4%} recall={gene[1]:.4%}")
    print(f"  genes TP={counts[0]:,} P(reference)={counts[1]:,} PP(predicted)={counts[2]:,}")
    if args.gffcompare:
        print(f"gffcompare CDS exon-level: F1={exon_gff[2]:.4%} precision={exon_gff[0]:.4%} recall={exon_gff[1]:.4%}")
        print(f"gffcompare CDS locus-level: F1={locus_gff[2]:.4%} precision={locus_gff[0]:.4%} recall={locus_gff[1]:.4%}")
        print(f"  filtered transcripts: reference={reference_count:,} query={query_count:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
