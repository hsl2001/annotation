#!/usr/bin/env python3
"""Accuracy metrics for reference and query GFF3 files.

Always reports nucleotide (bp) and exon-level F1. When a genome FASTA is supplied
(--genome), it additionally reports a protein-level gene F1 that mirrors the
BUSCO protein-mode protocol used to compare annotation tools such as anno and
AnnEvo: proteins are extracted with gffread, aligned with blastp, and a
predicted gene counts as a true positive only when its best hit against the
reference clears bit-score/e-value thresholds, sits in the same genomic region,
and covers more than a set fraction of both proteins. Run inside the micromamba
``anno`` environment so gffread/makeblastdb/blastp are on PATH:

    micromamba run -n anno python3 test_accuracy.py \
        -r reference.gff3 -q anno.gff3 -g genome.fasta.gz
    micromamba run -n anno python3 test_accuracy.py \
        -r reference.gff3 -q annevo.gff3 -g genome.fasta.gz
"""

import argparse
import gzip
import pathlib
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


def read_exons(path):
    transcripts = {}
    exon_groups = defaultdict(list)
    plain = defaultdict(list)
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected GFF3/GTF with 9 columns")
            feature = fields[2].lower()
            info = attributes(fields[8])
            identifier = info.get("ID") or info.get("transcript_id")
            parents = (info.get("Parent") or info.get("gene_id") or "").split(",")
            if feature in ("mrna", "transcript") and identifier:
                transcripts[identifier] = (parents[0] or identifier, fields[6])
            if feature != "exon":
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 1 or end < start or fields[6] not in ("+", "-"):
                raise ValueError(f"{path}:{line_number}: invalid exon interval or strand")
            parent = parents[0]
            interval = (fields[0], fields[6], start, end)
            if parent:
                exon_groups[parent].append(interval)
            else:
                plain[(fields[0], fields[6])].append(interval)

    selected = [interval for intervals in plain.values() for interval in intervals]
    by_gene = defaultdict(list)
    for transcript, intervals in exon_groups.items():
        locus, strand = transcripts.get(transcript, (transcript, intervals[0][1]))
        length = sum(end - start + 1 for _, _, start, end in intervals)
        by_gene[(locus, strand)].append((length, intervals))
    for candidates in by_gene.values():
        selected.extend(max(candidates, key=lambda item: item[0])[1])
    return set(selected)


def read_genes(path):
    """Map each coding gene to its longest-CDS transcript and CDS-spanning region.

    Returns gene_id -> (seqid, start, end, strand, transcript_id).
    """
    transcripts = {}
    cds = defaultdict(list)
    cds_meta = {}
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected GFF3/GTF with 9 columns")
            feature = fields[2].lower()
            info = attributes(fields[8])
            if feature in ("mrna", "transcript"):
                identifier = info.get("ID") or info.get("transcript_id")
                if identifier:
                    parent = (info.get("Parent") or info.get("gene_id") or identifier).split(",")[0]
                    transcripts[identifier] = (parent or identifier, fields[6], fields[0])
            elif feature == "cds":
                parent = (info.get("Parent") or info.get("transcript_id") or "").split(",")[0]
                if not parent:
                    continue
                cds[parent].append((int(fields[3]), int(fields[4])))
                cds_meta[parent] = (fields[0], fields[6])

    genes = {}
    for transcript, intervals in cds.items():
        gene, strand, seqid = transcripts.get(transcript, (None, None, None))
        if gene is None:
            seqid, strand = cds_meta[transcript]
            gene = transcript
        length = sum(end - start + 1 for start, end in intervals)
        start = min(start for start, _ in intervals)
        end = max(end for _, end in intervals)
        previous = genes.get(gene)
        if previous is None or length > previous[0]:
            genes[gene] = (length, seqid, start, end, strand, transcript)
    return {gene: value[1:] for gene, value in genes.items()}


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
        plain = workdir / "genome.fa"
        with gzip.open(genome, "rb") as source, open(plain, "wb") as target:
            shutil.copyfileobj(source, target)
        return plain
    link = workdir / "genome.fa"
    link.symlink_to(pathlib.Path(genome).resolve())
    return link


def read_fasta(path):
    name, chunks = None, []
    with open(path) as stream:
        for line in stream:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def extract_proteins(gff, genome, genes, out_fasta, workdir, tag):
    """Translate the longest-CDS transcript of each gene, header renamed to gene id."""
    raw = workdir / f"{tag}.raw.faa"
    subprocess.run(["gffread", str(gff), "-g", str(genome), "-y", str(raw)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
        raise ValueError("no protein sequences extracted; check genome and GFF3 coordinates")

    ref_region = {gene: (seqid, start, end) for gene, (seqid, start, end, *_) in ref_genes.items()}
    query_region = {gene: (seqid, start, end) for gene, (seqid, start, end, *_) in query_genes.items()}

    candidates = []
    for query, subject, hit_score, hit_evalue, qcov, scov in blast_hits(query_fasta, ref_fasta, workdir, evalue, threads):
        if hit_score < bitscore or hit_evalue > evalue:
            continue
        if qcov <= coverage or scov <= coverage:
            continue
        region_q, region_r = query_region.get(query), ref_region.get(subject)
        if region_q is None or region_r is None or not region_overlap(region_q, region_r):
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
    args = parser.parse_args()

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
    except (OSError, ValueError, subprocess.CalledProcessError, FileNotFoundError) as error:
        parser.error(str(error))
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)
    print(f"gene-level (protein): F1={gene[2]:.4%} precision={gene[0]:.4%} recall={gene[1]:.4%}")
    print(f"  genes TP={counts[0]:,} P(reference)={counts[1]:,} PP(predicted)={counts[2]:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
