#!/usr/bin/env python3
import argparse
import gzip
import subprocess
import sys
from collections import defaultdict

DEFAULT_GENOME = "data/Col-CC_v2_genome.fasta.gz"
DEFAULT_GFF = "data/TAIR12_1Feb26_HEADERED.gff3"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def read_intervals(path, kind):
    """Return {seqid: [(start, end), ...]} 1-based inclusive intervals."""
    by_seq = defaultdict(list)
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if kind == "gff":
                if len(fields) != 9 or fields[2].lower() != "exon":
                    continue
                seqid, start, end = fields[0], int(fields[3]), int(fields[4])
            else:
                if len(fields) < 3:
                    raise ValueError(f"{path}:{line_number}: expected BED fields")
                seqid, start, end = fields[0], int(fields[1]) + 1, int(fields[2])
            if end < start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            by_seq[seqid].append((start, end))
    return by_seq


def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def total_length(intervals):
    return sum(end - start + 1 for start, end in intervals)


def overlap_length(left, right):
    total = left_index = right_index = 0
    while left_index < len(left) and right_index < len(right):
        start = max(left[left_index][0], right[right_index][0])
        end = min(left[left_index][1], right[right_index][1])
        if start <= end:
            total += end - start + 1
        if left[left_index][1] < right[right_index][1]:
            left_index += 1
        else:
            right_index += 1
    return total


def base_scores(reference, prediction):
    true_positive = false_positive = false_negative = 0
    for seqid in set(reference) | set(prediction):
        ref = merge(reference.get(seqid, []))
        pred = merge(prediction.get(seqid, []))
        shared = overlap_length(ref, pred)
        true_positive += shared
        false_positive += total_length(pred) - shared
        false_negative += total_length(ref) - shared
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else 0.0
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1, true_positive, false_positive, false_negative


def main():
    parser = argparse.ArgumentParser(description="End-to-end bp-level exon F1 for anno")
    parser.add_argument("--genome", default=DEFAULT_GENOME, help="genome FASTA[.gz]")
    parser.add_argument("--gff", default=DEFAULT_GFF, help="reference GFF3 exon annotation")
    parser.add_argument("--bed", default="candidates.bed", help="prediction BED path")
    parser.add_argument("--binary", default="./anno", help="anno executable")
    parser.add_argument("--skip-run", action="store_true", help="reuse existing --bed")
    args = parser.parse_args()

    if not args.skip_run:
        with open(args.bed, "w") as bed:
            subprocess.run([args.binary, args.genome], stdout=bed, check=True)

    try:
        reference = read_intervals(args.gff, "gff")
        prediction = read_intervals(args.bed, "bed")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    precision, recall, f1, tp, fp, fn = base_scores(reference, prediction)
    print(f"bp-level: F1={f1:.4%} precision={precision:.4%} recall={recall:.4%}")
    print(f"  bases TP={tp:,} FP={fp:,} FN={fn:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
