#!/usr/bin/env python3
"""Overall exon and covered-base accuracy for reference and query GFF3 files."""

import argparse
import gzip
import sys
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


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


def main():
    parser = argparse.ArgumentParser(description="Overall exon and bp-level accuracy")
    parser.add_argument("-r", "--reference", required=True, help="reference GFF3/GTF(.gz)")
    parser.add_argument("-q", "--query", required=True, help="query GFF3/GTF(.gz)")
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
