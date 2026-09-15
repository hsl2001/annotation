#!/usr/bin/env python3
"""Exact exon and covered-base accuracy without gffcompare."""

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from typing import Any


def open_text(path):
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def read_exons(path):
    exons = set()
    duplicates = 0
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) in (3, 6):
                seqid, start, end = fields[:3]
                strand = fields[5] if len(fields) == 6 else "."
                if strand == ".":
                    strand = "+"
                bed = True
            elif len(fields) == 9:
                if fields[2].lower() != "exon":
                    continue
                seqid, start, end, strand = fields[0], fields[3], fields[4], fields[6]
                bed = False
            else:
                raise ValueError(f"{path}:{line_number}: expected GFF3 or BED fields")
            try:
                start, end = int(start), int(end)
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid exon coordinates") from error
            if bed:
                if start < 0 or end <= start:
                    raise ValueError(f"{path}:{line_number}: invalid BED interval")
                start += 1
            if start < 1 or end < start or strand not in ("+", "-"):
                raise ValueError(f"{path}:{line_number}: invalid exon interval or strand")
            exon = (seqid, strand, start, end)
            duplicates += exon in exons
            exons.add(exon)
    return exons, duplicates


def merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def union_by_group(exons):
    groups = defaultdict(list)
    for seqid, strand, start, end in exons:
        groups[(seqid, strand)].append((start, end))
    return {key: merge(intervals) for key, intervals in groups.items()}


def interval_length(intervals):
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


def metrics(reference, prediction) -> dict[str, Any]:
    true_positive = len(reference & prediction)
    false_negative = len(reference - prediction)
    false_positive = len(prediction - reference)
    ref_union, pred_union = union_by_group(reference), union_by_group(prediction)
    groups = set(ref_union) | set(pred_union)
    ref_bases = sum(interval_length(ref_union.get(group, [])) for group in groups)
    pred_bases = sum(interval_length(pred_union.get(group, [])) for group in groups)
    base_true_positive = sum(overlap_length(ref_union.get(group, []), pred_union.get(group, [])) for group in groups)
    return {
        "reference_exons": len(reference),
        "prediction_exons": len(prediction),
        "true_positive": true_positive,
        "false_positive": false_positive,
        "false_negative": false_negative,
        "reference_bases": ref_bases,
        "prediction_bases": pred_bases,
        "base_true_positive": base_true_positive,
        "base_false_positive": pred_bases - base_true_positive,
        "base_false_negative": ref_bases - base_true_positive,
        "groups": len(groups),
    }


def scores(values: dict[str, Any], prefix="") -> dict[str, float]:
    tp = values[prefix + "true_positive"]
    fp = values[prefix + "false_positive"]
    fn = values[prefix + "false_negative"]
    precision = tp / (tp + fp) if tp + fp else 0.0
    sensitivity = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else 0.0
    return {"precision": precision, "sensitivity": sensitivity, "f1": f1}


def report(reference, prediction) -> dict[str, Any]:
    result = metrics(reference, prediction)
    result["exon"] = scores(result)
    result["base"] = scores(result, "base_")
    by_strand: dict[str, dict[str, Any]] = {}
    for strand in ("+", "-"):
        ref = {exon for exon in reference if exon[1] == strand}
        pred = {exon for exon in prediction if exon[1] == strand}
        strand_result = metrics(ref, pred)
        strand_result["exon"] = scores(strand_result)
        strand_result["base"] = scores(strand_result, "base_")
        by_strand[strand] = strand_result
    result["strand"] = by_strand
    return result


def print_metric(label, result):
    exon, base = result["exon"], result["base"]
    print(f"{label}: exon Sn={exon['sensitivity']:.4%} Pr={exon['precision']:.4%} F1={exon['f1']:.4%}; "
          f"base Sn={base['sensitivity']:.4%} Pr={base['precision']:.4%} F1={base['f1']:.4%}")
    print(f"  exons ref={result['reference_exons']:,} pred={result['prediction_exons']:,} "
          f"TP={result['true_positive']:,} FP={result['false_positive']:,} FN={result['false_negative']:,}")
    print(f"  bases ref={result['reference_bases']:,} pred={result['prediction_bases']:,} "
          f"TP={result['base_true_positive']:,} FP={result['base_false_positive']:,} FN={result['base_false_negative']:,}")


def main():
    parser = argparse.ArgumentParser(description="Exact exon and covered-base accuracy without gffcompare")
    parser.add_argument("reference", help="reference GFF3/GTF, optionally .gz")
    parser.add_argument("prediction", help="prediction GFF3/GTF, optionally .gz")
    args = parser.parse_args()
    try:
        reference, reference_duplicates = read_exons(args.reference)
        prediction, prediction_duplicates = read_exons(args.prediction)
        result = report(reference, prediction)
        result["reference_duplicates"] = reference_duplicates
        result["prediction_duplicates"] = prediction_duplicates
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print_metric("overall", result)
    for strand in ("+", "-"):
        print_metric(f"strand {strand}", result["strand"][strand])
    print(f"deduplicated: reference={reference_duplicates:,} prediction={prediction_duplicates:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
