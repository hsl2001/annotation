#!/usr/bin/env python3
"""Held-out accuracy for anno: train on every chromosome except one, predict the held-out one.

Reports strand-aware CDS base-level precision/recall/F1 (the standard coding-base metric)
and exact CDS-segment F1 against all reference isoforms.
"""

import argparse
import gzip
import pathlib
import subprocess
import sys
import tempfile
from collections import defaultdict

DEFAULT_GENOME = "data/Col-CC_v2_genome.fasta.gz"
DEFAULT_GFF = "data/TAIR12_1Feb26_HEADERED.gff3"


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def split_fasta(genome, holdout, train_path, test_path):
    with open_text(genome) as source, open(train_path, "w") as train, open(test_path, "w") as test:
        target = None
        for line in source:
            if line.startswith(">"):
                target = test if line[1:].split()[0] == holdout else train
            if target is None:
                raise ValueError("FASTA does not start with a header")
            target.write(line)


def split_gff(gff, holdout, train_path):
    with open_text(gff) as source, open(train_path, "w") as train:
        for line in source:
            if line.startswith("#") or line.split("\t", 1)[0] != holdout:
                train.write(line)


def read_cds(path, seqid):
    intervals = defaultdict(list)
    with open_text(path) as stream:
        for line in stream:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) == 9 and fields[2] == "CDS" and fields[0] == seqid:
                intervals[fields[6]].append((int(fields[3]), int(fields[4])))
    return intervals


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
        start, end = max(left[i][0], right[j][0]), min(left[i][1], right[j][1])
        total += max(0, end - start + 1)
        if left[i][1] < right[j][1]:
            i += 1
        else:
            j += 1
    return total


def f1(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return 2 * precision * recall / (precision + recall) if precision + recall else 0.0, precision, recall


def evaluate(reference, prediction):
    tp = fp = fn = 0
    for strand in "+-":
        ref, pred = merge(reference[strand]), merge(prediction[strand])
        shared = overlap(ref, pred)
        tp += shared
        fp += sum(e - s + 1 for s, e in pred) - shared
        fn += sum(e - s + 1 for s, e in ref) - shared
    base = f1(tp, fp, fn)
    ref_set = {(s, *i) for s in "+-" for i in reference[s]}
    pred_set = {(s, *i) for s in "+-" for i in prediction[s]}
    exact = f1(len(ref_set & pred_set), len(pred_set - ref_set), len(ref_set - pred_set))
    return base, exact


def main():
    parser = argparse.ArgumentParser(description="Held-out CDS accuracy for anno")
    parser.add_argument("--genome", default=DEFAULT_GENOME)
    parser.add_argument("--gff", default=DEFAULT_GFF)
    parser.add_argument("--holdout", default="Chr5", help="chromosome excluded from training")
    parser.add_argument("--binary", default="./anno")
    parser.add_argument("--epochs", type=int, default=6)
    parser.add_argument("--workdir", type=pathlib.Path, help="keep intermediate files here")
    args = parser.parse_args()

    workdir = args.workdir or pathlib.Path(tempfile.mkdtemp(prefix="anno-eval-"))
    workdir.mkdir(parents=True, exist_ok=True)
    train_fa, test_fa, train_gff = workdir / "train.fa", workdir / "test.fa", workdir / "train.gff3"
    model, prediction = workdir / "anno.model", workdir / "prediction.gff3"
    try:
        split_fasta(args.genome, args.holdout, train_fa, test_fa)
        split_gff(args.gff, args.holdout, train_gff)
        with open(prediction, "w") as out:
            subprocess.run([args.binary, "-m", model, "-e", str(args.epochs), train_fa, train_gff],
                           stdout=subprocess.DEVNULL, check=True)
            subprocess.run([args.binary, "-m", model, test_fa], stdout=out, check=True)
        base, exact = evaluate(read_cds(args.gff, args.holdout), read_cds(prediction, args.holdout))
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    print(f"{args.holdout} CDS bp-level: F1={base[0]:.4%} precision={base[1]:.4%} recall={base[2]:.4%}")
    print(f"{args.holdout} CDS exact-segment: F1={exact[0]:.4%} precision={exact[1]:.4%} recall={exact[2]:.4%}")
    print(f"work directory: {workdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
