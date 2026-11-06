#!/usr/bin/env python3
"""Train full Nipponbare chr5 within 120 seconds and check self-prediction."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.accuracy import OUTPUT_LABELS, metrics

FASTA = ROOT / "data/truth.nip.chr5.fasta"
GFF = ROOT / "data/truth.nip.chr5.gff3"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--iterations", type=int, default=100)
    parser.add_argument("-p", "--threads", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()
    if args.iterations < 1 or args.threads < 1:
        raise SystemExit("iterations and threads must be positive")

    subprocess.run(["make", "-s"], cwd=ROOT, check=True)
    with tempfile.TemporaryDirectory(prefix="surf-integration-") as directory:
        model = Path(directory) / "model.crf"
        prediction = Path(directory) / "prediction.gff3"
        command = ["./surf", "-p", str(args.threads), "-m", str(model)]
        started = time.perf_counter()
        try:
            subprocess.run(
                command + ["-i", str(args.iterations), str(FASTA), str(GFF)],
                cwd=ROOT,
                check=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired:
            raise AssertionError("full chr5 training exceeded 120 seconds") from None
        elapsed = time.perf_counter() - started
        if elapsed > 120:
            raise AssertionError(f"full chr5 training took {elapsed:.2f}s (>120s)")
        print(f"Training: {elapsed:.2f}s, {args.threads} threads", flush=True)
        with prediction.open("w") as output:
            subprocess.run(command + [str(FASTA)], cwd=ROOT, check=True, stdout=output)
        result = metrics(GFF, prediction)
        gene_count = sum(
            line.split("\t")[2] == "gene"
            for line in prediction.read_text().splitlines()
            if line and not line.startswith("#")
        )

    for name, label in OUTPUT_LABELS:
        precision, recall, f1 = result[name]
        print(f"{label}: F1={f1:.2%} precision={precision:.2%} recall={recall:.2%}")
    exon_f1, nucleotide_f1 = result["exon-segment"][2], result["nt"][2]
    if gene_count < 10 or nucleotide_f1 < 0.95 or exon_f1 < 0.80:
        raise AssertionError(
            f"genes={gene_count}, exon F1={exon_f1:.2%}, nt F1={nucleotide_f1:.2%}"
        )


if __name__ == "__main__":
    main()
