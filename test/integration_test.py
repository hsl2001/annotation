#!/usr/bin/env python3
"""Train and self-predict the bundled 100k genome as an integration test."""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "data/truth.nip.chr5.100k.fasta"
GFF = ROOT / "data/truth.nip.chr5.100k.gff3"


def run(command: list[str], cwd: Path, capture: bool = False) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        capture_output=capture,
    )
    return result.stdout + result.stderr if capture else ""


def metric(report: str, name: str) -> float:
    match = re.search(rf"^{re.escape(name)}: F1=([0-9.]+)%", report, re.MULTILINE)
    if not match:
        raise AssertionError(f"missing {name} metric:\n{report}")
    return float(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--iterations",
        type=int,
        default=100,
        help="maximum CRF iterations (default: 100)",
    )
    args = parser.parse_args()
    if args.iterations < 1:
        raise SystemExit("iterations must be positive")

    run(["make", "-s"], ROOT)
    with tempfile.TemporaryDirectory(prefix="surf-integration-") as directory:
        work = Path(directory)
        model = work / "model.crf"
        prediction = work / "prediction.gff3"
        run(
            [
                "./surf",
                "-p",
                "1",
                "-i",
                str(args.iterations),
                "-m",
                str(model),
                str(FASTA),
                str(GFF),
            ],
            ROOT,
        )
        with prediction.open("w") as output:
            subprocess.run(
                ["./surf", "-m", str(model), str(FASTA)],
                cwd=ROOT,
                check=True,
                stdout=output,
                text=True,
            )
        report = run(
            [
                "python3",
                "scripts/accuracy.py",
                "-r",
                str(GFF),
                "-q",
                str(prediction),
            ],
            ROOT,
            capture=True,
        )
        gene_count = sum(
            line.split("\t")[2] == "gene"
            for line in prediction.read_text().splitlines()
            if line and not line.startswith("#")
        )

    exon_f1 = metric(report, "exon-segment")
    nucleotide_f1 = metric(report, "nt")
    if gene_count < 10 or nucleotide_f1 < 80.0:
        raise AssertionError(
            f"self-prediction is too weak "
            f"(genes={gene_count}, exon F1={exon_f1:.2f}%, "
            f"nt F1={nucleotide_f1:.2f}%):\n{report}"
        )
    print(report, end="")


if __name__ == "__main__":
    main()
