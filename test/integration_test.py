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
from scripts.accuracy import ORDER, calculate, read_gff3

g100K_FASTA = ROOT / "data/truth.nip.chr5.100k.fasta"
g100K_GFF = ROOT / "data/truth.nip.chr5.100k.gff3"
g3M_FASTA = ROOT / "data/truth.nip.chr5.3m.fasta"
g3M_GFF = ROOT / "data/truth.nip.chr5.3m.gff3"


def evaluate(reference: Path, prediction: Path) -> dict:
    reference_transcripts, _ = read_gff3(str(reference), "auto")
    prediction_transcripts, _ = read_gff3(str(prediction), "auto")
    return calculate(reference_transcripts, prediction_transcripts, 0.5)


def print_metrics(title: str, result: dict) -> None:
    print(f"\n{title}")
    for name in ORDER:
        values = result[name].as_dict()
        f1 = values["f1"]
        precision = values["precision"]
        recall = values["recall"]
        format_value = lambda value: "NA" if value is None else f"{value:.2%}"
        print(
            f"{name}: F1={format_value(f1)} "
            f"precision={format_value(precision)} recall={format_value(recall)}"
        )


def assert_thresholds(title: str, result: dict, nt_min: float, exon_min: float) -> None:
    nt_f1 = result["nt"].as_dict()["f1"]
    exon_f1 = result["exon-segment"].as_dict()["f1"]
    if nt_f1 is None or exon_f1 is None or nt_f1 < nt_min or exon_f1 < exon_min:
        raise AssertionError(
            f"{title}: nt F1={nt_f1}, exon-segment F1={exon_f1}; "
            f"required nt >= {nt_min:.2%}, exon >= {exon_min:.2%}"
        )


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
        g100k_prediction = Path(directory) / "nip-100k-prediction.gff3"
        g3m_prediction = Path(directory) / "nip-3m-prediction.gff3"
        command = ["./surf", "-p", str(args.threads), "-m", str(model)]
        started = time.perf_counter()
        try:
            subprocess.run(
                command + ["-i", str(args.iterations), str(g100K_FASTA), str(g100K_GFF)],
                cwd=ROOT,
                check=True,
                timeout=6000,
            )
        except subprocess.TimeoutExpired:
            raise AssertionError("100k training exceeded 6000 seconds") from None
        elapsed = time.perf_counter() - started
        if elapsed > 6000:
            raise AssertionError(f"100k training took {elapsed:.2f}s")
        print(f"Training: {elapsed:.2f}s, {args.threads} threads", flush=True)
        with g100k_prediction.open("w") as output:
            subprocess.run(
                command + [str(g100K_FASTA)], cwd=ROOT, check=True, stdout=output
            )
        g100k_result = evaluate(g100K_GFF, g100k_prediction)
        with g3m_prediction.open("w") as output:
            subprocess.run(
                command + [str(g3M_FASTA)], cwd=ROOT, check=True, stdout=output
            )
        g3m_result = evaluate(g3M_GFF, g3m_prediction)

    print_metrics("100k self-prediction", g100k_result)
    assert_thresholds("100k self-prediction", g100k_result, 0.85, 0.70)
    print_metrics("100k -> 3m generalization", g3m_result)
    assert_thresholds("100k -> 3m generalization", g3m_result, 0.80, 0.65)


if __name__ == "__main__":
    main()
