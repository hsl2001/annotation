#!/usr/bin/env python3
"""Tune the common-scale CWT exon caller without reference leakage.

The C caller detects local maxima independently in widths 4, 5, and 8, keeps
only exact position intersections, then groups those positions into candidate
intervals. This script mirrors that rule on the exported raw CWT matrix and
optimizes caller parameters against the supplied reference annotations. It is
parameter calibration, not model training, so the reported score is in-sample.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import optuna

DEFAULT_MATRIX = pathlib.Path("Col-CC")
DEFAULT_GFF = pathlib.Path("data/TAIR12_1Feb26_HEADERED.gff3")
DEFAULT_OUT = pathlib.Path("tuning.json")


@dataclass
class PeakData:
    name: str
    length: int
    peaks: list[tuple[np.ndarray, np.ndarray]]
    reference: list[tuple[int, int]]


def load_metadata(matrix_dir: pathlib.Path):
    with (matrix_dir / "contigs.tsv").open() as stream:
        if stream.readline().strip() != "# anno-cwt-v1":
            raise ValueError("Unsupported CWT index version")
        header = stream.readline().rstrip().split("\t")
        if len(header) != 2 or header[0] != "# kernel_widths":
            raise ValueError("Missing kernel width metadata")
        widths = [int(value) for value in header[1].split(",") if value]
        fields = stream.readline().rstrip().split("\t")
        rows = [dict(zip(fields, line.rstrip().split("\t"))) for line in stream if line.strip()]
    if widths != [4, 5, 8]:
        raise ValueError(f"Expected CWT widths [4, 5, 8], found {widths}")
    total_rows = sum(int(row["length"]) * 2 for row in rows)
    matrix_path = matrix_dir / "matrix.bin"
    expected_bytes = total_rows * len(widths) * np.dtype("<c8").itemsize
    if matrix_path.stat().st_size != expected_bytes:
        raise ValueError("matrix.bin size does not match contigs.tsv")
    matrix = np.memmap(matrix_path, dtype="<c8", mode="r",
                       shape=(total_rows, len(widths)))
    return matrix, rows


def read_reference(path: pathlib.Path):
    intervals: dict[str, list[tuple[int, int]]] = {}
    opener = path.open
    if str(path).endswith(".gz"):
        import gzip
        opener = lambda: gzip.open(path, "rt")
    with opener() as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2].lower() != "exon":
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid exon coordinates") from error
            intervals.setdefault(fields[0], []).append((start, end))
    return {name: merge(values) for name, values in intervals.items()}


def merge(intervals: Iterable[tuple[int, int]]):
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def build_peak_data(matrix_dir: pathlib.Path, gff: pathlib.Path):
    matrix, rows = load_metadata(matrix_dir)
    reference = read_reference(gff)
    data = []
    expected_offset = 0
    for row in rows:
        name = row["seqid"]
        length = int(row["length"])
        plus_offset = int(row["plus_offset"])
        if plus_offset != expected_offset:
            raise ValueError(f"Invalid matrix offset for {name}")
        expected_offset += length * 2
        values = np.asarray(matrix[plus_offset:plus_offset + length])
        power = values.real.astype(np.float64) ** 2 + values.imag.astype(np.float64) ** 2
        scale_peaks = []
        for scale in range(power.shape[1]):
            signal = power[:, scale]
            local = np.zeros(length, dtype=bool)
            local[1:-1] = (signal[1:-1] >= signal[:-2]) & (signal[1:-1] >= signal[2:])
            local[:4] = False
            local[-4:] = False
            positions = np.flatnonzero(local)
            scale_peaks.append((positions.astype(np.int32), signal[positions].astype(np.float32)))
        data.append(PeakData(name, length, scale_peaks, reference.get(name, [])))
    return data


def common_peaks(item: PeakData, threshold: float):
    selected = [positions[values >= threshold]
                for positions, values in item.peaks]
    if not all(len(values) for values in selected):
        return np.empty(0, dtype=np.int32)
    return np.intersect1d(np.intersect1d(selected[0], selected[1]), selected[2])


def predict(peaks: np.ndarray, window: int, min_peaks: int):
    cursor = 0
    result = []
    while cursor < len(peaks):
        last = cursor + 1
        while last < len(peaks) and peaks[last] - peaks[cursor] <= window:
            last += 1
        if last - cursor >= min_peaks:
            result.append((int(peaks[cursor]) + 1, int(peaks[last - 1])))
            cursor = last
        else:
            cursor += 1
    return result


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


def score(data: list[PeakData], threshold: float, window: int, min_peaks: int):
    true_positive = false_positive = false_negative = 0
    for item in data:
        prediction = predict(common_peaks(item, threshold), window, min_peaks)
        shared = overlap_length(item.reference, prediction)
        true_positive += shared
        false_positive += interval_length(prediction) - shared
        false_negative += interval_length(item.reference) - shared
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else 0.0
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"f1": f1, "precision": precision, "recall": recall,
            "tp": true_positive, "fp": false_positive, "fn": false_negative}


def main():
    parser = argparse.ArgumentParser(description="Optuna tuning for common 4/5/8 CWT peaks")
    parser.add_argument("--matrix", type=pathlib.Path, default=DEFAULT_MATRIX)
    parser.add_argument("--gff", type=pathlib.Path, default=DEFAULT_GFF)
    parser.add_argument("--trials", type=int, default=30)
    parser.add_argument("--timeout", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260915)
    parser.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    _, rows = load_metadata(args.matrix)
    if not rows:
        parser.error("No contigs found in CWT metadata")
    data = build_peak_data(args.matrix, args.gff)

    sampler = optuna.samplers.TPESampler(seed=args.seed)
    study = optuna.create_study(direction="maximize", sampler=sampler)

    def objective(trial):
        threshold = trial.suggest_float("power_threshold", 0.5, 2.0)
        window = trial.suggest_int("exon_window", 50, 400, step=10)
        min_peaks = trial.suggest_int("min_exon_peaks", 1, 25)
        result = score(data, threshold, window, min_peaks)
        trial.set_user_attr("precision", result["precision"])
        trial.set_user_attr("recall", result["recall"])
        return result["f1"]

    study.optimize(objective, n_trials=args.trials,
                   timeout=args.timeout or None, show_progress_bar=False)
    best = study.best_trial
    result = score(data, best.params["power_threshold"],
                   best.params["exon_window"], best.params["min_exon_peaks"])
    output = {
        "matrix": str(args.matrix), "gff": str(args.gff),
        "n_trials": len(study.trials), "seed": args.seed,
        "best_params": best.params,
        "anno_args": ["-t", str(best.params["power_threshold"]),
                  "-w", str(best.params["exon_window"]),
                  "-m", str(best.params["min_exon_peaks"])],
        "evaluation": result,
    }
    args.out.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
