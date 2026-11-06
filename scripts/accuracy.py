#!/usr/bin/env python3
"""Compare exon, locus, and nucleotide accuracy between two GFF3 files.

The report includes exon-segment, locus, and nucleotide metrics, plus CDS and
UTR segment metrics. Only the Python standard library is required.
All query features are scored, including those outside the reference bounds.
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path
from typing import TextIO

FEATURES: tuple[str, ...] = (
    "exon",
    "cds",
    "five_prime_utr",
    "three_prime_utr",
)

Span = tuple[int, int]
TranscriptModel = tuple[str, str, str, list[Span]]
Metric = tuple[float, float, float]

OUTPUT_LABELS: tuple[tuple[str, str], ...] = (
    ("exon-segment", "exon-segment"),
    ("locus", "locus"),
    ("nt", "nt"),
    ("cds-segment", "CDS-segment"),
    ("five_prime_utr-segment", "five_prime_UTR-segment"),
    ("three_prime_utr-segment", "three_prime_UTR-segment"),
)


def open_text(path: str | Path) -> TextIO:
    """Open a UTF-8 plain-text or gzip-compressed annotation."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def parse_attributes(text: str) -> dict[str, str]:
    """Parse GFF3 key/value attributes, ignoring fields without '='."""
    attributes = {}
    for field in text.split(";"):
        key, separator, value = field.strip().partition("=")
        if separator:
            attributes[key] = value
    return attributes


def read_models(path: str | Path) -> dict[str, dict[str, TranscriptModel]]:
    """Read per-transcript spans, keyed by feature and transcript ID.

    Each transcript record contains its gene ID, sequence ID, strand, and
    sorted 1-based inclusive intervals.
    """
    transcript_genes: dict[str, str] = {}
    spans_by_feature: dict[str, dict[str, set[Span]]] = {
        feature: defaultdict(set) for feature in FEATURES
    }
    locations: dict[str, tuple[str, str]] = {}

    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 columns")

            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid interval") from error
            if start < 1 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid interval")

            feature = fields[2].lower()
            attributes = parse_attributes(fields[8])
            if feature in ("mrna", "transcript") and "ID" in attributes:
                transcript_id = attributes["ID"]
                transcript_genes[transcript_id] = attributes.get(
                    "Parent", transcript_id
                ).split(",")[0]
                continue
            if feature not in spans_by_feature:
                continue

            sequence_id, strand = fields[0], fields[6]
            if strand not in ("+", "-"):
                raise ValueError(f"{path}:{line_number}: invalid strand")
            for transcript_id in filter(None, attributes.get("Parent", "").split(",")):
                spans_by_feature[feature][transcript_id].add((start, end))
                locations[transcript_id] = (sequence_id, strand)

    models: dict[str, dict[str, TranscriptModel]] = {}
    for feature, transcripts in spans_by_feature.items():
        models[feature] = {}
        for transcript_id, spans in transcripts.items():
            sequence_id, strand = locations[transcript_id]
            models[feature][transcript_id] = (
                transcript_genes.get(transcript_id, transcript_id),
                sequence_id,
                strand,
                sorted(spans),
            )
    return models


def scores(true_positive: int, false_positive: int, false_negative: int) -> Metric:
    """Return precision, recall, and F1 for the supplied counts."""
    precision = (
        true_positive / (true_positive + false_positive)
        if true_positive + false_positive
        else 0.0
    )
    recall = (
        true_positive / (true_positive + false_negative)
        if true_positive + false_negative
        else 0.0
    )
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def set_scores(reference, query) -> Metric:
    """Score exact set matches, treating each set member as one item."""
    return scores(
        len(reference & query),
        len(query - reference),
        len(reference - query),
    )


def segment_keys(models: dict[str, TranscriptModel]) -> set[tuple[str, str, int, int]]:
    """Return feature intervals keyed by sequence and strand."""
    return {
        (sequence_id, strand, start, end)
        for _, sequence_id, strand, spans in models.values()
        for start, end in spans
    }


def locus_keys(models: dict[str, TranscriptModel]) -> set[tuple[str, str, int, int]]:
    """Return the outer genomic span for each gene, sequence, and strand."""
    gene_spans: dict[tuple[str, str, str], tuple[int, int]] = {}
    for gene_id, sequence_id, strand, spans in models.values():
        key = (gene_id, sequence_id, strand)
        start = min(span_start for span_start, _ in spans)
        end = max(span_end for _, span_end in spans)
        previous = gene_spans.get(key, (start, end))
        gene_spans[key] = (min(start, previous[0]), max(end, previous[1]))
    return {
        (sequence_id, strand, start, end)
        for (_, sequence_id, strand), (start, end) in gene_spans.items()
    }


def merge_intervals(
    models: dict[str, TranscriptModel],
) -> dict[tuple[str, str], list[Span]]:
    """Union overlapping or adjacent intervals per sequence and strand."""
    grouped: dict[tuple[str, str], list[Span]] = defaultdict(list)
    for _, sequence_id, strand, spans in models.values():
        grouped[sequence_id, strand].extend(spans)

    merged: dict[tuple[str, str], list[Span]] = {}
    for key, spans in grouped.items():
        intervals: list[list[int]] = []
        for start, end in sorted(spans):
            if intervals and start <= intervals[-1][1] + 1:
                intervals[-1][1] = max(intervals[-1][1], end)
            else:
                intervals.append([start, end])
        merged[key] = [(start, end) for start, end in intervals]
    return merged


def covered_bases(intervals: dict[tuple[str, str], list[Span]]) -> int:
    """Count bases covered by a disjoint interval mapping."""
    return sum(end - start + 1 for spans in intervals.values() for start, end in spans)


def nucleotide_scores(
    reference: dict[str, TranscriptModel], query: dict[str, TranscriptModel]
) -> Metric:
    """Score exon overlap at nucleotide resolution."""
    reference_intervals = merge_intervals(reference)
    query_intervals = merge_intervals(query)
    overlap = 0

    for key, reference_spans in reference_intervals.items():
        query_spans = query_intervals.get(key, [])
        ref_index = query_index = 0
        while ref_index < len(reference_spans) and query_index < len(query_spans):
            ref_start, ref_end = reference_spans[ref_index]
            query_start, query_end = query_spans[query_index]
            overlap += max(
                0,
                min(ref_end, query_end) - max(ref_start, query_start) + 1,
            )
            if ref_end < query_end:
                ref_index += 1
            else:
                query_index += 1

    return scores(
        overlap,
        covered_bases(query_intervals) - overlap,
        covered_bases(reference_intervals) - overlap,
    )


def metrics(reference_path: str | Path, query_path: str | Path) -> dict[str, Metric]:
    """Calculate metrics over all features in the two annotations."""
    reference, query = read_models(reference_path), read_models(query_path)
    result = {
        f"{feature}-segment": set_scores(
            segment_keys(reference[feature]), segment_keys(query[feature])
        )
        for feature in FEATURES
    }
    result["locus"] = set_scores(
        locus_keys(reference["exon"]), locus_keys(query["exon"])
    )
    result["nt"] = nucleotide_scores(reference["exon"], query["exon"])
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-r", "--reference", required=True, help="reference GFF3(.gz)")
    parser.add_argument("-q", "--query", required=True, help="query GFF3(.gz)")
    args = parser.parse_args(argv)

    try:
        result = metrics(args.reference, args.query)
    except (OSError, ValueError) as error:
        parser.error(str(error))

    for name, label in OUTPUT_LABELS:
        precision, recall, f1 = result[name]
        print(f"{label}: F1={f1:.2%} precision={precision:.2%} recall={recall:.2%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
