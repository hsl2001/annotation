#!/usr/bin/env python3
"""Exon-segment, locus and nucleotide accuracy of a query GFF3 against a reference.

Required metrics: exon-segment, locus, nt.  Extra: CDS / five_prime_UTR /
three_prime_UTR segments.  Standard library only.
"""
import argparse
from collections import defaultdict
import gzip
import sys

FEATURES = ("exon", "cds", "five_prime_utr", "three_prime_utr")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def attributes(text):
    pairs = (field.strip().partition("=") for field in text.split(";"))
    return {key: value for key, sep, value in pairs if sep}


def read_models(path):
    """Return {feature: {transcript: (gene, seqid, strand, [(start, end)])}}."""
    genes, parts = {}, {feature: defaultdict(set) for feature in FEATURES}
    where = {}
    with open_text(path) as stream:
        for number, line in enumerate(stream, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected 9 columns")
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"{path}:{number}: invalid interval")
            kind, info = fields[2].lower(), attributes(fields[8])
            if kind in ("mrna", "transcript") and "ID" in info:
                genes[info["ID"]] = info.get("Parent", info["ID"]).split(",")[0]
            elif kind in parts:
                if fields[6] not in "+-" or len(fields[6]) != 1:
                    raise ValueError(f"{path}:{number}: invalid strand")
                for parent in filter(None, info.get("Parent", "").split(",")):
                    parts[kind][parent].add((start, end))
                    where[parent] = (fields[0], fields[6])
    return {feature: {t: (genes.get(t, t),) + where[t] + (sorted(s),)
                      for t, s in models.items()}
            for feature, models in parts.items()}


def scores(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def set_scores(reference, query):
    return scores(len(reference & query), len(query - reference), len(reference - query))


def segments(models):
    return {(seqid, strand, s, e) for _, seqid, strand, spans in models.values()
            for s, e in spans}


def loci(models):
    span = {}
    for gene, seqid, strand, spans in models.values():
        key = (gene, seqid, strand)
        lo, hi = spans[0][0], max(e for _, e in spans)
        old = span.get(key, (lo, hi))
        span[key] = (min(lo, old[0]), max(hi, old[1]))
    return {(seqid, strand) + bounds for (_, seqid, strand), bounds in span.items()}


def merged(models):
    """Union of intervals per (seqid, strand), sorted and non-overlapping."""
    grouped = defaultdict(list)
    for _, seqid, strand, spans in models.values():
        grouped[seqid, strand].extend(spans)
    result = {}
    for key, spans in grouped.items():
        out = []
        for s, e in sorted(spans):
            if out and s <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        result[key] = out
    return result


def nucleotide_scores(reference, query):
    ref, qry = merged(reference), merged(query)
    size = lambda groups: sum(e - s + 1 for spans in groups.values() for s, e in spans)
    overlap = 0
    for key, a in ref.items():
        b, i, j = qry.get(key, []), 0, 0
        while i < len(a) and j < len(b):
            overlap += max(0, min(a[i][1], b[j][1]) - max(a[i][0], b[j][0]) + 1)
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
    return scores(overlap, size(qry) - overlap, size(ref) - overlap)


def metrics(reference_path, query_path):
    reference, query = read_models(reference_path), read_models(query_path)
    result = {f"{f}-segment": set_scores(segments(reference[f]), segments(query[f]))
              for f in FEATURES}
    result["locus"] = set_scores(loci(reference["exon"]), loci(query["exon"]))
    result["nt"] = nucleotide_scores(reference["exon"], query["exon"])
    return result


def main():
    parser = argparse.ArgumentParser(description=f"{__doc__ or ''}".splitlines()[0])
    parser.add_argument("-r", "--reference", required=True, help="reference GFF3(.gz)")
    parser.add_argument("-q", "--query", required=True, help="query GFF3(.gz)")
    parser.add_argument("-g", "--genome", help="genome FASTA (accepted, unused)")
    args = parser.parse_args()
    try:
        result = metrics(args.reference, args.query)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    for name, label in (("exon-segment", "exon-segment"), ("locus", "locus"), ("nt", "nt"),
                        ("cds-segment", "CDS-segment"),
                        ("five_prime_utr-segment", "five_prime_UTR-segment"),
                        ("three_prime_utr-segment", "three_prime_UTR-segment")):
        precision, recall, f1 = result[name]
        print(f"{label}: F1={f1:.2%} precision={precision:.2%} recall={recall:.2%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
