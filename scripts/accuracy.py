#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict, deque
from collections.abc import Hashable
from dataclasses import dataclass, field
from typing import TypeVar
from urllib.parse import unquote

Span = tuple[int, int]
Loc = tuple[str, str]
IntervalKey = tuple[str, str, int, int]
Key = TypeVar("Key", bound=Hashable)


FEATURE_TYPES = {
    "exon": "exon",
    "cds": "CDS",
    "five_prime_utr": "five_prime_UTR",
    "5_prime_utr": "five_prime_UTR",
    "5utr": "five_prime_UTR",
    "three_prime_utr": "three_prime_UTR",
    "3_prime_utr": "three_prime_UTR",
    "3utr": "three_prime_UTR",
    "utr": "UTR_untyped",
    "5'-utr": "five_prime_UTR",
    "3'-utr": "three_prime_UTR",
    "start_codon": "start_codon",
    "stop_codon": "stop_codon",
}
TRANSCRIPT_TYPES = {
    "mrna", "transcript", "primary_transcript", "ncrna", "lnc_rna",
    "trna", "rrna", "mirna", "snrna", "snorna", "pseudogenic_transcript",
}
FEATURES = ("exon", "CDS", "five_prime_UTR", "three_prime_UTR", "UTR")


@dataclass
class Transcript:
    tid: str
    gid: str
    seqid: str
    strand: str
    raw: dict[str, list[Span]] = field(default_factory=lambda: defaultdict(list))
    spans: dict[str, list[Span]] = field(default_factory=dict)


@dataclass(frozen=True)
class Counts:
    tp: int
    fp: int
    fn: int

    @property
    def reference(self) -> int:
        return self.tp + self.fn

    @property
    def prediction(self) -> int:
        return self.tp + self.fp

    def as_dict(self) -> dict:
        if not self.reference:
            return {"tp": self.tp, "fp": self.fp, "fn": self.fn,
                    "reference": 0, "prediction": self.prediction,
                    "precision": None, "recall": None, "f1": None}
        p = self.tp / self.prediction if self.prediction else 0.0
        r = self.tp / self.reference
        f1 = 2 * p * r / (p + r) if p + r else 0.0
        return {"tp": self.tp, "fp": self.fp, "fn": self.fn,
                "reference": self.reference, "prediction": self.prediction,
                "precision": p, "recall": r, "f1": f1}


def parse_attributes(s: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for piece in s.split(";"):
        if "=" not in piece:
            continue
        key, raw_value = piece.strip().split("=", 1)
        out[key] = raw_value
    return out


def attr_value(attrs: dict[str, str], name: str) -> str:
    return unquote(attrs.get(name, ""))


def parents(attrs: dict[str, str]) -> list[str]:
    return [unquote(p) for p in attrs.get("Parent", "").split(",") if p]


def merge(spans: list[Span]) -> list[Span]:
    result: list[list[int]] = []
    for a, b in sorted(spans):
        if result and a <= result[-1][1] + 1:
            result[-1][1] = max(b, result[-1][1])
        else:
            result.append([a, b])
    return [(a, b) for a, b in result]


def finalize_transcript(t: Transcript, utr_source: str) -> dict[str, int]:
    raw = t.raw
    cds = merge(raw.get("CDS", []))
    has_exon_rows = bool(raw.get("exon"))
    explicit_5 = merge(raw.get("five_prime_UTR", []))
    explicit_3 = merge(raw.get("three_prime_UTR", []))
    generic_utr = merge(raw.get("UTR_untyped", []))
    codons = merge(raw.get("start_codon", []) + raw.get("stop_codon", []))
    exons = merge(raw["exon"] if has_exon_rows else
                  cds + codons + explicit_5 + explicit_3 + generic_utr)

    inferred5: list[Span] = []
    inferred3: list[Span] = []
    untyped_unclassified = 0
    if cds:
        coding_left = min(a for a, _ in cds + codons)
        coding_right = max(b for _, b in cds + codons)

        def side_intervals(spans: list[Span]) -> tuple[list[Span], list[Span], int]:
            left: list[Span] = []
            right: list[Span] = []
            ambiguous = 0
            for a, b in spans:
                if a < coding_left:
                    left.append((a, min(b, coding_left - 1)))
                if b > coding_right:
                    right.append((max(a, coding_right + 1), b))
                if max(a, coding_left) <= min(b, coding_right):
                    ambiguous += 1
            if t.strand == "-":
                return right, left, ambiguous
            return left, right, ambiguous

        g5, g3, untyped_unclassified = side_intervals(generic_utr)
        explicit_5 += g5
        explicit_3 += g3
        if utr_source == "auto":
            inferred5, inferred3, _ = side_intervals(exons)
    else:
        untyped_unclassified = len(generic_utr)

    five = merge(explicit_5 + inferred5)
    three = merge(explicit_3 + inferred3)
    t.spans = {
        "exon": exons,
        "CDS": cds,
        "five_prime_UTR": five,
        "three_prime_UTR": three,
        "UTR": merge(five + three),
    }
    return {
        "transcripts_with_explicit_exons": int(has_exon_rows),
        "transcripts_with_CDS": int(bool(cds)),
        "explicit_5prime_rows": len(raw.get("five_prime_UTR", [])),
        "explicit_3prime_rows": len(raw.get("three_prime_UTR", [])),
        "generic_UTR_rows": len(raw.get("UTR_untyped", [])),
        "inferred_5prime_segments": len(inferred5),
        "inferred_3prime_segments": len(inferred3),
        "unclassified_generic_UTR_segments": untyped_unclassified,
    }


def read_gff3(path: str, utr_source: str) -> tuple[list[Transcript], dict]:
    op = gzip.open if path.endswith(".gz") else open
    transcripts: dict[str, tuple[str, Loc]] = {}
    gene_ids: set[str] = set()
    child_rows: list[tuple[int, str, Loc, Span, list[str]]] = []
    seqids: set[str] = set()

    with op(path, "rt", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_no}: expected 9 tab-separated fields")
            seqid, _, raw_type, start_str, end_str, _, strand, _, attr_str = fields
            try:
                start, end = int(start_str), int(end_str)
            except ValueError as exc:
                raise ValueError(f"{path}:{line_no}: invalid feature coordinates") from exc
            if start < 1 or end < start:
                raise ValueError(f"{path}:{line_no}: invalid interval {start}-{end}")
            feature_type = raw_type.lower()
            is_relevant = (feature_type in FEATURE_TYPES or
                           feature_type in TRANSCRIPT_TYPES or
                           feature_type in ("gene", "pseudogene"))
            if not is_relevant:
                continue
            if strand not in ("+", "-"):
                if feature_type in FEATURE_TYPES:
                    raise ValueError(f"{path}:{line_no}: feature {raw_type} needs + or - strand")
                continue
            seqids.add(seqid)
            attrs = parse_attributes(attr_str)
            fid = attr_value(attrs, "ID")
            loc = (seqid, strand)
            if feature_type in ("gene", "pseudogene") and fid:
                gene_ids.add(fid)
            elif feature_type in TRANSCRIPT_TYPES and fid:
                ps = parents(attrs)
                gid = ps[0] if ps else fid
                if fid in transcripts and transcripts[fid] != (gid, loc):
                    raise ValueError(f"{path}:{line_no}: conflicting transcript ID {fid!r}")
                transcripts[fid] = (gid, loc)
            elif feature_type in FEATURE_TYPES:
                ps = parents(attrs)
                if not ps:
                    raise ValueError(f"{path}:{line_no}: {raw_type} lacks Parent=; "
                                     "please convert to hierarchical GFF3 first")
                child_rows.append((line_no, FEATURE_TYPES[feature_type], loc,
                                   (start, end), ps))

    models: dict[str, Transcript] = {}
    for line_no, kind, loc, span, ps in child_rows:
        for pid in ps:
            if pid in gene_ids and pid not in transcripts:
                tid, gid = f"__direct_gene__:{pid}", pid
            elif pid in transcripts:
                tid, gid = pid, transcripts[pid][0]
                if transcripts[pid][1] != loc:
                    raise ValueError(f"{path}:{line_no}: {pid!r} location/strand "
                                     "conflicts with parent transcript")
            else:
                tid, gid = pid, pid
            if tid not in models:
                models[tid] = Transcript(tid, gid, loc[0], loc[1])
            t = models[tid]
            if (t.seqid, t.strand) != loc:
                raise ValueError(f"{path}:{line_no}: transcript {tid!r} spans "
                                 "different contigs or strands")
            t.raw[kind].append(span)

    stats: Counter[str] = Counter()
    for t in models.values():
        stats.update(finalize_transcript(t, utr_source))
    stats["transcripts_scored"] = len(models)
    stats["unique_seqids"] = len(seqids)
    return list(models.values()), dict(stats)


def filter_protein_coding(transcripts: list[Transcript]) -> list[Transcript]:
    return [t for t in transcripts if t.spans["CDS"]]


def segment_keys(transcripts: list[Transcript], kind: str) -> set[IntervalKey]:
    return {(t.seqid, t.strand, a, b) for t in transcripts
            for a, b in t.spans[kind]}


def exact_counts(ref: set[Key] | Counter[Key], pred: set[Key] | Counter[Key]) -> Counts:
    if isinstance(ref, set) and isinstance(pred, set):
        tp = len(ref & pred)
        return Counts(tp, len(pred) - tp, len(ref) - tp)
    if isinstance(ref, Counter) and isinstance(pred, Counter):
        tp = sum((ref & pred).values())
        return Counts(tp, sum(pred.values()) - tp, sum(ref.values()) - tp)
    raise TypeError("ref and pred must both be sets or both be Counters")


def union_coverage(transcripts: list[Transcript], kind: str) -> dict[Loc, list[Span]]:
    by_loc: dict[Loc, list[Span]] = defaultdict(list)
    for t in transcripts:
        by_loc[t.seqid, t.strand].extend(t.spans[kind])
    return {loc: merge(spans) for loc, spans in by_loc.items()}


def base_counts(ref: list[Transcript], pred: list[Transcript], kind: str) -> Counts:
    a = union_coverage(ref, kind)
    b = union_coverage(pred, kind)
    nref = sum(y - x + 1 for spans in a.values() for x, y in spans)
    npred = sum(y - x + 1 for spans in b.values() for x, y in spans)
    overlap = 0
    for loc, xspans in a.items():
        yspans = b.get(loc, [])
        i = j = 0
        while i < len(xspans) and j < len(yspans):
            x0, x1 = xspans[i]
            y0, y1 = yspans[j]
            overlap += max(0, min(x1, y1) - max(x0, y0) + 1)
            if x1 <= y1:
                i += 1
            else:
                j += 1
    return Counts(overlap, npred - overlap, nref - overlap)


def intron_keys(transcripts: list[Transcript]) -> set[IntervalKey]:
    out: set[IntervalKey] = set()
    for t in transcripts:
        for (_, prev_end), (next_start, _) in zip(t.spans["exon"], t.spans["exon"][1:]):
            if prev_end + 1 < next_start:
                out.add((t.seqid, t.strand, prev_end + 1, next_start - 1))
    return out


def junction_sites(introns: set[IntervalKey], site: str) -> set[tuple[str, str, int]]:
    output = set()
    for seqid, strand, a, b in introns:
        pos = (a if strand == "+" else b) if site == "donor" else (b if strand == "+" else a)
        output.add((seqid, strand, pos))
    return output


def chain_counts(ref: list[Transcript], pred: list[Transcript], kind: str) -> Counts:
    def chains(xs: list[Transcript]) -> Counter:
        return Counter((t.seqid, t.strand, tuple(t.spans[kind])) for t in xs
                       if t.spans[kind])
    return exact_counts(chains(ref), chains(pred))


def full_transcript_counts(ref: list[Transcript], pred: list[Transcript]) -> Counts:
    kinds = ("exon", "CDS", "five_prime_UTR", "three_prime_UTR")

    def models(xs: list[Transcript]) -> Counter:
        return Counter((t.seqid, t.strand, *(tuple(t.spans[k]) for k in kinds))
                       for t in xs if t.spans["exon"])

    return exact_counts(models(ref), models(pred))


def gene_loci(transcripts: list[Transcript]) -> list[IntervalKey]:
    genes: dict[tuple[str, str, str], Span] = {}
    for t in transcripts:
        spans = t.spans["exon"]
        if not spans:
            continue
        key = (t.gid, t.seqid, t.strand)
        a, b = spans[0][0], spans[-1][1]
        if key in genes:
            c, d = genes[key]
            genes[key] = (min(a, c), max(b, d))
        else:
            genes[key] = (a, b)
    return [(seqid, strand, a, b)
            for (_, seqid, strand), (a, b) in genes.items()]


def locus_overlap_counts(ref: list[IntervalKey], pred: list[IntervalKey],
                         min_reciprocal: float) -> Counts:
    rlocs: dict[Loc, list[tuple[int, int, int]]] = defaultdict(list)
    qlocs: dict[Loc, list[tuple[int, int, int]]] = defaultdict(list)
    for i, (seqid, strand, a, b) in enumerate(ref):
        rlocs[seqid, strand].append((a, b, i))
    for j, (seqid, strand, a, b) in enumerate(pred):
        qlocs[seqid, strand].append((a, b, j))

    adj: list[list[int]] = [[] for _ in ref]
    for loc, rs in rlocs.items():
        qs = sorted(qlocs.get(loc, []))
        rs = sorted(rs)
        left = 0
        for a, b, i in rs:
            while left < len(qs) and qs[left][1] < a:
                left += 1
            pos = left
            while pos < len(qs):
                qa, qb, j = qs[pos]
                if qa > b:
                    break
                overlap = min(b, qb) - max(a, qa) + 1
                if (overlap > 0 and overlap / (b - a + 1) >= min_reciprocal
                        and overlap / (qb - qa + 1) >= min_reciprocal):
                    adj[i].append(j)
                pos += 1

    pair_r = [-1] * len(ref)
    pair_q = [-1] * len(pred)
    dist = [-1] * len(ref)

    def bfs() -> bool:
        q = deque()
        for u in range(len(ref)):
            dist[u] = 0 if pair_r[u] == -1 else -1
            if pair_r[u] == -1:
                q.append(u)
        reachable = False
        while q:
            u = q.popleft()
            for v in adj[u]:
                mate = pair_q[v]
                if mate == -1:
                    reachable = True
                elif dist[mate] == -1:
                    dist[mate] = dist[u] + 1
                    q.append(mate)
        return reachable

    def dfs(u: int) -> bool:
        for v in adj[u]:
            mate = pair_q[v]
            if mate == -1 or (dist[mate] == dist[u] + 1 and dfs(mate)):
                pair_r[u] = v
                pair_q[v] = u
                return True
        dist[u] = -1
        return False

    matched = 0
    while bfs():
        for u in range(len(ref)):
            if pair_r[u] == -1 and dfs(u):
                matched += 1
    return Counts(matched, len(pred) - matched, len(ref) - matched)


def calculate(ref: list[Transcript], pred: list[Transcript],
              min_reciprocal: float) -> dict[str, Counts]:
    result: dict[str, Counts] = {}
    for kind in FEATURES:
        result[f"{kind}-segment"] = exact_counts(segment_keys(ref, kind), segment_keys(pred, kind))
        result[f"{kind}-nt"] = base_counts(ref, pred, kind)
    result["nt"] = result["exon-nt"]
    ref_introns, pred_introns = intron_keys(ref), intron_keys(pred)
    result["intron-segment"] = exact_counts(ref_introns, pred_introns)
    result["donor-site"] = exact_counts(junction_sites(ref_introns, "donor"),
                                        junction_sites(pred_introns, "donor"))
    result["acceptor-site"] = exact_counts(junction_sites(ref_introns, "acceptor"),
                                           junction_sites(pred_introns, "acceptor"))
    result["transcript-exon-chain"] = chain_counts(ref, pred, "exon")
    result["transcript-CDS-chain"] = chain_counts(ref, pred, "CDS")
    result["transcript-full-chain"] = full_transcript_counts(ref, pred)
    rl, ql = gene_loci(ref), gene_loci(pred)
    result["locus"] = exact_counts(Counter(rl), Counter(ql))
    result["locus-overlap"] = locus_overlap_counts(rl, ql, min_reciprocal)
    return result


ORDER = (
    "locus", "locus-overlap",
    "nt", "exon-nt", "CDS-nt", "five_prime_UTR-nt", "three_prime_UTR-nt", "UTR-nt",
    "exon-segment", "CDS-segment", "five_prime_UTR-segment", "three_prime_UTR-segment",
    "UTR-segment", "intron-segment", "donor-site", "acceptor-site",
    "transcript-exon-chain", "transcript-CDS-chain", "transcript-full-chain",
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", required=True)
    parser.add_argument("-q", "--query", required=True)
    parser.add_argument("--utr-source", choices=("auto", "explicit"), default="auto")
    parser.add_argument("--protein-coding-only", action="store_true")
    parser.add_argument("--min-reciprocal-overlap", type=float, default=0.5)
    parser.add_argument("--format", choices=("text", "tsv", "json"), default="text")
    args = parser.parse_args(argv)
    if not 0 < args.min_reciprocal_overlap <= 1:
        parser.error("--min-reciprocal-overlap must be >0 and <=1")
    try:
        ref, ref_meta = read_gff3(args.reference, args.utr_source)
        pred, pred_meta = read_gff3(args.query, args.utr_source)
        if args.protein_coding_only:
            ref = filter_protein_coding(ref)
            pred = filter_protein_coding(pred)
        result = calculate(ref, pred, args.min_reciprocal_overlap)
    except (ValueError, OSError, RecursionError) as exc:
        parser.error(str(exc))

    if args.format == "json":
        print(json.dumps({"settings": {"utr_source": args.utr_source,
                                        "protein_coding_only": args.protein_coding_only,
                                        "min_reciprocal_overlap": args.min_reciprocal_overlap},
                          "reference": {**ref_meta, "transcripts_scored": len(ref)},
                          "query": {**pred_meta, "transcripts_scored": len(pred)},
                          "metrics": {k: result[k].as_dict() for k in ORDER}}, indent=2))
    elif args.format == "tsv":
        print("metric\tf1\tprecision\trecall\ttp\tfp\tfn\treference\tprediction")
        for name in ORDER:
            d = result[name].as_dict()
            def num(x: float | None) -> str:
                return "NA" if x is None else f"{x:.8f}"
            print(f"{name}\t{num(d['f1'])}\t{num(d['precision'])}\t{num(d['recall'])}\t"
                  f"{d['tp']}\t{d['fp']}\t{d['fn']}\t{d['reference']}\t{d['prediction']}")
    else:
        headers = ("metric", "F1", "precision", "recall", "TP", "FP", "FN")
        rows = []
        for name in ORDER:
            d = result[name].as_dict()
            def pct(x: float | None) -> str:
                return "NA" if x is None else f"{x:.2%}"
            rows.append((name, pct(d["f1"]), pct(d["precision"]), pct(d["recall"]),
                         str(d["tp"]), str(d["fp"]), str(d["fn"])))
        widths = [max(len(row[i]) for row in [headers, *rows]) for i in range(len(headers))]
        print("  ".join(value.ljust(widths[i]) for i, value in enumerate(headers)))
        print("  ".join("-" * width for width in widths))
        for row in rows:
            print("  ".join(value.ljust(widths[i]) for i, value in enumerate(row)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
