#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "matplotlib",
#     "numpy",
# ]
# ///

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import csv
from dataclasses import dataclass
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any, Callable, Iterator, Sequence, TextIO

import matplotlib
matplotlib.use("Agg")
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
import numpy as np
from numpy.typing import NDArray

plt.rcParams["svg.fonttype"] = "none"

Interval = tuple[str, int, int, str]
KmerHit = tuple[str, str]
Arr = NDArray[Any]

PLOT_NAMES: tuple[str, ...] = (
    "exons", "exons-mutated", "exon-kmer", "control", "control-aggt",
    "control-aggt-p1", "first-exons", "last-exons", "control-no-repeat",
    "control-aggt-no-repeat", "control-aggt-p1-no-repeat",
)
PLOT_CHOICES: tuple[str, ...] = ("all", *PLOT_NAMES)
REPEAT_TYPES: frozenset[str] = frozenset(("repeat_region", "mobile_element"))
BOUNDARIES: tuple[str, str] = ("start", "end")
CONTROL_REPLICATES, ABS_RADIUS = 2, 300
KMER_FLANK, KMER_SCALES, KMER_RADIUS = 100, (3, 4, 5, 6), 10
KMER_WIDTH = 2 * KMER_RADIUS + 1
MEME_TIME, MEME_MOTIFS, MEME_MINW, MEME_MAXW, MEME_ORDER = 14400, 3, 6, 15, 2
SINGLE_STRAND_DNA_ALPH = (
    'ALPHABET "DNA"\n'
    'A "Adenine" CC0000\n'
    'C "Cytosine" 0000CC\n'
    'G "Guanine" FFB300\n'
    'T "Thymine" 008000\n'
    'N "Any base" = ACGT\n'
)
DNA_COMPLEMENT = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")


def open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def fasta_records(path: Path, wanted: set[str] | None = None) -> Iterator[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as stream:
        for line in stream:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    yield name, "".join(chunks)
                name, chunks = (line[1:].split() or [""])[0], []
            elif name is not None and (wanted is None or name in wanted):
                chunks.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        yield name, "".join(chunks)


def gff_records(path: Path) -> Iterator[list[str]]:
    with open_text(path) as stream:
        for line_no, line in enumerate(stream, 1):
            if line.startswith("#"):
                if line.rstrip() == "##FASTA":
                    break
                continue
            if not (s := line.rstrip("\r\n")):
                continue
            if len(row := s.split("\t")) != 9:
                raise ValueError(f"Expected 9 GFF/GTF columns at line {line_no}")
            yield row


def gff_attributes(text: str) -> dict[str, str]:
    return {k.strip(): v for f in text.split(";") if "=" in f for k, v in (f.split("=", 1),)}


def sort_intervals(intervals: Sequence[Interval] | set[Interval]) -> list[Interval]:
    return sorted(intervals, key=lambda x: (x[2] - x[1] + 1, x))


def merge_spans(spans: list[tuple[int, int]]) -> tuple[Arr, Arr]:
    if not spans:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    spans.sort()
    merged: list[list[int]] = [[spans[0][0], spans[0][1]]]
    for s, e in spans[1:]:
        if s <= merged[-1][1] + 1:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    arr = np.asarray(merged, dtype=np.int64)
    return np.ascontiguousarray(arr[:, 0]), np.ascontiguousarray(arr[:, 1])


@dataclass(frozen=True)
class Contig:
    name: str
    length: int
    plus: int
    minus: int
    parameter: str = "Scale"

    def offset(self, strand: str) -> int:
        return self.plus if strand == "+" else self.minus

    def bounds(self, start: int, end: int, strand: str) -> tuple[int, int]:
        if strand not in ("+", "-") or start < 1 or end < start or end > self.length:
            raise ValueError(f"Invalid interval: {start}-{end} ({strand})")
        first = start - 1 if strand == "+" else self.length - end
        return first, first + end - start + 1


class CwtMatrix:
    def __init__(self, matrix: np.memmap, contigs: dict[str, Contig], widths: list[int]) -> None:
        self.matrix, self.contigs, self.widths = matrix, contigs, widths
        self.parameter = next(iter(contigs.values())).parameter
        self.f32: Arr = np.ndarray(matrix.shape + (2,), dtype=np.float32, buffer=matrix)

    def close(self) -> None:
        del self.f32
        if (mm := getattr(self.matrix, "_mmap", None)) is not None:
            mm.close()
        del self.matrix

    def power_at(self, start: int, end: int) -> Arr:
        sl = self.f32[start:end]
        return np.square(sl[..., 0]) + np.square(sl[..., 1])

    def validate(self, name: str, start: int, end: int, strand: str = "+", kind: str = "Contig") -> Interval:
        if name not in self.contigs:
            raise ValueError(f"{kind} absent from CWT matrix: {name}")
        self.contigs[name].bounds(start, end, strand)
        return (name, start, end, strand)

    @classmethod
    def from_dir(cls, directory: Path) -> CwtMatrix:
        with open(directory / "contigs.tsv") as stream:
            if (version := stream.readline().strip()) not in ("# anno-cwt-v1", "# anno-cwt-v2"):
                raise ValueError("Unsupported CWT index version")
            fields = stream.readline().rstrip().split("\t")
            parameter = "Scale" if version == "# anno-cwt-v2" else "Kernel width"
            if fields[0] != ("# scales" if version == "# anno-cwt-v2" else "# kernel_widths"):
                raise ValueError("Missing CWT parameters")
            widths = [int(v) for v in ",".join(fields[1:]).split(",") if v]
            contigs: dict[str, Contig] = {}
            expected = 0
            for row in csv.DictReader(stream, delimiter="\t"):
                name = row["seqid"]
                length, plus, minus = int(row["length"]), int(row["plus_offset"]), int(row["minus_offset"])
                if name in contigs or length <= 0 or plus != expected or minus != plus + length:
                    raise ValueError(f"Invalid contig offsets: {name}")
                contigs[name] = Contig(name, length, plus, minus, parameter)
                expected += 2 * length
        if not widths or not contigs:
            raise ValueError("Missing CWT metadata")
        path = directory / "matrix.bin"
        if not path.is_file() or path.stat().st_size != expected * len(widths) * np.dtype("<c8").itemsize:
            raise ValueError("Expected a little-endian complex64 matrix matching contigs.tsv")
        return cls(np.memmap(path, dtype="<c8", mode="r", shape=(expected, len(widths))), contigs, widths)

    @classmethod
    def load_or_create(cls, binary: Path, fasta: Path, target_dir: Path, threads: int) -> CwtMatrix:
        if (target_dir / "matrix.bin").is_file() and (target_dir / "contigs.tsv").is_file():
            return cls.from_dir(target_dir)
        target_dir.parent.mkdir(parents=True, exist_ok=True)
        exe = str(binary.resolve()) if binary.exists() else str(binary)
        subprocess.run([exe, "-p", str(max(1, threads)), "cwt", str(fasta), str(target_dir)], check=True)
        return cls.from_dir(target_dir)

    def aggregate_boundaries(self, intervals: Sequence[Interval], rows: int, radius: int = ABS_RADIUS) -> tuple[Arr, Arr]:
        ordered = sort_intervals(intervals)
        n, w_len, n_scales = len(ordered), 2 * radius + 1, len(self.widths)
        d_rows = min(rows, n)
        edges = np.linspace(0, n, d_rows + 1).astype(np.int64)
        d_tot = np.zeros((d_rows, 2, w_len, n_scales), dtype=np.float64)
        d_cnt = np.zeros((d_rows, 2, w_len, 1), dtype=np.int64)
        for r in range(d_rows):
            n_int = [0, 0]
            for name, start, end, strand in ordered[int(edges[r]) : int(edges[r + 1])]:
                c = self.contigs[name]
                off, (first, last) = c.offset(strand), c.bounds(start, end, strand)
                for b_idx, center in ((0, first), (1, last - 1)):
                    left, right = center - radius, center + radius + 1
                    if 0 <= left and right <= c.length:
                        d_tot[r, b_idx] += self.power_at(off + left, off + right)
                        n_int[b_idx] += 1
                    elif (cl := max(0, left)) < (cr := min(c.length, right)):
                        dst = slice(cl - left, cl - left + cr - cl)
                        d_tot[r, b_idx, dst] += self.power_at(off + cl, off + cr)
                        d_cnt[r, b_idx, dst] += 1
            d_cnt[r, 0] += n_int[0]
            d_cnt[r, 1] += n_int[1]
        d_tot_t, d_cnt_t = d_tot.transpose(0, 1, 3, 2), d_cnt.transpose(0, 1, 3, 2)
        tot, cnt = d_tot_t.sum(axis=0), d_cnt_t.sum(axis=0)
        mean = np.divide(tot, cnt, out=np.full(tot.shape, np.nan, dtype=np.float64), where=cnt > 0)
        np.divide(d_tot_t, d_cnt_t, out=d_tot_t, where=d_cnt_t > 0)
        d_tot_t[ np.broadcast_to(d_cnt_t == 0, d_tot_t.shape) ] = np.nan
        return d_tot_t.astype(np.float32), mean

    def kmer_maxima(
        self, intervals: Sequence[Interval], sequences: dict[str, str], mean_profile: Arr,
        bounds: dict[Interval, tuple[int, int]] | None = None,
    ) -> tuple[dict[str, list[KmerHit]], dict[str, int], int]:
        scale_idx = (
            [self.widths.index(s) for s in KMER_SCALES if s in self.widths]
            if self.parameter == "Scale"
            else list(range(min(len(self.widths), len(KMER_SCALES))))
        )
        if not scale_idx:
            raise ValueError("CWT matrix has no scales available for exon k-mer maxima")
        peaks: dict[str, int] = {}
        for b_idx, (boundary, sl_range, base_off) in enumerate((
            ("start", slice(ABS_RADIUS - KMER_FLANK, ABS_RADIUS), -KMER_FLANK),
            ("end", slice(ABS_RADIUS + 1, ABS_RADIUS + KMER_FLANK + 1), 1),
        )):
            curve = np.nanmean(mean_profile[b_idx, scale_idx, sl_range], axis=0)
            inner = curve[1:-1]
            local = (inner > curve[:-2]) & (inner >= curve[2:])
            rel = 1 + int(np.argmax(np.where(local, inner, -np.inf))) if np.any(local) else int(np.nanargmax(curve))
            peaks[boundary] = base_off + rel
            print(f"{boundary}\t{peaks[boundary]:+d}")
        windows: dict[str, list[KmerHit]] = {"start": [], "end": []}
        skipped = 0
        for iv in intervals:
            name, start, end, strand = iv
            c = self.contigs[name]
            length, (first, last) = c.length, c.bounds(start, end, strand)
            lo, hi = bounds.get(iv, (0, length)) if bounds is not None else (0, length)
            seq = sequences[name]
            for boundary, anchor, code, w_min, w_max in (
                ("start", first, "U", lo, first),
                ("end", last - 1, "D", last, hi),
            ):
                dist = peaks[boundary]
                m = anchor + dist
                ws, we, gp = (m - KMER_RADIUS, m + KMER_RADIUS + 1, m + 1) if strand == "+" else (length - m - KMER_RADIUS - 1, length - m + KMER_RADIUS, length - m)
                if m - KMER_RADIUS < w_min or m + KMER_RADIUS + 1 > w_max or ws < 0 or we > length:
                    skipped += 1
                    continue
                win = seq[ws:we] if strand == "+" else seq[ws:we].translate(DNA_COMPLEMENT)[::-1]
                if len(win) != KMER_WIDTH or any(ch not in "ACGT" for ch in win):
                    skipped += 1
                    continue
                windows[boundary].append((f"{name}:{start}-{end}:{strand}:{code}:{gp}:dist={dist:+d}bp", win))
        return windows, peaks, skipped


class GenomeData:
    def __init__(self, gff: Path, contigs: dict[str, Contig], sequences: dict[str, str] | None = None) -> None:
        self.gff, self.contigs = gff, contigs
        self.names = list(self.contigs)
        self.name_to_id = {n: i for i, n in enumerate(self.names)}
        self.c_lens = np.array([self.contigs[n].length for n in self.names], dtype=np.int64)
        self.c_offs = np.concatenate(([0], np.cumsum(self.c_lens[:-1])))
        self.exon_sets, self.kmer_bounds, self.has_repeats, self._excl_spans, self.splice_sites = self._parse_gff()
        self.sequences: dict[str, str] = sequences or {}
        self._aggt_cache: tuple[list[tuple[Arr, Arr]], list[Arr]] | None = None

    def _parse_gff(self) -> tuple[dict[str, list[Interval]], dict[Interval, tuple[int, int]], bool, dict[bool, tuple[Arr, Arr]], dict[str, list[int]]]:
        intervals, gene_ids, rrna_ids = set(), set(), set()
        parents_by_id: dict[str, list[str]] = {}
        exon_rows: list[tuple[list[str], list[str]]] = []
        rep_spans: list[tuple[int, int]] = []
        for row in gff_records(self.gff):
            ftype = row[2]
            if ftype in REPEAT_TYPES and row[0] in self.contigs:
                name, s, e = row[0], int(row[3]), int(row[4])
                if s < 1 or e < s or e > self.contigs[name].length:
                    raise ValueError(f"Invalid repeat interval: {name}:{s}-{e}")
                off = int(self.c_offs[self.name_to_id[name]])
                rep_spans.append((off + s, off + e))
            attrs = gff_attributes(row[8])
            parents = [p for p in attrs.get("Parent", "").split(",") if p]
            is_rrna = (
                ftype == "rRNA"
                or attrs.get("gbkey") == "rRNA"
                or attrs.get("biotype") == "rRNA"
                or attrs.get("gene_biotype") == "rRNA"
                or any(attrs.get(k, "").endswith("ribosomal RNA") for k in ("product", "description"))
            )
            if ident := attrs.get("ID"):
                parents_by_id[ident] = parents
                if ftype == "gene":
                    gene_ids.add(ident)
                if is_rrna:
                    rrna_ids.add(ident)
            if ftype == "exon" and not is_rrna:
                exon_rows.append((row, parents))

        def resolve_ancestors(parent: str) -> tuple[list[str], bool]:
            if parent not in anc_cache:
                seen, pending, genes, rrna = set(), [parent], [], False
                while pending:
                    if (cur := pending.pop()) not in seen:
                        seen.add(cur)
                        if cur in rrna_ids:
                            rrna = True
                        if cur in gene_ids:
                            genes.append(cur)
                        else:
                            pending.extend(parents_by_id.get(cur, ()))
                anc_cache[parent] = (genes, rrna)
            return anc_cache[parent]

        anc_cache: dict[str, tuple[list[str], bool]] = {}
        tx_exons: dict[str, list[Interval]] = defaultdict(list)
        for row, parents in exon_rows:
            resolved = [resolve_ancestors(p) for p in parents]
            if any(rrna for _, rrna in resolved):
                continue
            name, start, end, strand = row[0], int(row[3]), int(row[4]), row[6]
            if name not in self.contigs:
                raise ValueError(f"Exon contig absent from CWT matrix: {name}")
            self.contigs[name].bounds(start, end, strand)
            iv: Interval = (name, start, end, strand)
            intervals.add(iv)
            for p, (genes, _) in zip(parents, resolved):
                if genes:
                    tx_exons[p].append(iv)
        if not intervals:
            raise ValueError("No exon features found; CDS is not silently substituted for exon")

        sites: dict[str, set[int]] = defaultdict(set)
        first_exons, last_exons, multi_exons = set(), set(), set()
        kmer_bounds: dict[Interval, tuple[int, int]] = {}
        for raw_exons in tx_exons.values():
            exons = sorted(set(raw_exons), key=lambda x: (x[1], x[2]))
            name, _, _, strand = exons[0]
            plus, c = strand == "+", self.contigs[name]
            ln = c.length
            sites[name].update(s - 1 for s in (exons[0][1], exons[0][1] + 1, exons[-1][2] - 1, exons[-1][2]) if 1 <= s <= ln)
            for (_, _, e0, _), ( _, s1, _, _) in zip(exons, exons[1:]):
                if (hi := s1 - 1) >= (lo := e0 + 1):
                    sites[name].update(s - 1 for s in (lo, lo + 1, hi - 1, hi) if lo <= s <= hi and 1 <= s <= ln)
            f_pos = exons[0][1] if plus else exons[-1][2]
            l_pos = exons[-1][2] if plus else exons[0][1]
            first_exons.update(x for x in exons if (x[1] if plus else x[2]) == f_pos)
            last_exons.update(x for x in exons if (x[2] if plus else x[1]) == l_pos)
            if len(exons) < 2:
                continue
            multi_exons.update(exons)
            ordered_tx = exons if plus else exons[::-1]
            bounds_tx = [c.bounds(x[1], x[2], strand) for x in ordered_tx]
            for i, iv in enumerate(ordered_tx):
                lo = max((bounds_tx[j][1] for j in range(i) if bounds_tx[j][1] <= bounds_tx[i][0]), default=0)
                hi = min((bounds_tx[j][0] for j in range(i + 1, len(ordered_tx)) if bounds_tx[j][0] >= bounds_tx[i][1]), default=ln)
                prev_lo, prev_hi = kmer_bounds.get(iv, (0, ln))
                kmer_bounds[iv] = (max(prev_lo, lo), min(prev_hi, hi))
        exon_spans = [(int(self.c_offs[self.name_to_id[n]]) + s, int(self.c_offs[self.name_to_id[n]]) + e) for n, s, e, _ in intervals]
        excl = {False: merge_spans(list(exon_spans))}
        if rep_spans:
            excl[True] = merge_spans(exon_spans + rep_spans)
        return (
            {
                "exons": sort_intervals(intervals),
                "multi-exons": sort_intervals(multi_exons),
                "first-exons": sort_intervals(first_exons),
                "last-exons": sort_intervals(last_exons),
            },
            kmer_bounds,
            bool(rep_spans),
            excl,
            {name: sorted(pos) for name, pos in sites.items()},
        )

    def load_sequences(self, fasta: Path) -> dict[str, str]:
        if not self.sequences:
            wanted = set(self.contigs)
            seqs = {name: seq.upper() for name, seq in fasta_records(fasta, wanted)}
            if missing := wanted - seqs.keys():
                raise ValueError(f"FASTA missing sequences for CWT contigs: {', '.join(sorted(missing))}")
            for name, seq in seqs.items():
                if len(seq) != self.contigs[name].length:
                    raise ValueError(f"FASTA length for {name} ({len(seq)}) differs from CWT matrix ({self.contigs[name].length})")
            self.sequences = seqs
        return self.sequences

    def _aggt_tables(self) -> tuple[list[tuple[Arr, Arr]], list[Arr]]:
        if self._aggt_cache is None:
            codes = {b: ord(b) for b in "ACGT"}
            c_pos: list[list[Arr]] = [[], []]
            c_cid: list[list[Arr]] = [[], []]
            e_pos: list[list[Arr]] = [[], []]
            for cid, name in enumerate(self.names):
                off = int(self.c_offs[cid])
                arr = np.frombuffer(self.sequences[name].encode("ascii"), dtype=np.uint8)
                if len(arr) >= 3:
                    l, r = arr[:-1], arr[1:]
                    for s_idx, (a0, a1, d0, d1) in enumerate((("A", "G", "G", "T"), ("A", "C", "C", "T"))):
                        st = (np.nonzero((l == codes[a0]) & (r == codes[a1]))[0] + 3).astype(np.int32)
                        st = st[st <= len(arr)]
                        c_pos[s_idx].append(st)
                        c_cid[s_idx].append(np.full(len(st), cid, dtype=np.int32))
                        e_pos[s_idx].append(np.nonzero((l == codes[d0]) & (r == codes[d1]))[0].astype(np.int64) + off)
            cands = [(np.concatenate(c_pos[i]) if c_pos[i] else np.empty(0, dtype=np.int32),
                      np.concatenate(c_cid[i]) if c_cid[i] else np.empty(0, dtype=np.int32)) for i in (0, 1)]
            ends = [np.concatenate(e_pos[i]) if e_pos[i] else np.empty(0, dtype=np.int64) for i in (0, 1)]
            if cands[0][0].size == 0 and cands[1][0].size == 0:
                raise ValueError("No AG/GT-flanked non-exon sites available in the provided FASTA")
            self._aggt_cache = (cands, ends)
        return self._aggt_cache

    def sample_controls(
        self, seed: int, aggt: bool = False, length_offset: int = 0, repeat_free: bool = False, attempts: int = 200
    ) -> list[Interval]:
        cands, ends = self._aggt_tables() if aggt else ([], [])
        ex_s, ex_e = self._excl_spans[repeat_free and self.has_repeats]
        weights, rng = np.cumsum(self.c_lens, dtype=np.float64), np.random.default_rng(seed)
        needs = np.array([e - s + 1 + length_offset for _, s, e, _ in self.exon_sets["exons"]], dtype=np.int64)
        res_cid, res_st, res_en, res_str = np.full(len(needs), -1, dtype=np.int32), np.empty_like(needs), np.empty_like(needs), np.empty(len(needs), dtype=np.int8)
        rem = np.arange(len(needs), dtype=np.int64)
        for _ in range(attempts):
            if not len(rem):
                break
            strands = (rng.random(len(rem)) >= 0.5).astype(np.int8)
            if aggt:
                batches = []
                for s_idx, (pos_arr, cid_arr) in enumerate(cands):
                    if len(pos_arr) and len(sel := rem[strands == s_idx]):
                        idx = rng.integers(0, len(pos_arr), size=len(sel))
                        cid, st = cid_arr[idx], pos_arr[idx].astype(np.int64)
                        en, e_arr = st + needs[sel] - 1, ends[s_idx]
                        ge = self.c_offs[cid] + en
                        pos_idx = np.minimum(np.searchsorted(e_arr, ge), max(0, len(e_arr) - 1))
                        cond = (en <= self.c_lens[cid] - 2) & (len(e_arr) > 0) & (e_arr[pos_idx] == ge if len(e_arr) else False)
                        batches.append((sel, cid, st, en, s_idx, cond))
            else:
                cid = np.minimum(np.searchsorted(weights, rng.random(len(rem)) * weights[-1], side="right"), len(self.names) - 1)
                spans, req = self.c_lens[cid], needs[rem]
                st = np.floor(rng.random(len(rem)) * np.maximum(1, spans - req + 1)).astype(np.int64) + 1
                batches = [(rem, cid, st, st + req - 1, strands, spans >= req)]
            for sel, cid, st, en, s_val, cond in batches:
                g_st, g_en = self.c_offs[cid] + st, self.c_offs[cid] + en
                k = np.searchsorted(ex_s, g_en, side="right") - 1
                ok = cond & ((k < 0) | (ex_e[np.maximum(k, 0)] < g_st))
                hit = sel[ok]
                res_cid[hit], res_st[hit], res_en[hit] = cid[ok], st[ok], en[ok]
                res_str[hit] = s_val if isinstance(s_val, int) else s_val[ok]
            rem = rem[res_cid[rem] < 0]
        if not len(valid := np.flatnonzero(res_cid >= 0)):
            raise ValueError(f"Could not place any {'AG/GT-flanked ' if aggt else ''}non-exon control at the requested lengths")
        return [(self.names[res_cid[i]], int(res_st[i]), int(res_en[i]), "+" if res_str[i] == 0 else "-") for i in valid]

    def write_mutated_fasta(self, fasta: Path, seed: int, output: Path) -> Path:
        if not self.splice_sites:
            raise ValueError("No splice sites found")
        seqs = self.load_sequences(fasta)
        rng = np.random.default_rng(seed)
        with open(output, "w") as out:
            for name in self.names:
                buf = bytearray(seqs[name], "ascii")
                for i in self.splice_sites.get(name, ()):
                    opts = "ACGT".replace(c, "") if (c := chr(buf[i]).upper()) in "ACGT" else "ACGT"
                    buf[i] = ord(opts[int(rng.integers(len(opts)))])
                txt = buf.decode("ascii")
                out.write(f">{name}\n" + "\n".join(txt[i : i + 60] for i in range(0, len(txt), 60)) + "\n")
        return output


def _render_abs_task(task: tuple[str, Arr, Arr, list[int], str, int, Path]) -> None:
    boundary, img_vals, profile, widths, param, n_iv, out = task
    pos, ext = np.arange(-ABS_RADIUS, ABS_RADIUS + 1), (-ABS_RADIUS - 0.5, ABS_RADIUS + 0.5, 0, n_iv)
    fig, axes = plt.subplots(2, len(widths), figsize=(4.1 * len(widths), 6), gridspec_kw={"height_ratios": [4, 1.3]}, layout="constrained", squeeze=False)
    for s, w in enumerate(widths):
        ax, prof = axes[:, s]
        panel = img_vals[:, s]
        vmax = max(float(np.quantile(panel[np.isfinite(panel)], 0.995)), 1e-6)
        fig.colorbar(ax.imshow(panel, aspect="auto", origin="lower", interpolation="nearest", extent=ext, cmap="viridis", vmin=0, vmax=vmax), ax=ax, pad=0.01)
        ax.set(title=f"{param} {w} bp", ylabel="Interval rank" if s == 0 else "", ylim=ext[2:])
        prof.plot(pos, profile[s], color="#256c87", linewidth=1.3)
        prof.set(ylabel="Mean power" if s == 0 else "", ylim=(0.6, 1.2), xlabel=f"Position from {boundary} boundary (bp)")
        for p in (ax, prof):
            p.set_xticks((-ABS_RADIUS, 0, ABS_RADIUS), (f"-{ABS_RADIUS}", "0", f"+{ABS_RADIUS}"), fontsize=8)
            p.set_xlim(ext[:2])
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{boundary}.png", dpi=160)
    plt.close(fig)


class CwtPlotter:
    @classmethod
    def run_meme(cls, fasta: Path, output: Path) -> list[dict[str, Any]]:
        output.mkdir(parents=True, exist_ok=True)
        alph_file = output.parent / "dna-single.alph"
        if not alph_file.is_file():
            alph_file.write_text(SINGLE_STRAND_DNA_ALPH)
        env = {**os.environ, "TMPDIR": str(output.resolve())}
        subprocess.run(
            ["streme", "--verbosity", "1", "--oc", str(output.resolve()), "--alph", str(alph_file.resolve()),
             "--seed", "0", "--p", str(fasta.resolve()), "--minw", str(MEME_MINW), "--maxw", str(MEME_MAXW),
             "--time", str(MEME_TIME), "--order", str(MEME_ORDER), "--nmotifs", str(MEME_MOTIFS),
             "--align", "center", "--no-pgc"],
            check=True, env=env,
        )
        motifs: list[dict[str, Any]] = []
        if (html_path := output / "streme.html").is_file():
            html = html_path.read_text()
            if "// @JSON_VAR data" in html:
                start = html.index("// @JSON_VAR data")
                js = html[html.index("=", start) + 1 : html.index(";\n    </script>", start)].strip()
                data = json.loads(js)
                letters = [el.get("symbol", "") for el in data.get("alphabet", {}).get("symbols", []) if not el.get("equals")]
                idx = [letters.index(b) for b in "ACGT"] if set("ACGT").issubset(letters) else [0, 1, 2, 3]
                for m in data.get("motifs", []):
                    pwm = np.asarray(m["pwm"], dtype=np.float64)[:, idx]
                    motifs.append({"label": m.get("alt") or "STREME",
                                   "name": m.get("id", "").split("-", 1)[-1] or m.get("seed", ""),
                                   "width": int(m.get("width", len(pwm))),
                                   "sites": int(m.get("npassing", m.get("total_sites", 0))),
                                   "e_value": m.get("test_evalue") or m.get("train_pvalue") or "n/a",
                                   "probabilities": pwm})
        return motifs

    @staticmethod
    def render_meme_logo(axis: Any, motif: dict[str, Any], title: str) -> None:
        probs, pos = motif["probabilities"], np.arange(1, len(motif["probabilities"]) + 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            entropy = -np.sum(np.where(probs > 0, probs * np.log2(probs), 0.0), axis=1, keepdims=True)
        heights = probs * np.clip(2.0 - entropy, 0.0, 2.0)
        colors, font = ("#2ca25f", "#3182bd", "#f28e2b", "#d73027"), FontProperties(weight="bold")
        glyphs = {b: TextPath((0, 0), b, size=1, prop=font) for b in "ACGT"}
        bounds = {b: g.get_extents() for b, g in glyphs.items()}
        for col, x in enumerate(pos):
            bottom = 0.0
            for idx in np.argsort(heights[col]):
                if (h := float(heights[col, idx])) > 0:
                    box = bounds[base := "ACGT"[idx]]
                    tr = Affine2D().translate(-box.x0, -box.y0).scale(0.82 / box.width, h / box.height).translate(x - 0.41, bottom) + axis.transData
                    axis.add_patch(PathPatch(glyphs[base], transform=tr, color=colors[idx]))
                    bottom += h
        ticks = pos[::max(1, len(pos) // 8)]
        axis.set(ylabel="Information (bits)", ylim=(0, 2.04), yticks=(0, 1, 2), xlim=(0.5, len(pos) + 0.5), xlabel="Motif position")
        axis.set_xticks(ticks if ticks[-1] == pos[-1] else np.append(ticks, pos[-1]))
        axis.grid(axis="y", color="#888888", alpha=0.2, linewidth=0.6)
        axis.set_axisbelow(True)
        axis.set_title(title, fontsize=9)

    @classmethod
    def render_exon_kmer(cls, windows: dict[str, list[KmerHit]], peaks: dict[str, int], skipped: int, out: Path) -> None:
        out.mkdir(parents=True, exist_ok=True)
        fastas: dict[str, Path] = {}
        for b in BOUNDARIES:
            if not windows[b]:
                raise ValueError(f"No valid {KMER_WIDTH}-mer CWT maxima for exon {b}")
            with (fa := out / f"{b}-maxima-{KMER_WIDTH}mer.fasta").open("w") as stream:
                stream.writelines(f">{w[0]}\n{w[1]}\n" for w in windows[b])
            fastas[b] = fa
            print(f"Running STREME on {len(windows[b]):,} exon {b} maxima at {peaks[b]:+d} bp", file=sys.stderr, flush=True)
        with ThreadPoolExecutor(max_workers=2) as pool:
            futs = {b: pool.submit(cls.run_meme, fastas[b], out / f"meme-{b}") for b in BOUNDARIES}
            results = {b: futs[b].result() for b in BOUNDARIES}
        if skipped:
            print(f"Skipped {skipped:,} maxima with incomplete or ambiguous {KMER_WIDTH}-mers", file=sys.stderr, flush=True)
        with (out / "maxima-distances.tsv").open("w") as stream:
            stream.write("flank\tsequence_id\tdistance_bp\n")
            for b in BOUNDARIES:
                stream.writelines(f"{b}\t{w[0]}\t{peaks[b]}\n" for w in windows[b])

        fig, axes = plt.subplots(MEME_MOTIFS, 2, figsize=(17, 4.2 * MEME_MOTIFS), layout="constrained", squeeze=False)
        titles = {"start": f"Exon start upstream (peak {peaks['start']:+d} bp, n={len(windows['start']):,})",
                  "end": f"Exon end downstream (peak {peaks['end']:+d} bp, n={len(windows['end']):,})"}
        for col, b in enumerate(BOUNDARIES):
            for row in range(MEME_MOTIFS):
                ax = axes[row, col]
                if row < len(results[b]):
                    m = results[b][row]
                    cls.render_meme_logo(ax, m, f"{titles[b]}\n{m['label']}: {m['name']}\nwidth={m['width']}, sites={m['sites']:,}, E={m['e_value']}")
                else:
                    ax.axis("off")
                    ax.text(0.5, 0.5, f"{titles[b]}\nNo additional motif", ha="center", va="center", transform=ax.transAxes)
        fig.savefig(out / "motifs.png", dpi=180)
        plt.close(fig)

    @staticmethod
    def resample_mean(values: Arr, bins: int) -> Arr:
        arr = np.asarray(values, dtype=np.float64)
        if arr.ndim != 2 or not len(arr) or bins < 1:
            raise ValueError("Resampling requires nonempty position-by-scale values and positive bins")
        edges = np.linspace(0.0, len(arr), bins + 1)
        integral = np.concatenate((np.zeros((1, arr.shape[1])), np.cumsum(arr, axis=0)), axis=0)
        lower = np.floor(edges).astype(np.int64)
        return np.diff(integral[lower] + (edges - lower)[:, None] * arr[np.minimum(lower, len(arr) - 1)], axis=0) / (len(arr) / bins)

    @classmethod
    def render_roi(cls, intervals: Sequence[Interval], cwt: CwtMatrix, out: Path, flank: int, pixels: int) -> None:
        out.mkdir(parents=True, exist_ok=True)
        for name, start, end, strand in intervals:
            c = cwt.contigs[name]
            first, last = c.bounds(start, end, strand)
            before, after = min(flank, first), min(flank, c.length - last)
            s_idx, e_idx = c.offset(strand) + first - before, c.offset(strand) + last + after
            if len(values := cwt.matrix[s_idx:e_idx]) > 100000:
                raise ValueError("Use a region of at most 100,000 bp (including flanks)")
            bins, body_end = min(pixels, len(values)), before + last - first
            panels = (
                (cls.resample_mean(values.real, bins).T, "Real (bin mean)", "RdBu_r"),
                (cls.resample_mean(values.imag, bins).T, "Imaginary (bin mean)", "RdBu_r"),
                (np.log1p(cls.resample_mean(cwt.power_at(s_idx, e_idx), bins)).T, "log(1 + mean power)", "viridis"),
            )
            fig, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True, layout="constrained")
            for ax, (panel, title, cmap) in zip(axes, panels):
                vmax = max(float(np.nanmax(np.abs(panel))), 1e-6)
                fig.colorbar(ax.imshow(panel, aspect="auto", origin="lower", interpolation="nearest",
                                       extent=(0, len(values), -0.5, len(cwt.widths) - 0.5), cmap=cmap,
                                       vmin=0 if cmap == "viridis" else -vmax, vmax=vmax), ax=ax, pad=0.01)
                ax.set_yticks(range(len(cwt.widths)), cwt.widths)
                ax.set_ylabel(f"{c.parameter} (bp)")
                ax.set_title(title, loc="left", fontsize=11)
            l_b, r_b = (("start", start), ("end", end)) if strand == "+" else (("end", end), ("start", start))
            ticks: list[tuple[int, str, str]] = ([(0, f"5' flank\n{before} bp", "left")] if before else []) + [
                (before, f"{l_b[0]}\n{l_b[1]:,}", "right"), (body_end, f"{r_b[0]}\n{r_b[1]:,}", "left"),
            ] + ([(len(values), f"3' flank\n{after} bp", "right")] if after else [])
            axes[-1].set_xticks([p for p, _, _ in ticks], [lbl for _, lbl, _ in ticks])
            for lbl, (_, _, align) in zip(axes[-1].get_xticklabels(), ticks):
                lbl.set_horizontalalignment(align)
            axes[-1].tick_params(axis="x", labelsize=8, pad=6)
            axes[-1].set_xlabel("5' to 3' position; start/end mark genomic exon boundaries")
            fig.suptitle(f"{name}:{start:,}-{end:,} ({strand}) | CWT")
            fig.savefig(out / f"{name}_{start:06d}_{end:06d}_{'plus' if strand == '+' else 'minus'}.png", dpi=160)
            plt.close(fig)

    @classmethod
    def render_group(cls, names: list[str], cwt: CwtMatrix, genome: GenomeData, fasta: Path, args: argparse.Namespace) -> list[Any]:
        if (rep_names := [n for n in names if "no-repeat" in n]) and not genome.has_repeats:
            print("Skipping repeat-free plots: no repeat annotations on selected contigs", file=sys.stderr, flush=True)
            names = [n for n in names if n not in rep_names]
        if not names:
            return []
        if any(k in n for n in names for k in ("exon-kmer", "control-aggt")):
            genome.load_sequences(fasta)

        exons = genome.exon_sets["exons"]
        exons_agg = cwt.aggregate_boundaries(exons, args.rows) if "exons" in names else None
        kmer_cache: tuple[dict[str, list[KmerHit]], dict[str, int], int] | None = None
        if "exon-kmer" in names:
            multi_exons = genome.exon_sets["multi-exons"]
            if not multi_exons:
                raise ValueError("No multi-exon features found for exon-kmer")
            multi_agg = cwt.aggregate_boundaries(multi_exons, args.rows)
            kmer_cache = cwt.kmer_maxima(multi_exons, genome.sequences, multi_agg[1], genome.kmer_bounds)
        abs_tasks: list[Any] = []

        def queue_abs(agg: tuple[Arr, Arr], n_iv: int, out: Path) -> None:
            img = np.log1p(agg[0])
            abs_tasks.extend(
                (b, np.ascontiguousarray(img[:, i]), np.ascontiguousarray(agg[1][i]), cwt.widths, cwt.parameter, n_iv, out)
                for i, b in enumerate(BOUNDARIES)
            )

        for name in names:
            out_dir = args.out / name
            if name == "exon-kmer":
                assert kmer_cache is not None
                cls.render_exon_kmer(kmer_cache[0], kmer_cache[1], kmer_cache[2], out_dir)
            else:
                group = name.removesuffix("-mutated")
                base_group, repeat_free, aggt = group.removesuffix("-no-repeat"), group.endswith("-no-repeat"), name.startswith("control-aggt")
                for rep in range(1, CONTROL_REPLICATES + 1) if aggt else (0,):
                    sub_out = out_dir / f"replicate-{rep}" if aggt else out_dir
                    off = int(base_group.endswith("-p1"))
                    selected = (
                        genome.sample_controls(args.seed + off * CONTROL_REPLICATES + rep - 1, aggt=True, length_offset=off, repeat_free=repeat_free)
                        if aggt
                        else (genome.sample_controls(args.seed, repeat_free=repeat_free) if base_group == "control" else genome.exon_sets[group])
                    )
                    if repeat_free:
                        print(f"Sampled {len(selected):,}/{len(exons):,} {group} {sub_out.name}; unmatched lengths skipped after 200 attempts", file=sys.stderr, flush=True)
                    queue_abs(exons_agg if (group == "exons" and exons_agg is not None) else cwt.aggregate_boundaries(selected, args.rows), len(selected), sub_out)
        return abs_tasks

    @classmethod
    def run(cls, args: argparse.Namespace) -> None:
        args.out.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".cwt-work-", dir=args.out, ignore_cleanup_errors=True) as tmp:
            work = Path(tmp)
            if args.regions:
                if args.gff:
                    raise ValueError("-r/--region cannot be combined with -g/--gff")
                cwt = CwtMatrix.load_or_create(args.anno, args.fasta, args.cwt or (work / "original"), args.threads)
                try:
                    ivs: list[Interval] = []
                    for spec in args.regions:
                        if (bed_path := Path(spec)).is_file():
                            with open_text(bed_path) as stream:
                                for no, line in enumerate(stream, 1):
                                    if (s := line.rstrip("\r\n")) and not s.startswith(("#", "track", "browser")):
                                        if len(r := s.split("\t")) < 3:
                                            raise ValueError(f"Expected at least 3 BED columns at line {no}")
                                        ivs.append(cwt.validate(r[0], int(r[1]) + 1, int(r[2]), r[5] if len(r) > 5 and r[5] in ("+", "-") else "+", "BED contig"))
                        else:
                            if len(p := spec.split(":")) not in (2, 3) or "-" not in p[1]:
                                raise ValueError(f"Invalid region spec {spec!r}; expected SEQID:START-END[:STRAND] or BED file")
                            s0, e0 = p[1].split("-", 1)
                            ivs.append(cwt.validate(p[0], int(s0.replace(",", "")), int(e0.replace(",", "")), p[2] if len(p) == 3 else "+"))
                    if not ivs:
                        raise ValueError("No regions found")
                    cls.render_roi(ivs, cwt, args.out, args.flank, args.pixels)
                finally:
                    cwt.close()
                return

            if not args.gff:
                raise ValueError("-P/--plots requires -g/--gff")
            selected = list(PLOT_NAMES) if "all" in args.plots else [n for n in PLOT_NAMES if n in set(args.plots)]
            orig, mut = [n for n in selected if not n.endswith("-mutated")], [n for n in selected if n.endswith("-mutated")]
            tasks: list[Any] = []
            genome: GenomeData | None = None
            if orig:
                cwt = CwtMatrix.load_or_create(args.anno, args.fasta, args.cwt or (work / "original"), args.threads)
                try:
                    genome = GenomeData(args.gff, cwt.contigs)
                    print(f"Plotting {', '.join(orig)}", file=sys.stderr, flush=True)
                    tasks.extend(cls.render_group(orig, cwt, genome, args.fasta, args))
                finally:
                    cwt.close()
                del cwt
            if mut:
                if genome is None:
                    seqs = {name: seq.upper() for name, seq in fasta_records(args.fasta)}
                    if not seqs:
                        raise ValueError("No sequences found in assembly FASTA")
                    genome = GenomeData(args.gff, {n: Contig(n, len(s), 0, len(s)) for n, s in seqs.items()}, seqs)
                mut_fa = genome.write_mutated_fasta(args.fasta, args.seed, work / "mutated.fa")
                print(f"Plotting {', '.join(mut)}", file=sys.stderr, flush=True)
                cwt_mut = CwtMatrix.load_or_create(args.anno, mut_fa, work / "mutated", args.threads)
                try:
                    tasks.extend(cls.render_group(mut, cwt_mut, genome, mut_fa, args))
                finally:
                    cwt_mut.close()
                del cwt_mut
            del genome
        if tasks:
            if (workers := min(len(tasks), args.threads)) > 1:
                with ProcessPoolExecutor(max_workers=workers) as pool:
                    list(pool.map(_render_abs_task, tasks))
            else:
                for task in tasks:
                    _render_abs_task(task)


def int_min(minimum: int, label: str) -> Callable[[str], int]:
    def parse(text: str) -> int:
        if (val := int(text)) < minimum:
            raise argparse.ArgumentTypeError(f"Expected a {label} integer")
        return val
    return parse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(epilog=f"Available plots: {', '.join(PLOT_CHOICES)}")
    pos_int, nonneg_int = int_min(1, "positive"), int_min(0, "nonnegative")
    parser.add_argument("fasta", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("-P", "--plots", nargs="+", choices=PLOT_CHOICES, metavar="PLOT")
    mode.add_argument("-r", "--region", dest="regions", nargs="+", action="extend", metavar="REGION")
    parser.add_argument("-g", "--gff", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=Path("plots"))
    parser.add_argument("-c", "--cwt", type=Path)
    parser.add_argument("--anno", type=Path, default=Path("./bin/cpu/anno"))
    parser.add_argument("-p", "--threads", type=pos_int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("-n", "--rows", type=pos_int, default=1500)
    parser.add_argument("-s", "--seed", type=int, default=42)
    parser.add_argument("-f", "--flank", type=nonneg_int, default=10)
    parser.add_argument("--pixels", type=pos_int, default=1600)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        CwtPlotter.run(args)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
