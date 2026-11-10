#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import math
import os
import subprocess
import sys
import tempfile
from collections import defaultdict, deque
from collections.abc import Callable, Iterable, Iterator, Sequence
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any, Literal, TextIO

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
from numpy.typing import DTypeLike, NDArray

matplotlib.use("Agg")

plt.rcParams["svg.fonttype"] = "none"


Interval = tuple[str, int, int, str]
KmerHit = tuple[str, str]
Arr = NDArray[Any]
ArrayPair = tuple[Arr, Arr]
ExclusionSpans = dict[tuple[bool, bool], ArrayPair]
ExonSets = dict[str, list[Interval]]
KmerBounds = dict[Interval, tuple[int, int]]
SpliceSites = dict[str, list[int]]
ParsedAnnotations = tuple[ExonSets, KmerBounds, bool, ExclusionSpans, SpliceSites]
MotifTables = tuple[list[ArrayPair], list[Arr]]
AbsRenderTask = tuple[str, Arr, Arr, list[int], int, Path]

PLOT_NAMES: tuple[str, ...] = (
    "exons",
    "introns",
    "exons-mutated",
    "exon-kmer",
    "control",
    "control-aggt",
    "control-aggt-p1",
    "single-exons",
    "first-exons",
    "last-exons",
    "first-last-exons",
    "control-no-repeat",
    "control-aggt-no-repeat",
    "control-aggt-p1-no-repeat",
    "control-gtag",
    "control-gtag-no-repeat",
)
PLOT_CHOICES: tuple[str, ...] = ("all", *PLOT_NAMES)
REPEAT_TYPES: frozenset[str] = frozenset(("repeat_region", "mobile_element"))
BOUNDARIES: tuple[str, str] = ("start", "end")
CONTROL_REPLICATES = 2
ABS_RADIUS = 600
KMER_FLANK = 100
KMER_RADIUS = 10
KMER_WIDTH = 2 * KMER_RADIUS + 1
MEME_TIME = 14400
MEME_MOTIFS = 3
MEME_MINW = 6
MEME_MAXW = 15
MEME_ORDER = 2
SINGLE_STRAND_DNA_ALPH = (
    'ALPHABET "DNA"\n'
    'A "Adenine" CC0000\n'
    'C "Cytosine" 0000CC\n'
    'G "Guanine" FFB300\n'
    'T "Thymine" 008000\n'
    'N "Any base" = ACGT\n'
)
DNA_COMPLEMENT = str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN")

PERIODS: tuple[float, ...] = (2.0, 3.0, 4.0)
W0: float = 6.0
MORLET_LIMIT: float = 8.0

NODES = np.array(
    [
        0.09501250983763744,
        0.2816035507792589,
        0.4580167776572274,
        0.6178762444026438,
        0.7554044083550030,
        0.8656312023878318,
        0.9445750230732326,
        0.9894009349916499,
    ]
)

WEIGHTS = np.array(
    [
        0.1894506104550685,
        0.1826034150449236,
        0.1691565193950025,
        0.1495959888165767,
        0.1246289712555339,
        0.0951585116824928,
        0.0622535239386479,
        0.0271524594117541,
    ]
)

BASE = np.array(
    [
        (1.0 + 1.0j) / math.sqrt(2.0),
        (-1.0 - 1.0j) / math.sqrt(2.0),
        (-1.0 + 1.0j) / math.sqrt(2.0),
        (1.0 - 1.0j) / math.sqrt(2.0),
        0.0j,
    ],
    dtype=np.complex64,
)

LUT = np.full(256, 4, np.uint8)
for _code, _val in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    LUT[_code] = _val


def parallel_map(
    function: Callable[[Any], Any], items: Iterable[Any], threads: int = 1
) -> Iterator[Any]:
    """Ordered reduction with at most one outstanding chunk per worker."""
    if threads < 1:
        raise ValueError("threads must be positive")
    if threads == 1:
        yield from map(function, items)
        return
    iterator = iter(items)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        pending: deque[Future[Any]] = deque()
        for _ in range(threads):
            try:
                item = next(iterator)
            except StopIteration:
                break
            pending.append(executor.submit(function, item))
        while pending:
            yield pending.popleft().result()
            try:
                item = next(iterator)
            except StopIteration:
                continue
            pending.append(executor.submit(function, item))


def info(message: str) -> None:
    logging.getLogger(__name__).info(message)


def new_map(
    path: str | Path,
    shape: tuple[int, ...],
    dtype: DTypeLike = np.float32,
) -> np.memmap:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return np.lib.format.open_memmap(path, mode="w+", dtype=dtype, shape=shape)


@dataclass(frozen=True)
class Record:
    sid: str
    codes: NDArray[np.uint8]
    energy: NDArray[np.float32]


def open_text(path: str | Path) -> TextIO:
    """Open a UTF-8 plain-text or gzip-compressed file."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def read_fasta(path: str | Path) -> dict[str, NDArray[np.uint8]]:
    """Read FASTA records as arrays of encoded base channels."""
    records: dict[str, NDArray[np.uint8]] = {}
    name: str | None = None
    lines: list[str] = []

    def save_record() -> None:
        if name is None:
            return
        if name in records:
            raise ValueError(f"FASTA duplicated ID: {name}")
        sequence = "".join(lines).encode("ascii", errors="replace")
        codes = LUT[np.frombuffer(sequence, dtype=np.uint8)]
        if not len(codes):
            raise ValueError(f"FASTA sequence is empty: {name}")
        records[name] = codes

    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith(">"):
                save_record()
                header = line[1:].split()
                if not header:
                    raise ValueError(f"FASTA header has no ID at line {line_number}")
                name, lines = header[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(
                        f"FASTA sequence occurs before its header at line {line_number}"
                    )
                lines.append("".join(line.split()))

    save_record()
    if not records:
        raise ValueError("FASTA contains no records")
    return records


@dataclass(frozen=True)
class CWTSpec:
    periods: tuple[float, ...] = PERIODS
    w0: float = W0

    def __post_init__(self) -> None:
        if not self.periods or any(
            not math.isfinite(period) or period < 2 for period in self.periods
        ):
            raise ValueError("CWT periods must be at least 2 bp")
        if not math.isfinite(self.w0) or self.w0 <= 0:
            raise ValueError("w0 must be a finite positive value")


def cwt_morlet_integral(lower: float, upper: float, w0: float = W0) -> complex:
    lower = max(lower, -MORLET_LIMIT)
    upper = min(upper, MORLET_LIMIT)
    if lower >= upper:
        return 0.0

    squared_frequency = w0 * w0
    norm_squared = -math.expm1(-squared_frequency) + 2.0 * math.exp(
        -0.75 * squared_frequency
    ) * math.expm1(-0.25 * squared_frequency)
    norm = 1.0 / math.sqrt(math.sqrt(math.pi) * norm_squared)
    correction = math.exp(-0.5 * squared_frequency)
    segments = math.ceil((upper - lower) * 2.0)
    step = (upper - lower) / segments
    integral = 0.0

    for seg in range(segments):
        middle = lower + (seg + 0.5) * step
        for node, weight in zip(NODES, WEIGHTS):
            left = middle - node * step / 2.0
            right = middle + node * step / 2.0

            val_l = math.exp(-0.5 * left**2) * (np.exp(-1j * w0 * left) - correction)
            val_r = math.exp(-0.5 * right**2) * (np.exp(-1j * w0 * right) - correction)
            integral += weight * (step / 2.0) * (val_l + val_r)

    return norm * integral


def build_kernel(dilation: float, w0: float = W0) -> np.ndarray:
    """Integrate the corrected Morlet wavelet over bases at the given dilation."""
    radius = math.ceil(MORLET_LIMIT * dilation + 0.5)
    width = 2 * radius + 1
    kernel = np.zeros(width, dtype=np.complex128)

    for tap in range(width):
        disp = tap - radius
        kernel[tap] = math.sqrt(dilation) * cwt_morlet_integral(
            (disp - 0.5) / dilation, (disp + 0.5) / dilation, w0
        )
    return kernel


def compute_energy(
    codes: NDArray[np.uint8],
    path: str | Path,
    spec: CWTSpec,
    chunk: int,
    threads: int,
) -> np.memmap:
    n = len(codes)
    kernels = [build_kernel(p, spec.w0) for p in spec.periods]
    widths = [len(k) for k in kernels]
    num_scales = len(spec.periods)

    result = new_map(path, (n,))

    def compute_chunk(start: int) -> None:
        end = min(n, start + chunk)
        chunk_len = end - start
        powers = np.zeros((num_scales, chunk_len), dtype=np.float32)

        for scale_idx, kernel in enumerate(kernels):
            width = widths[scale_idx]
            radius = width // 2

            input_start = start - radius if start > radius else 0
            input_end = min(n, end + radius)
            input_length = input_end - input_start

            required = input_length + width - 1
            fft_size = 1 << (required - 1).bit_length()

            signal_buf = np.zeros(fft_size, dtype=np.complex128)
            signal_buf[:input_length] = BASE[codes[input_start:input_end]]

            kernel_rev = kernel[::-1]
            kernel_buf = np.zeros(fft_size, dtype=np.complex128)
            kernel_buf[:width] = kernel_rev

            sig_fft = np.fft.fft(signal_buf)
            ker_fft = np.fft.fft(kernel_buf)
            conv = np.fft.ifft(sig_fft * ker_fft)

            offset = radius - input_start
            indices = np.arange(start, end) + offset
            coeffs = conv[indices]

            powers[scale_idx] = (coeffs.real**2 + coeffs.imag**2).astype(np.float32)

        result[start:end] = powers.mean(axis=0)

    for _ in parallel_map(compute_chunk, range(0, n, chunk), threads):
        pass

    result.flush()
    return result


def digest(metadata: Any) -> str:
    return hashlib.sha256(json.dumps(metadata, sort_keys=True).encode()).hexdigest()


def prepare_records(
    fasta: Path,
    spec: CWTSpec,
    folder: Path,
    chunk: int,
    threads: int,
) -> list[Record]:
    records: list[Record] = []
    for sid, codes in read_fasta(fasta).items():
        key = digest(sid)[:16]
        info(f"CWT: {sid} ({len(codes):,} bp)")
        energy = compute_energy(
            codes, Path(folder) / f"{key}.npy", spec, chunk, threads
        )
        records.append(Record(sid, codes, energy))
    return records


def fasta_records(
    path: Path, wanted: set[str] | None = None
) -> Iterator[tuple[str, str]]:
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
    """Parse GFF3 key/value attributes, ignoring fields without '='."""
    attributes = {}
    for field in text.split(";"):
        key, separator, value = field.partition("=")
        if separator:
            attributes[key.strip()] = value
    return attributes


def sort_intervals(intervals: Sequence[Interval] | set[Interval]) -> list[Interval]:
    return sorted(intervals, key=lambda x: (x[2] - x[1] + 1, x))


def merge_spans(spans: list[tuple[int, int]]) -> tuple[Arr, Arr]:
    if not spans:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    spans.sort()
    merged: list[list[int]] = [[spans[0][0], spans[0][1]]]
    for s, e in spans[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    arr = np.asarray(merged, dtype=np.int64)
    return np.ascontiguousarray(arr[:, 0]), np.ascontiguousarray(arr[:, 1])


@dataclass(frozen=True)
class Contig:
    length: int
    plus: int
    minus: int

    def offset(self, strand: str) -> int:
        return self.plus if strand == "+" else self.minus

    def bounds(self, start: int, end: int, strand: str) -> tuple[int, int]:
        if strand not in ("+", "-") or start < 1 or end < start or end > self.length:
            raise ValueError(f"Invalid interval: {start}-{end} ({strand})")
        first = start - 1 if strand == "+" else self.length - end
        return first, first + end - start + 1


class CwtMatrix:
    def __init__(
        self, matrix: np.memmap, contigs: dict[str, Contig], widths: list[int]
    ) -> None:
        self.matrix, self.contigs, self.widths = matrix, contigs, widths

    def close(self) -> None:
        if (mm := getattr(self.matrix, "_mmap", None)) is not None:
            mm.close()
        del self.matrix

    def power_at(self, start: int, end: int) -> Arr:
        return self.matrix[start:end, None]

    def validate(
        self,
        name: str,
        start: int,
        end: int,
        strand: str = "+",
        kind: str = "Contig",
    ) -> Interval:
        if name not in self.contigs:
            raise ValueError(f"{kind} absent from CWT matrix: {name}")
        self.contigs[name].bounds(start, end, strand)
        return (name, start, end, strand)

    @classmethod
    def load(cls, fasta: Path, target_dir: Path, threads: int, chunk: int) -> CwtMatrix:
        target_dir.mkdir(parents=True, exist_ok=True)
        spec = CWTSpec()
        records = prepare_records(fasta, spec, target_dir / "energy", chunk, threads)
        contigs: dict[str, Contig] = {}
        expected = 0
        for record in records:
            length = len(record.codes)
            contigs[record.sid] = Contig(length, expected, expected + length)
            expected += 2 * length

        matrix = np.lib.format.open_memmap(
            target_dir / "surf-cwt.npy",
            mode="w+",
            dtype=np.float32,
            shape=(expected,),
        )
        for record in records:
            contig = contigs[record.sid]
            matrix[contig.plus : contig.minus] = record.energy
            matrix[contig.minus : contig.minus + len(record.codes)] = record.energy[
                ::-1
            ]
        matrix.flush()
        return cls(matrix, contigs, [1])

    def aggregate_boundaries(
        self,
        intervals: Sequence[Interval],
        rows: int,
        radius: int = ABS_RADIUS,
    ) -> ArrayPair:
        ordered = sort_intervals(intervals)
        interval_count = len(ordered)
        window_length = 2 * radius + 1
        scale_count = len(self.widths)
        row_count = min(rows, interval_count)
        row_edges = np.linspace(0, interval_count, row_count + 1).astype(np.int64)
        totals = np.zeros((row_count, 2, window_length, scale_count), dtype=np.float64)
        counts = np.zeros((row_count, 2, window_length, 1), dtype=np.int64)

        for row_index in range(row_count):
            complete_intervals = [0, 0]
            row_intervals = ordered[
                int(row_edges[row_index]) : int(row_edges[row_index + 1])
            ]
            for name, start, end, strand in row_intervals:
                contig = self.contigs[name]
                offset = contig.offset(strand)
                first, last = contig.bounds(start, end, strand)
                for boundary_index, center in ((0, first), (1, last - 1)):
                    left = center - radius
                    right = center + radius + 1
                    if 0 <= left and right <= contig.length:
                        totals[row_index, boundary_index] += self.power_at(
                            offset + left, offset + right
                        )
                        complete_intervals[boundary_index] += 1
                    else:
                        clipped_left = max(0, left)
                        clipped_right = min(contig.length, right)
                        if clipped_left >= clipped_right:
                            continue
                        destination = slice(
                            clipped_left - left,
                            clipped_left - left + clipped_right - clipped_left,
                        )
                        totals[row_index, boundary_index, destination] += self.power_at(
                            offset + clipped_left, offset + clipped_right
                        )
                        counts[row_index, boundary_index, destination] += 1

            counts[row_index, 0] += complete_intervals[0]
            counts[row_index, 1] += complete_intervals[1]

        row_totals = totals.transpose(0, 1, 3, 2)
        row_counts = counts.transpose(0, 1, 3, 2)
        total = row_totals.sum(axis=0)
        count = row_counts.sum(axis=0)
        mean_profile = np.divide(
            total,
            count,
            out=np.full(total.shape, np.nan, dtype=np.float64),
            where=count > 0,
        )
        np.divide(row_totals, row_counts, out=row_totals, where=row_counts > 0)
        missing = np.broadcast_to(row_counts == 0, row_totals.shape)
        row_totals[missing] = np.nan
        return row_totals.astype(np.float32), mean_profile

    def kmer_maxima(
        self,
        intervals: Sequence[Interval],
        sequences: dict[str, str],
        mean_profile: Arr,
        bounds: dict[Interval, tuple[int, int]] | None = None,
    ) -> tuple[dict[str, list[KmerHit]], dict[str, int], int]:
        scale_indices = list(range(len(self.widths)))
        if not scale_indices:
            raise ValueError("CWT matrix has no scales available for exon k-mer maxima")

        peaks: dict[str, int] = {}
        boundary_slices = (
            ("start", slice(ABS_RADIUS - KMER_FLANK, ABS_RADIUS), -KMER_FLANK),
            ("end", slice(ABS_RADIUS + 1, ABS_RADIUS + KMER_FLANK + 1), 1),
        )
        for boundary_index, (boundary, profile_slice, base_offset) in enumerate(
            boundary_slices
        ):
            curve = np.nanmean(
                mean_profile[boundary_index, scale_indices, profile_slice], axis=0
            )
            interior = curve[1:-1]
            local_maxima = (interior > curve[:-2]) & (interior >= curve[2:])
            if np.any(local_maxima):
                relative_peak = 1 + int(
                    np.argmax(np.where(local_maxima, interior, -np.inf))
                )
            else:
                relative_peak = int(np.nanargmax(curve))
            peaks[boundary] = base_offset + relative_peak
            print(f"{boundary}\t{peaks[boundary]:+d}")

        windows: dict[str, list[KmerHit]] = {"start": [], "end": []}
        skipped = 0
        for interval in intervals:
            name, start, end, strand = interval
            contig = self.contigs[name]
            length = contig.length
            first, last = contig.bounds(start, end, strand)
            window_min, window_max = (
                bounds.get(interval, (0, length)) if bounds is not None else (0, length)
            )
            sequence = sequences[name]
            boundaries = (
                ("start", first, "U", window_min, first),
                ("end", last - 1, "D", last, window_max),
            )
            for boundary, anchor, code, window_min, window_max in boundaries:
                distance = peaks[boundary]
                motif_position = anchor + distance
                if strand == "+":
                    window_start = motif_position - KMER_RADIUS
                    window_end = motif_position + KMER_RADIUS + 1
                    genomic_position = motif_position + 1
                else:
                    window_start = length - motif_position - KMER_RADIUS - 1
                    window_end = length - motif_position + KMER_RADIUS
                    genomic_position = length - motif_position

                if (
                    motif_position - KMER_RADIUS < window_min
                    or motif_position + KMER_RADIUS + 1 > window_max
                    or window_start < 0
                    or window_end > length
                ):
                    skipped += 1
                    continue

                kmer = sequence[window_start:window_end]
                if strand == "-":
                    kmer = kmer.translate(DNA_COMPLEMENT)[::-1]
                if len(kmer) != KMER_WIDTH or any(base not in "ACGT" for base in kmer):
                    skipped += 1
                    continue
                header = (
                    f"{name}:{start}-{end}:{strand}:{code}:{genomic_position}:"
                    f"dist={distance:+d}bp"
                )
                windows[boundary].append((header, kmer))

        return windows, peaks, skipped


class GenomeData:
    def __init__(
        self,
        gff: Path,
        contigs: dict[str, Contig],
        sequences: dict[str, str] | None = None,
    ) -> None:
        self.gff = gff
        self.contigs = contigs
        self.names = list(contigs)
        self.name_to_id = {name: index for index, name in enumerate(self.names)}
        self.c_lens = np.array(
            [contigs[name].length for name in self.names], dtype=np.int64
        )
        self.c_offs = np.concatenate(([0], np.cumsum(self.c_lens[:-1])))
        (
            self.exon_sets,
            self.kmer_bounds,
            self.has_repeats,
            self._excl_spans,
            self.splice_sites,
        ) = self._parse_gff()
        self.sequences = sequences or {}
        self._motif_cache: dict[bool, MotifTables] = {}

    def _global_offset(self, contig: str) -> int:
        return int(self.c_offs[self.name_to_id[contig]])

    def _parse_gff(self) -> ParsedAnnotations:
        intervals: set[Interval] = set()
        gene_ids: set[str] = set()
        rrna_ids: set[str] = set()
        parents_by_id: dict[str, list[str]] = {}
        exon_rows: list[tuple[list[str], list[str]]] = []
        rep_spans: list[tuple[int, int]] = []
        for row in gff_records(self.gff):
            ftype = row[2]
            if ftype in REPEAT_TYPES and row[0] in self.contigs:
                name, s, e = row[0], int(row[3]), int(row[4])
                if s < 1 or e < s or e > self.contigs[name].length:
                    raise ValueError(f"Invalid repeat interval: {name}:{s}-{e}")
                off = self._global_offset(name)
                rep_spans.append((off + s, off + e))
            attrs = gff_attributes(row[8])
            parents = [p for p in attrs.get("Parent", "").split(",") if p]
            is_rrna = (
                ftype == "rRNA"
                or attrs.get("gbkey") == "rRNA"
                or attrs.get("biotype") == "rRNA"
                or attrs.get("gene_biotype") == "rRNA"
                or any(
                    attrs.get(key, "").endswith("ribosomal RNA")
                    for key in ("product", "description")
                )
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
            raise ValueError(
                "No exon features found; CDS is not silently substituted for exon"
            )

        sites: dict[str, set[int]] = defaultdict(set)
        single_exons: set[Interval] = set()
        first_exons: set[Interval] = set()
        last_exons: set[Interval] = set()
        first_last_start: set[Interval] = set()
        first_last_end: set[Interval] = set()
        multi_exons: set[Interval] = set()
        introns: set[Interval] = set()
        kmer_bounds: dict[Interval, tuple[int, int]] = {}
        for raw_exons in tx_exons.values():
            exons = sorted(set(raw_exons), key=lambda x: (x[1], x[2]))
            name, _, _, strand = exons[0]
            plus, c = strand == "+", self.contigs[name]
            ln = c.length
            terminal_sites = (
                exons[0][1],
                exons[0][1] + 1,
                exons[-1][2] - 1,
                exons[-1][2],
            )
            sites[name].update(site - 1 for site in terminal_sites if 1 <= site <= ln)
            for (_, _, e0, _), (_, s1, _, _) in pairwise(exons):
                if (hi := s1 - 1) >= (lo := e0 + 1):
                    introns.add((name, lo, hi, strand))
                    intron_sites = (lo, lo + 1, hi - 1, hi)
                    sites[name].update(
                        site - 1
                        for site in intron_sites
                        if lo <= site <= hi and 1 <= site <= ln
                    )
            f_pos = exons[0][1] if plus else exons[-1][2]
            l_pos = exons[-1][2] if plus else exons[0][1]
            first_exons.update(x for x in exons if (x[1] if plus else x[2]) == f_pos)
            last_exons.update(x for x in exons if (x[2] if plus else x[1]) == l_pos)
            if len(exons) < 2:
                single_exons.update(exons)
                continue
            first_last_start.update(
                x for x in exons if (x[1] if plus else x[2]) == f_pos
            )
            first_last_end.update(x for x in exons if (x[2] if plus else x[1]) == l_pos)
            multi_exons.update(exons)
            ordered_tx = exons if plus else exons[::-1]
            bounds_tx = [c.bounds(x[1], x[2], strand) for x in ordered_tx]
            for i, iv in enumerate(ordered_tx):
                lo = max(
                    (
                        bounds_tx[j][1]
                        for j in range(i)
                        if bounds_tx[j][1] <= bounds_tx[i][0]
                    ),
                    default=0,
                )
                hi = min(
                    (
                        bounds_tx[j][0]
                        for j in range(i + 1, len(ordered_tx))
                        if bounds_tx[j][0] >= bounds_tx[i][1]
                    ),
                    default=ln,
                )
                prev_lo, prev_hi = kmer_bounds.get(iv, (0, ln))
                kmer_bounds[iv] = (max(prev_lo, lo), min(prev_hi, hi))
        exon_spans = [
            (self._global_offset(name) + start, self._global_offset(name) + end)
            for name, start, end, _ in intervals
        ]
        intron_spans = [
            (self._global_offset(name) + start, self._global_offset(name) + end)
            for name, start, end, _ in introns
        ]
        excl: ExclusionSpans = {
            (False, False): merge_spans(exon_spans),
            (False, True): merge_spans(exon_spans + intron_spans),
        }
        if rep_spans:
            excl[True, False] = merge_spans(exon_spans + rep_spans)
            excl[True, True] = merge_spans(exon_spans + intron_spans + rep_spans)
        return (
            {
                "exons": sort_intervals(intervals),
                "introns": sort_intervals(introns),
                "multi-exons": sort_intervals(multi_exons),
                "single-exons": sort_intervals(single_exons),
                "first-exons": sort_intervals(first_exons),
                "last-exons": sort_intervals(last_exons),
                "first-last-start": sort_intervals(first_last_start),
                "first-last-end": sort_intervals(first_last_end),
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
                missing_names = ", ".join(sorted(missing))
                raise ValueError(
                    f"FASTA missing sequences for CWT contigs: {missing_names}"
                )
            for name, seq in seqs.items():
                expected_length = self.contigs[name].length
                if len(seq) != expected_length:
                    raise ValueError(
                        f"FASTA length for {name} ({len(seq)}) differs from "
                        f"CWT matrix ({expected_length})"
                    )
            self.sequences = seqs
        return self.sequences

    def _motif_tables(self, gtag: bool = False) -> MotifTables:
        """Cache splice-motif starts and ends for both genomic strands."""
        if gtag:
            motifs = (("G", "T", "A", "G"), ("C", "T", "A", "C"))
            start_shift, end_shift = 1, 2
        else:
            motifs = (("A", "G", "G", "T"), ("A", "C", "C", "T"))
            start_shift, end_shift = 3, 0

        if gtag not in self._motif_cache:
            base_codes = {base: ord(base) for base in "ACGT"}
            candidate_positions: list[list[Arr]] = [[], []]
            candidate_contigs: list[list[Arr]] = [[], []]
            end_positions: list[list[Arr]] = [[], []]

            for contig_index, name in enumerate(self.names):
                offset = self._global_offset(name)
                bases = np.frombuffer(
                    self.sequences[name].encode("ascii"), dtype=np.uint8
                )
                if len(bases) < 3:
                    continue

                left_bases, right_bases = bases[:-1], bases[1:]
                for strand_index, (
                    start_left,
                    start_right,
                    end_left,
                    end_right,
                ) in enumerate(motifs):
                    starts = (
                        np.nonzero(
                            (left_bases == base_codes[start_left])
                            & (right_bases == base_codes[start_right])
                        )[0]
                        + start_shift
                    ).astype(np.int32)
                    starts = starts[starts <= len(bases)]
                    candidate_positions[strand_index].append(starts)
                    candidate_contigs[strand_index].append(
                        np.full(len(starts), contig_index, dtype=np.int32)
                    )
                    ends = np.nonzero(
                        (left_bases == base_codes[end_left])
                        & (right_bases == base_codes[end_right])
                    )[0].astype(np.int64)
                    end_positions[strand_index].append(ends + end_shift + offset)

            candidates = [
                (
                    np.concatenate(candidate_positions[index])
                    if candidate_positions[index]
                    else np.empty(0, dtype=np.int32),
                    np.concatenate(candidate_contigs[index])
                    if candidate_contigs[index]
                    else np.empty(0, dtype=np.int32),
                )
                for index in range(2)
            ]
            ends = [
                np.concatenate(end_positions[index])
                if end_positions[index]
                else np.empty(0, dtype=np.int64)
                for index in range(2)
            ]
            if all(candidate[0].size == 0 for candidate in candidates):
                motif_name = "GT/AG" if gtag else "AG/GT-flanked"
                raise ValueError(
                    f"No {motif_name} sites available in the provided FASTA"
                )
            self._motif_cache[gtag] = (candidates, ends)

        return self._motif_cache[gtag]

    def sample_controls(
        self,
        seed: int,
        aggt: bool = False,
        length_offset: int = 0,
        repeat_free: bool = False,
        attempts: int = 200,
        gtag: bool = False,
    ) -> list[Interval]:
        use_motifs = aggt or gtag
        candidates, end_sites = self._motif_tables(gtag) if use_motifs else ([], [])
        end_margin = 0 if gtag else 2
        excluded_starts, excluded_ends = self._excl_spans[
            repeat_free and self.has_repeats, gtag
        ]
        genome_weights = np.cumsum(self.c_lens, dtype=np.float64)
        rng = np.random.default_rng(seed)
        feature_type = "introns" if gtag else "exons"
        required_lengths = np.array(
            [
                end - start + 1 + length_offset
                for _, start, end, _ in self.exon_sets[feature_type]
            ],
            dtype=np.int64,
        )

        result_contigs = np.full(len(required_lengths), -1, dtype=np.int32)
        result_starts = np.empty_like(required_lengths)
        result_ends = np.empty_like(required_lengths)
        result_strands = np.empty(len(required_lengths), dtype=np.int8)
        remaining = np.arange(len(required_lengths), dtype=np.int64)

        for _ in range(attempts):
            if not len(remaining):
                break
            sampled_strands = (rng.random(len(remaining)) >= 0.5).astype(np.int8)

            if use_motifs:
                batches = []
                for strand_index, (positions, contig_ids) in enumerate(candidates):
                    selected = remaining[sampled_strands == strand_index]
                    if not len(positions) or not len(selected):
                        continue

                    sample_indices = rng.integers(0, len(positions), size=len(selected))
                    contig_ids = contig_ids[sample_indices]
                    starts = positions[sample_indices].astype(np.int64)
                    ends = starts + required_lengths[selected] - 1
                    end_positions = end_sites[strand_index]
                    global_ends = self.c_offs[contig_ids] + ends
                    end_indices = np.minimum(
                        np.searchsorted(end_positions, global_ends),
                        max(0, len(end_positions) - 1),
                    )
                    end_matches = (
                        end_positions[end_indices] == global_ends
                        if len(end_positions)
                        else False
                    )
                    valid = (
                        (ends <= self.c_lens[contig_ids] - end_margin)
                        & (len(end_positions) > 0)
                        & end_matches
                    )
                    batches.append(
                        (selected, contig_ids, starts, ends, strand_index, valid)
                    )
            else:
                contig_ids = np.minimum(
                    np.searchsorted(
                        genome_weights,
                        rng.random(len(remaining)) * genome_weights[-1],
                        side="right",
                    ),
                    len(self.names) - 1,
                )
                contig_lengths = self.c_lens[contig_ids]
                lengths = required_lengths[remaining]
                starts = (
                    np.floor(
                        rng.random(len(remaining))
                        * np.maximum(1, contig_lengths - lengths + 1)
                    ).astype(np.int64)
                    + 1
                )
                batches = [
                    (
                        remaining,
                        contig_ids,
                        starts,
                        starts + lengths - 1,
                        sampled_strands,
                        contig_lengths >= lengths,
                    )
                ]

            for selected, contig_ids, starts, ends, strands, valid in batches:
                global_starts = self.c_offs[contig_ids] + starts
                global_ends = self.c_offs[contig_ids] + ends
                excluded_index = (
                    np.searchsorted(excluded_starts, global_ends, side="right") - 1
                )
                outside_exclusions = (excluded_index < 0) | (
                    excluded_ends[np.maximum(excluded_index, 0)] < global_starts
                )
                accepted = valid & outside_exclusions
                hits = selected[accepted]
                result_contigs[hits] = contig_ids[accepted]
                result_starts[hits] = starts[accepted]
                result_ends[hits] = ends[accepted]
                result_strands[hits] = (
                    strands if isinstance(strands, int) else strands[accepted]
                )

            remaining = remaining[result_contigs[remaining] < 0]

        valid_indices = np.flatnonzero(result_contigs >= 0)
        if not len(valid_indices):
            feature_name = (
                "GT/AG non-genic"
                if gtag
                else "AG/GT-flanked non-exon"
                if aggt
                else "non-exon"
            )
            raise ValueError(
                f"Could not place any {feature_name} control at the requested lengths"
            )

        return [
            (
                self.names[result_contigs[index]],
                int(result_starts[index]),
                int(result_ends[index]),
                "+" if result_strands[index] == 0 else "-",
            )
            for index in valid_indices
        ]

    def write_mutated_fasta(self, fasta: Path, seed: int, output: Path) -> Path:
        if not self.splice_sites:
            raise ValueError("No splice sites found")
        seqs = self.load_sequences(fasta)
        rng = np.random.default_rng(seed)
        with open(output, "w") as out:
            for name in self.names:
                buf = bytearray(seqs[name], "ascii")
                for index in self.splice_sites.get(name, ()):
                    current_base = chr(buf[index]).upper()
                    alternatives = (
                        "ACGT".replace(current_base, "")
                        if current_base in "ACGT"
                        else "ACGT"
                    )
                    replacement_index = int(rng.integers(len(alternatives)))
                    buf[index] = ord(alternatives[replacement_index])

                sequence = buf.decode("ascii")
                wrapped_sequence = "\n".join(
                    sequence[index : index + 60]
                    for index in range(0, len(sequence), 60)
                )
                out.write(f">{name}\n{wrapped_sequence}\n")
        return output


def _render_abs_task(task: AbsRenderTask) -> None:
    boundary, image_values, profile, widths, interval_count, output = task
    positions = np.arange(-ABS_RADIUS, ABS_RADIUS + 1)
    extent = (-ABS_RADIUS - 0.5, ABS_RADIUS + 0.5, 0, interval_count)
    figure, axes = plt.subplots(
        2,
        len(widths),
        figsize=(4.1 * len(widths), 6),
        gridspec_kw={"height_ratios": (4, 1.2)},
        layout="constrained",
        squeeze=False,
    )

    for scale_index in range(len(widths)):
        image_axis, energy_axis = axes[:, scale_index]
        panel = image_values[:, scale_index]
        finite_values = panel[np.isfinite(panel)]
        maximum = max(float(np.quantile(finite_values, 0.995)), 1e-6)
        image = image_axis.imshow(
            panel,
            aspect="auto",
            origin="lower",
            interpolation="nearest",
            extent=extent,
            cmap="viridis",
            vmin=0,
            vmax=maximum,
        )
        figure.colorbar(image, ax=image_axis, pad=0.01)
        image_axis.set(
            title="Surf CWT energy",
            ylabel="Interval rank" if scale_index == 0 else "",
            ylim=extent[2:],
        )

        energy_axis.plot(
            positions, profile[scale_index], color="#256c87", linewidth=1.1
        )
        energy_axis.set_ylabel("Mean energy" if scale_index == 0 else "")
        energy_axis.set_xlabel(f"Position from {boundary} boundary (bp)")
        energy_axis.grid(axis="y", color="#888888", alpha=0.2, linewidth=0.6)
        for axis in (image_axis, energy_axis):
            axis.set_xticks(
                (-ABS_RADIUS, 0, ABS_RADIUS),
                (f"-{ABS_RADIUS}", "0", f"+{ABS_RADIUS}"),
                fontsize=8,
            )
            axis.set_xlim(extent[:2])

    output.mkdir(parents=True, exist_ok=True)
    figure.savefig(output / f"{boundary}.png", dpi=160)
    plt.close(figure)


class CwtPlotter:
    @classmethod
    def run_meme(cls, fasta: Path, output: Path) -> list[dict[str, Any]]:
        output.mkdir(parents=True, exist_ok=True)
        alphabet_path = output.parent / "dna-single.alph"
        if not alphabet_path.is_file():
            alphabet_path.write_text(SINGLE_STRAND_DNA_ALPH)

        environment = {**os.environ, "TMPDIR": str(output.resolve())}
        command = [
            "streme",
            "--verbosity",
            "1",
            "--oc",
            str(output.resolve()),
            "--alph",
            str(alphabet_path.resolve()),
            "--seed",
            "0",
            "--p",
            str(fasta.resolve()),
            "--minw",
            str(MEME_MINW),
            "--maxw",
            str(MEME_MAXW),
            "--time",
            str(MEME_TIME),
            "--order",
            str(MEME_ORDER),
            "--nmotifs",
            str(MEME_MOTIFS),
            "--align",
            "center",
            "--no-pgc",
        ]
        subprocess.run(command, check=True, env=environment)

        motifs: list[dict[str, Any]] = []
        html_path = output / "streme.html"
        if html_path.is_file():
            html = html_path.read_text()
            marker = "// @JSON_VAR data"
            if marker in html:
                marker_index = html.index(marker)
                json_start = html.index("=", marker_index) + 1
                json_end = html.index(";\n    </script>", marker_index)
                data = json.loads(html[json_start:json_end].strip())
                alphabet = data.get("alphabet", {}).get("symbols", [])
                letters = [
                    symbol.get("symbol", "")
                    for symbol in alphabet
                    if not symbol.get("equals")
                ]
                base_indices = (
                    [letters.index(base) for base in "ACGT"]
                    if set("ACGT").issubset(letters)
                    else [0, 1, 2, 3]
                )

                for motif in data.get("motifs", []):
                    probabilities = np.asarray(motif["pwm"], dtype=np.float64)[
                        :, base_indices
                    ]
                    motifs.append(
                        {
                            "label": motif.get("alt") or "STREME",
                            "name": motif.get("id", "").split("-", 1)[-1]
                            or motif.get("seed", ""),
                            "width": int(motif.get("width", len(probabilities))),
                            "sites": int(
                                motif.get("npassing", motif.get("total_sites", 0))
                            ),
                            "e_value": motif.get("test_evalue")
                            or motif.get("train_pvalue")
                            or "n/a",
                            "probabilities": probabilities,
                        }
                    )
        return motifs

    @staticmethod
    def render_meme_logo(axis: Any, motif: dict[str, Any], title: str) -> None:
        probabilities = motif["probabilities"]
        positions = np.arange(1, len(probabilities) + 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            entropy = -np.sum(
                np.where(
                    probabilities > 0,
                    probabilities * np.log2(probabilities),
                    0.0,
                ),
                axis=1,
                keepdims=True,
            )
        heights = probabilities * np.clip(2.0 - entropy, 0.0, 2.0)
        colors = ("#2ca25f", "#3182bd", "#f28e2b", "#d73027")
        font = FontProperties(weight="bold")
        glyphs = {base: TextPath((0, 0), base, size=1, prop=font) for base in "ACGT"}
        glyph_bounds = {base: glyph.get_extents() for base, glyph in glyphs.items()}

        for column, position in enumerate(positions):
            bottom = 0.0
            for base_index in np.argsort(heights[column]):
                height = float(heights[column, base_index])
                if not height > 0:
                    continue

                base = "ACGT"[base_index]
                bounds = glyph_bounds[base]
                transform = (
                    Affine2D()
                    .translate(-bounds.x0, -bounds.y0)
                    .scale(0.82 / bounds.width, height / bounds.height)
                    .translate(position - 0.41, bottom)
                    + axis.transData
                )
                axis.add_patch(
                    PathPatch(
                        glyphs[base],
                        transform=transform,
                        color=colors[base_index],
                    )
                )
                bottom += height

        ticks = positions[:: max(1, len(positions) // 8)]
        if ticks[-1] != positions[-1]:
            ticks = np.append(ticks, positions[-1])
        axis.set(
            ylabel="Information (bits)",
            ylim=(0, 2.04),
            yticks=(0, 1, 2),
            xlim=(0.5, len(positions) + 0.5),
            xlabel="Motif position",
        )
        axis.set_xticks(ticks)
        axis.grid(axis="y", color="#888888", alpha=0.2, linewidth=0.6)
        axis.set_axisbelow(True)
        axis.set_title(title, fontsize=9)

    @classmethod
    def render_exon_kmer(
        cls,
        windows: dict[str, list[KmerHit]],
        peaks: dict[str, int],
        skipped: int,
        output: Path,
    ) -> None:
        output.mkdir(parents=True, exist_ok=True)
        fasta_paths: dict[str, Path] = {}
        for boundary in BOUNDARIES:
            boundary_windows = windows[boundary]
            if not boundary_windows:
                raise ValueError(
                    f"No valid {KMER_WIDTH}-mer CWT maxima for exon {boundary}"
                )

            fasta_path = output / f"{boundary}-maxima-{KMER_WIDTH}mer.fasta"
            with fasta_path.open("w") as stream:
                stream.writelines(
                    f">{header}\n{sequence}\n" for header, sequence in boundary_windows
                )
            fasta_paths[boundary] = fasta_path
            print(
                f"Running STREME on {len(boundary_windows):,} exon {boundary} "
                f"maxima at {peaks[boundary]:+d} bp",
                file=sys.stderr,
                flush=True,
            )

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = {
                boundary: pool.submit(
                    cls.run_meme,
                    fasta_paths[boundary],
                    output / f"meme-{boundary}",
                )
                for boundary in BOUNDARIES
            }
            results = {boundary: futures[boundary].result() for boundary in BOUNDARIES}

        if skipped:
            print(
                f"Skipped {skipped:,} maxima with incomplete or ambiguous "
                f"{KMER_WIDTH}-mers",
                file=sys.stderr,
                flush=True,
            )
        with (output / "maxima-distances.tsv").open("w") as stream:
            stream.write("flank\tsequence_id\tdistance_bp\n")
            for boundary in BOUNDARIES:
                stream.writelines(
                    f"{boundary}\t{header}\t{peaks[boundary]}\n"
                    for header, _ in windows[boundary]
                )

        figure, axes = plt.subplots(
            MEME_MOTIFS,
            2,
            figsize=(17, 4.2 * MEME_MOTIFS),
            layout="constrained",
            squeeze=False,
        )
        titles = {
            "start": (
                f"Exon start upstream (peak {peaks['start']:+d} bp, "
                f"n={len(windows['start']):,})"
            ),
            "end": (
                f"Exon end downstream (peak {peaks['end']:+d} bp, "
                f"n={len(windows['end']):,})"
            ),
        }
        for column, boundary in enumerate(BOUNDARIES):
            for row in range(MEME_MOTIFS):
                axis = axes[row, column]
                if row < len(results[boundary]):
                    motif = results[boundary][row]
                    motif_title = (
                        f"{titles[boundary]}\n"
                        f"{motif['label']}: {motif['name']}\n"
                        f"width={motif['width']}, sites={motif['sites']:,}, "
                        f"E={motif['e_value']}"
                    )
                    cls.render_meme_logo(axis, motif, motif_title)
                else:
                    axis.axis("off")
                    axis.text(
                        0.5,
                        0.5,
                        f"{titles[boundary]}\nNo additional motif",
                        ha="center",
                        va="center",
                        transform=axis.transAxes,
                    )
        figure.savefig(output / "motifs.png", dpi=180)
        plt.close(figure)

    @staticmethod
    def resample_mean(values: Arr, bins: int) -> Arr:
        array = np.asarray(values, dtype=np.float64)
        if array.ndim != 2 or not len(array) or bins < 1:
            raise ValueError(
                "Resampling requires nonempty position-by-scale values "
                "and positive bins"
            )

        edges = np.linspace(0.0, len(array), bins + 1)
        cumulative_integral = np.concatenate(
            (
                np.zeros((1, array.shape[1])),
                np.cumsum(array, axis=0),
            ),
            axis=0,
        )
        lower_indices = np.floor(edges).astype(np.int64)
        interpolated_integral = cumulative_integral[lower_indices] + (
            (edges - lower_indices)[:, None]
            * array[np.minimum(lower_indices, len(array) - 1)]
        )
        return np.diff(interpolated_integral, axis=0) / (len(array) / bins)

    @classmethod
    def render_roi(
        cls,
        intervals: Sequence[Interval],
        cwt: CwtMatrix,
        output: Path,
        flank: int,
        pixels: int,
    ) -> None:
        output.mkdir(parents=True, exist_ok=True)
        for name, start, end, strand in intervals:
            contig = cwt.contigs[name]
            first, last = contig.bounds(start, end, strand)
            upstream_flank = min(flank, first)
            downstream_flank = min(flank, contig.length - last)
            start_index = contig.offset(strand) + first - upstream_flank
            end_index = contig.offset(strand) + last + downstream_flank
            values = cwt.matrix[start_index:end_index]
            if len(values) > 100000:
                raise ValueError(
                    "Use a region of at most 100,000 bp (including flanks)"
                )

            bin_count = min(pixels, len(values))
            body_end = upstream_flank + last - first
            panel = np.log1p(
                cls.resample_mean(cwt.power_at(start_index, end_index), bin_count)
            ).T
            figure, axis = plt.subplots(figsize=(13, 4), layout="constrained")
            maximum = max(float(np.nanmax(panel)), 1e-6)
            image = axis.imshow(
                panel,
                aspect="auto",
                origin="lower",
                interpolation="nearest",
                extent=(0, len(values), -0.5, len(cwt.widths) - 0.5),
                cmap="viridis",
                vmin=0,
                vmax=maximum,
            )
            figure.colorbar(image, ax=axis, pad=0.01)
            axis.set_yticks(
                range(len(cwt.widths)),
                [str(width) for width in cwt.widths],
            )
            axis.set_ylabel("Surf CWT energy\n(mean over periods)")
            axis.set_title("log(1 + mean Surf CWT energy)", loc="left", fontsize=11)

            if strand == "+":
                left_boundary = ("start", start)
                right_boundary = ("end", end)
            else:
                left_boundary = ("end", end)
                right_boundary = ("start", start)

            ticks: list[tuple[int, str, Literal["left", "right"]]] = []
            if upstream_flank:
                ticks.append((0, f"5' flank\n{upstream_flank} bp", "left"))
            ticks.extend(
                (
                    (
                        upstream_flank,
                        f"{left_boundary[0]}\n{left_boundary[1]:,}",
                        "right",
                    ),
                    (
                        body_end,
                        f"{right_boundary[0]}\n{right_boundary[1]:,}",
                        "left",
                    ),
                )
            )
            if downstream_flank:
                ticks.append((len(values), f"3' flank\n{downstream_flank} bp", "right"))
            axis.set_xticks(
                [position for position, _, _ in ticks],
                [label for _, label, _ in ticks],
            )
            for label, (_, _, alignment) in zip(axis.get_xticklabels(), ticks):
                label.set_horizontalalignment(alignment)
            axis.tick_params(axis="x", labelsize=8, pad=6)
            axis.set_xlabel("5' to 3' position; start/end mark genomic exon boundaries")
            figure.suptitle(f"{name}:{start:,}-{end:,} ({strand}) | Surf CWT")
            strand_name = "plus" if strand == "+" else "minus"
            figure.savefig(
                output / f"{name}_{start:06d}_{end:06d}_{strand_name}.png",
                dpi=160,
            )
            plt.close(figure)

    @classmethod
    def render_group(
        cls,
        names: list[str],
        cwt: CwtMatrix,
        genome: GenomeData,
        fasta: Path,
        args: argparse.Namespace,
    ) -> list[AbsRenderTask]:
        repeat_free_names = [name for name in names if "no-repeat" in name]
        if repeat_free_names and not genome.has_repeats:
            print(
                "Skipping repeat-free plots: no repeat annotations on selected contigs",
                file=sys.stderr,
                flush=True,
            )
            names = [name for name in names if name not in repeat_free_names]
        if not names:
            return []

        sequence_dependent_plots = ("exon-kmer", "control-aggt", "control-gtag")
        if any(
            prefix in plot_name
            for plot_name in names
            for prefix in sequence_dependent_plots
        ):
            genome.load_sequences(fasta)

        exons = genome.exon_sets["exons"]
        exon_aggregate = (
            cwt.aggregate_boundaries(exons, args.rows) if "exons" in names else None
        )
        kmer_cache: tuple[dict[str, list[KmerHit]], dict[str, int], int] | None = None
        if "exon-kmer" in names:
            multi_exons = genome.exon_sets["multi-exons"]
            if not multi_exons:
                raise ValueError("No multi-exon features found for exon-kmer")
            multi_exon_aggregate = cwt.aggregate_boundaries(multi_exons, args.rows)
            kmer_cache = cwt.kmer_maxima(
                multi_exons,
                genome.sequences,
                multi_exon_aggregate[1],
                genome.kmer_bounds,
            )

        absolute_tasks: list[AbsRenderTask] = []

        def queue_absolute_plot(
            aggregate: tuple[Arr, Arr],
            interval_count: int,
            output: Path,
            boundary_indices: Sequence[int] = range(len(BOUNDARIES)),
        ) -> None:
            log_image = np.log1p(aggregate[0])
            absolute_tasks.extend(
                (
                    boundary,
                    np.ascontiguousarray(log_image[:, index]),
                    np.ascontiguousarray(aggregate[1][index]),
                    cwt.widths,
                    interval_count,
                    output,
                )
                for index, boundary in enumerate(BOUNDARIES)
                if index in boundary_indices
            )

        for plot_name in names:
            output = args.out / plot_name
            if plot_name == "exon-kmer":
                assert kmer_cache is not None
                cls.render_exon_kmer(*kmer_cache, output)
                continue

            if plot_name == "first-last-exons":
                first_exons = genome.exon_sets["first-last-start"]
                last_exons = genome.exon_sets["first-last-end"]
                if not first_exons or not last_exons:
                    raise ValueError(
                        "No multi-exon features found for first-last-exons"
                    )
                queue_absolute_plot(
                    cwt.aggregate_boundaries(first_exons, args.rows),
                    len(first_exons),
                    output,
                    (0,),
                )
                queue_absolute_plot(
                    cwt.aggregate_boundaries(last_exons, args.rows),
                    len(last_exons),
                    output,
                    (1,),
                )
                continue

            group = plot_name.removesuffix("-mutated")
            base_group = group.removesuffix("-no-repeat")
            repeat_free = group.endswith("-no-repeat")
            is_aggt_control = plot_name.startswith("control-aggt")
            is_gtag_control = plot_name.startswith("control-gtag")
            replicates = (
                range(1, CONTROL_REPLICATES + 1)
                if is_aggt_control or is_gtag_control
                else (0,)
            )

            for replicate in replicates:
                replicate_output = (
                    output / f"replicate-{replicate}"
                    if is_aggt_control or is_gtag_control
                    else output
                )
                length_offset = int(base_group.endswith("-p1"))
                if is_gtag_control:
                    selected = genome.sample_controls(
                        args.seed + replicate - 1,
                        gtag=True,
                        repeat_free=repeat_free,
                    )
                elif is_aggt_control:
                    selected = genome.sample_controls(
                        args.seed + length_offset * CONTROL_REPLICATES + replicate - 1,
                        aggt=True,
                        length_offset=length_offset,
                        repeat_free=repeat_free,
                    )
                elif base_group == "control":
                    selected = genome.sample_controls(
                        args.seed, repeat_free=repeat_free
                    )
                else:
                    selected = genome.exon_sets[group]

                if not selected:
                    raise ValueError(f"No {group} features found")
                if repeat_free or is_gtag_control:
                    total = len(
                        genome.exon_sets["introns"] if is_gtag_control else exons
                    )
                    print(
                        f"Sampled {len(selected):,}/{total:,} {group} "
                        f"{replicate_output.name}; unmatched lengths skipped "
                        "after 200 attempts",
                        file=sys.stderr,
                        flush=True,
                    )

                aggregate = (
                    exon_aggregate
                    if group == "exons" and exon_aggregate is not None
                    else cwt.aggregate_boundaries(selected, args.rows)
                )
                queue_absolute_plot(
                    aggregate,
                    len(selected),
                    replicate_output,
                )

        return absolute_tasks

    @classmethod
    def run(cls, args: argparse.Namespace) -> None:
        args.out.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(
            prefix=".cwt-work-",
            dir=args.out,
            ignore_cleanup_errors=True,
        ) as temporary_directory:
            work = Path(temporary_directory)

            if args.regions:
                if args.gff:
                    raise ValueError("-r/--region cannot be combined with -g/--gff")
                cwt = CwtMatrix.load(
                    args.fasta,
                    args.cwt or (work / "original"),
                    args.threads,
                    args.chunk,
                )
                try:
                    intervals: list[Interval] = []
                    for region_spec in args.regions:
                        bed_path = Path(region_spec)
                        if bed_path.is_file():
                            with open_text(bed_path) as stream:
                                for line_number, line in enumerate(stream, 1):
                                    stripped = line.rstrip("\r\n")
                                    if not stripped or stripped.startswith(
                                        ("#", "track", "browser")
                                    ):
                                        continue

                                    fields = stripped.split("\t")
                                    if len(fields) < 3:
                                        raise ValueError(
                                            "Expected at least 3 BED columns "
                                            f"at line {line_number}"
                                        )
                                    strand = (
                                        fields[5]
                                        if len(fields) > 5 and fields[5] in ("+", "-")
                                        else "+"
                                    )
                                    intervals.append(
                                        cwt.validate(
                                            fields[0],
                                            int(fields[1]) + 1,
                                            int(fields[2]),
                                            strand,
                                            "BED contig",
                                        )
                                    )
                            continue

                        region_fields = region_spec.split(":")
                        if (
                            len(region_fields) not in (2, 3)
                            or "-" not in region_fields[1]
                        ):
                            raise ValueError(
                                f"Invalid region spec {region_spec!r}; expected "
                                "SEQID:START-END[:STRAND] or BED file"
                            )
                        start_text, end_text = region_fields[1].split("-", 1)
                        strand = region_fields[2] if len(region_fields) == 3 else "+"
                        intervals.append(
                            cwt.validate(
                                region_fields[0],
                                int(start_text.replace(",", "")),
                                int(end_text.replace(",", "")),
                                strand,
                            )
                        )

                    if not intervals:
                        raise ValueError("No regions found")
                    cls.render_roi(
                        intervals,
                        cwt,
                        args.out,
                        args.flank,
                        args.pixels,
                    )
                finally:
                    cwt.close()
                return

            if not args.gff:
                raise ValueError("-P/--plots requires -g/--gff")

            requested_plots = set(args.plots)
            selected_plots = (
                list(PLOT_NAMES)
                if "all" in requested_plots
                else [name for name in PLOT_NAMES if name in requested_plots]
            )
            original_plots = [
                name for name in selected_plots if not name.endswith("-mutated")
            ]
            mutated_plots = [
                name for name in selected_plots if name.endswith("-mutated")
            ]
            tasks: list[AbsRenderTask] = []
            genome: GenomeData | None = None

            if original_plots:
                cwt = CwtMatrix.load(
                    args.fasta,
                    args.cwt or (work / "original"),
                    args.threads,
                    args.chunk,
                )
                try:
                    genome = GenomeData(args.gff, cwt.contigs)
                    print(
                        f"Plotting {', '.join(original_plots)}",
                        file=sys.stderr,
                        flush=True,
                    )
                    tasks.extend(
                        cls.render_group(original_plots, cwt, genome, args.fasta, args)
                    )
                finally:
                    cwt.close()
                del cwt

            if mutated_plots:
                if genome is None:
                    sequences = {
                        name: sequence.upper()
                        for name, sequence in fasta_records(args.fasta)
                    }
                    if not sequences:
                        raise ValueError("No sequences found in assembly FASTA")
                    contigs = {
                        name: Contig(len(sequence), 0, len(sequence))
                        for name, sequence in sequences.items()
                    }
                    genome = GenomeData(args.gff, contigs, sequences)

                mutated_fasta = genome.write_mutated_fasta(
                    args.fasta,
                    args.seed,
                    work / "mutated.fa",
                )
                print(
                    f"Plotting {', '.join(mutated_plots)}",
                    file=sys.stderr,
                    flush=True,
                )
                mutated_cwt = CwtMatrix.load(
                    mutated_fasta,
                    work / "mutated",
                    args.threads,
                    args.chunk,
                )
                try:
                    tasks.extend(
                        cls.render_group(
                            mutated_plots,
                            mutated_cwt,
                            genome,
                            mutated_fasta,
                            args,
                        )
                    )
                finally:
                    mutated_cwt.close()
                del mutated_cwt
            del genome

        if tasks:
            worker_count = min(len(tasks), args.threads)
            if worker_count > 1:
                with ProcessPoolExecutor(max_workers=worker_count) as pool:
                    list(pool.map(_render_abs_task, tasks))
            else:
                for task in tasks:
                    _render_abs_task(task)


def int_min(minimum: int, label: str) -> Callable[[str], int]:
    def parse(text: str) -> int:
        value = int(text)
        if value < minimum:
            raise argparse.ArgumentTypeError(f"Expected a {label} integer")
        return value

    return parse


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        epilog=f"Available plots: {', '.join(PLOT_CHOICES)}"
    )
    positive_integer = int_min(1, "positive")
    nonnegative_integer = int_min(0, "nonnegative")

    parser.add_argument("fasta", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-P",
        "--plots",
        nargs="+",
        choices=PLOT_CHOICES,
        metavar="PLOT",
    )
    mode.add_argument(
        "-r",
        "--region",
        dest="regions",
        nargs="+",
        action="extend",
        metavar="REGION",
    )
    parser.add_argument("-g", "--gff", type=Path)
    parser.add_argument("-o", "--out", type=Path, default=Path("plots"))
    parser.add_argument("-c", "--cwt", type=Path)
    parser.add_argument(
        "--chunk",
        type=positive_integer,
        default=65536,
        help="Bases per CWT energy chunk",
    )
    parser.add_argument(
        "-p",
        "--threads",
        type=positive_integer,
        default=max(1, os.cpu_count() or 1),
    )
    parser.add_argument("-n", "--rows", type=positive_integer, default=1500)
    parser.add_argument("-s", "--seed", type=int, default=42)
    parser.add_argument("-f", "--flank", type=nonnegative_integer, default=10)
    parser.add_argument("--pixels", type=positive_integer, default=1600)
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        CwtPlotter.run(args)
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        json.JSONDecodeError,
        subprocess.CalledProcessError,
    ) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()
