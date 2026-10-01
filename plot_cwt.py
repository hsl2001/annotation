import argparse
import csv
import gzip
import pathlib
import subprocess
import sys
import tempfile
from collections import defaultdict
from itertools import chain

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

PLOT_NAMES = ("exons", "exons-mutated", "control", "control-aggt",
              "abs-exons", "abs-exons-mutated", "abs-control", "abs-control-aggt")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def fasta_records(path, wanted=None):
    name, chunks = None, []
    with open_text(path) as stream:
        for line in stream:
            if line.startswith(">"):
                if name is not None and (wanted is None or name in wanted):
                    yield name, "".join(chunks)
                header = line[1:].split()
                name = header[0] if header else ""
                chunks = []
            elif name is not None and (wanted is None or name in wanted):
                chunks.append(line.strip())
    if name is not None and (wanted is None or name in wanted):
        yield name, "".join(chunks)


def gff_exons(path):
    with open_text(path) as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE), 1):
            if not row or row[0].startswith("#"):
                if row and row[0] == "##FASTA":
                    break
                continue
            if len(row) != 9:
                raise ValueError(f"Expected 9 GFF/GTF columns at line {line_number}")
            if row[2] == "exon":
                yield row


def load_matrix(directory):
    with open(directory / "contigs.tsv") as stream:
        version = stream.readline().strip()
        if version not in ("# anno-cwt-v1", "# anno-cwt-v2"):
            raise ValueError("Unsupported CWT index version")
        fields = stream.readline().rstrip().split("\t")
        parameter = "Scale" if version == "# anno-cwt-v2" else "Kernel width"
        expected_field = "# scales" if version == "# anno-cwt-v2" else "# kernel_widths"
        if fields[0] != expected_field:
            raise ValueError("Missing CWT parameters")
        widths = [int(value) for value in ",".join(fields[1:]).split(",") if value]
        contigs = {}
        expected = 0
        for row in csv.DictReader(stream, delimiter="\t"):
            name = row["seqid"]
            length, plus, minus = (int(row[key]) for key in ("length", "plus_offset", "minus_offset"))
            if name in contigs or length <= 0 or plus != expected or minus != plus + length:
                raise ValueError(f"Invalid contig offsets: {name}")
            contigs[name] = {"length": length, "+": plus, "-": minus, "parameter": parameter}
            expected += 2 * length
    if not widths or not contigs:
        raise ValueError("Missing CWT metadata")
    path = directory / "matrix.bin"
    expected_bytes = expected * len(widths) * np.dtype("<c8").itemsize
    if not path.is_file() or path.stat().st_size != expected_bytes:
        raise ValueError("Expected a little-endian complex64 matrix matching contigs.tsv")
    matrix = np.memmap(path, dtype="<c8", mode="r", shape=(expected, len(widths)))
    return matrix, contigs, widths


def run_anno_cwt(binary, fasta, output):
    executable = str(binary.resolve()) if binary.exists() else str(binary)
    matrix_dir = output / "cwt"
    subprocess.run([executable, str(fasta), str(matrix_dir)], check=True)
    return matrix_dir


def oriented_bounds(contig, start, end, strand):
    if strand not in ("+", "-") or start < 1 or end < start:
        raise ValueError(f"Invalid interval: {start}-{end} ({strand})")
    length = contig["length"]
    span = end - start + 1
    if end > length:
        if span > length:
            raise ValueError(f"Circular interval exceeds contig length: {start}-{end} ({strand})")
        start = (start - 1) % length + 1
        end = start + span - 1
        if end > length:
            raise ValueError(f"Circular interval crosses contig origin: {start}-{end} ({strand})")
    first = start - 1 if strand == "+" else contig["length"] - end
    return first, first + end - start + 1


def power(values):
    return np.square(values.real, dtype=np.float64) + np.square(values.imag, dtype=np.float64)


def resample_mean(values, bins):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2 or not len(values) or bins < 1:
        raise ValueError("Resampling requires nonempty position-by-scale values and positive bins")
    edges = np.linspace(0.0, len(values), bins + 1)
    integral = np.concatenate((np.zeros((1, values.shape[1])), np.cumsum(values, axis=0)), axis=0)
    lower = np.floor(edges).astype(np.int64)
    cumulative = integral[lower] + (edges - lower)[:, None] * values[np.minimum(lower, len(values) - 1)]
    return np.diff(cumulative, axis=0) / (len(values) / bins)


def exon_power(matrix, contig, start, end, strand, bins, flank):
    first, last = oriented_bounds(contig, start, end, strand)
    offset = contig[strand]
    body = power(matrix[offset + first:offset + last])
    result = np.full((bins + 2 * flank, matrix.shape[1]), np.nan, dtype=np.float32)
    result[flank:flank + bins] = resample_mean(body, bins)
    before, after = min(flank, first), min(flank, contig["length"] - last)
    if before:
        result[flank - before:flank] = power(matrix[offset + first - before:offset + first])
    if after:
        result[flank + bins:flank + bins + after] = power(matrix[offset + last:offset + last + after])
    if not np.isfinite(body).all():
        raise ValueError("Nonfinite CWT coefficients in exon")
    return result.T


def read_exons(path, contigs):
    intervals = set()
    for row in gff_exons(path):
        name, start, end, strand = row[0], int(row[3]), int(row[4]), row[6]
        if name not in contigs:
            raise ValueError(f"Exon contig absent from CWT matrix: {name}")
        oriented_bounds(contigs[name], start, end, strand)
        intervals.add((name, start, end, strand))
    if not intervals:
        raise ValueError("No exon features found; CDS is not silently substituted for exon")
    return sorted(intervals, key=lambda item: (item[2] - item[1] + 1, item))


def aggregate_intervals(intervals, rows, sample_power, stem, action):
    intervals = sorted(intervals, key=lambda item: (item[2] - item[1] + 1, item))
    display_rows = min(rows, len(intervals))
    edges = np.linspace(0, len(intervals), display_rows + 1).astype(np.int64)
    samples = (sample_power(interval) for interval in intervals)
    first = next(samples)
    display_total = np.zeros((display_rows, *first.shape), dtype=np.float64)
    display_counts = np.zeros_like(display_total, dtype=np.int64)
    total = np.zeros(first.shape, dtype=np.float64)
    counts = np.zeros_like(total, dtype=np.int64)
    display_row = 0
    for row, values in enumerate(chain((first,), samples)):
        while display_row + 1 < display_rows and row >= edges[display_row + 1]:
            display_row += 1
        valid = np.isfinite(values)
        values = np.where(valid, values, 0.0)
        total += values
        counts += valid
        display_total[display_row] += values
        display_counts[display_row] += valid
        if (row + 1) % 10000 == 0:
            print(f"{action} {row + 1:,}/{len(intervals):,} {stem}", file=sys.stderr, flush=True)
    display = np.divide(display_total, display_counts, out=np.full_like(display_total, np.nan),
                        where=display_counts > 0).astype(np.float32)
    mean = np.divide(total, counts, out=np.full_like(total, np.nan), where=counts > 0)
    return display, mean


def heatmap_limit(values):
    return max(float(np.quantile(values[np.isfinite(values)], 0.995)), 1e-6)


def render_power_panels(display, mean, contigs, widths, path, *, positions, extent, ticks,
                        xlabel, height, maximum=None, color_label):
    image_values = np.log1p(display)
    if maximum is None:
        maximum = heatmap_limit(image_values)
    parameter = next(iter(contigs.values())).get("parameter", "Kernel width")
    figure, axes = plt.subplots(2, len(widths), figsize=(4.1 * len(widths), height),
                                gridspec_kw={"height_ratios": [4, 1.3]}, layout="constrained", squeeze=False)
    for scale, width in enumerate(widths):
        axis, profile = axes[:, scale]
        image = axis.imshow(image_values[:, scale], aspect="auto", origin="lower", interpolation="nearest",
                            extent=extent, cmap="viridis", vmin=0, vmax=maximum)
        axis.set_title(f"{parameter} {width} bp")
        axis.set_ylabel("Interval rank" if scale == 0 else "")
        profile.plot(positions, mean[scale], color="#256c87", linewidth=1.3)
        profile.set_ylabel("Mean power" if scale == 0 else "")
        profile.set_ylim(-0.1, 2.1)
        for panel in (axis, profile):
            panel.set_xticks(*ticks, fontsize=8)
            panel.set_xlim(extent[:2])
        profile.set_xlabel(xlabel)
    figure.colorbar(image, ax=list(axes[0]), label=color_label, shrink=0.65, pad=0.01)
    figure.savefig(path, dpi=160)
    plt.close(figure)


def render_intervals(intervals, matrix, contigs, widths, out, bins, flank, rows, stem):
    display, mean = aggregate_intervals(
        intervals, rows,
        lambda interval: exon_power(matrix, contigs[interval[0]], interval[1], interval[2],
                                    interval[3], bins, flank),
        stem, "Normalized")
    columns = bins + 2 * flank
    ticks, labels = [flank, flank + bins / 2, flank + bins], ["0%", "50%", "100%"]
    if flank:
        ticks, labels = [0, *ticks, columns], [f"-{flank} bp", *labels, f"+{flank} bp"]
    render_power_panels(
        display, mean, contigs, widths, out / f"{stem}.png",
        positions=np.arange(columns) + 0.5, extent=(0, columns, 0, len(intervals)),
        ticks=(ticks, labels), xlabel="5' to 3': normalized body + genomic flanks", height=9,
        color_label="log(1 + mean power), common scale; top 0.5% clipped")


def boundary_power(matrix, contig, interval, boundary, radius=100):
    name, start, end, strand = interval
    first, last = oriented_bounds(contig, start, end, strand)
    center = first if boundary == "start" else last - 1
    left, right = center - radius, center + radius + 1
    window = np.full((2 * radius + 1, matrix.shape[1]), np.nan, dtype=np.float32)
    clipped_left, clipped_right = max(0, left), min(contig["length"], right)
    if clipped_left < clipped_right:
        destination = clipped_left - left
        window[destination:destination + clipped_right - clipped_left] = power(
            matrix[contig[strand] + clipped_left:contig[strand] + clipped_right])
    return window


def render_absolute_intervals(intervals, matrix, contigs, widths, out, rows, stem):
    display, profiles = aggregate_intervals(
        intervals, rows,
        lambda interval: np.stack([boundary_power(matrix, contigs[interval[0]], interval, boundary).T
                                   for boundary in ("start", "end")]),
        stem, "Aligned")
    maximum = heatmap_limit(np.log1p(display))
    for boundary_index, boundary in enumerate(("start", "end")):
        render_power_panels(
            display[:, boundary_index], profiles[boundary_index], contigs, widths, out / f"{boundary}.png",
            positions=np.arange(-100, 101), extent=(-100.5, 100.5, 0, len(intervals)),
            ticks=((-100, 0, 100), ("-100", "0", "+100")),
            xlabel=f"Position from {boundary} boundary (bp)", height=6, maximum=maximum,
            color_label="log(1 + power), common scale; top 0.5% clipped")


def read_genome(path):
    """Load the assembly as ordered, mutable sequences keyed by contig name."""
    order, sequences = [], {}
    for name, sequence in fasta_records(path):
        if name in sequences:
            raise ValueError(f"Duplicate FASTA contig: {name}")
        order.append(name)
        sequences[name] = bytearray(sequence, "ascii")
    if not order:
        raise ValueError("No sequences found in assembly FASTA")
    return order, sequences


def exon_parents(attributes_text):
    for field in attributes_text.split(";"):
        key, separator, value = field.partition("=")
        if separator and key.strip() == "Parent":
            return [parent for parent in value.split(",") if parent]
    return []


def splice_positions(path, contigs):
    """Return intron-boundary and transcript-terminal exon positions for mutation."""
    transcripts = defaultdict(list)
    for row in gff_exons(path):
        name, start, end = row[0], int(row[3]), int(row[4])
        if name not in contigs:
            raise ValueError(f"Exon contig absent from assembly FASTA: {name}")
        for parent in exon_parents(row[8]):
            transcripts[name, parent].append((start, end))
    positions = defaultdict(set)
    introns = 0
    for (name, _), exons in transcripts.items():
        exons.sort()
        length = contigs[name]
        for site in (exons[0][0], exons[0][0] + 1, exons[-1][1] - 1, exons[-1][1]):
            if 1 <= site <= length:
                positions[name].add(site - 1)
        for (_, prev_end), (next_start, _) in zip(exons, exons[1:]):
            low, high = prev_end + 1, next_start - 1
            if high < low:
                continue
            introns += 1
            for site in (low, low + 1, high - 1, high):
                if low <= site <= high and 1 <= site <= length:
                    positions[name].add(site - 1)
    if not introns:
        raise ValueError("No introns found; splice donors/acceptors require multi-exon transcripts")
    return positions


def mutate_splice_sites(sequences, positions, seed):
    rng = np.random.default_rng(seed)
    alphabet = "ACGT"
    for name, sites in positions.items():
        sequence = sequences[name]
        for index in sorted(sites):
            current = chr(sequence[index]).upper()
            options = alphabet.replace(current, "") if current in alphabet else alphabet
            sequence[index] = ord(options[int(rng.integers(len(options)))])


def write_genome(path, order, sequences, width=60):
    with open(path, "w") as stream:
        for name in order:
            stream.write(f">{name}\n")
            sequence = sequences[name]
            for offset in range(0, len(sequence), width):
                stream.write(sequence[offset:offset + width].decode("ascii"))
                stream.write("\n")


def make_mutated_fasta(fasta, gff, seed, output):
    order, sequences = read_genome(fasta)
    contigs = {name: len(sequence) for name, sequence in sequences.items()}
    positions = splice_positions(gff, contigs)
    mutate_splice_sites(sequences, positions, seed)
    write_genome(output, order, sequences)
    return output


def plot_exon_length_distribution(exons, out):
    lengths = np.array([end - start + 1 for _, start, end, _ in exons])
    minimum, maximum = int(lengths.min()), int(lengths.max())
    edges = np.array([minimum * 0.5, maximum * 1.5]) if minimum == maximum else np.geomspace(minimum, maximum, 46)
    counts, edges = np.histogram(lengths, bins=edges)
    figure, axis = plt.subplots(figsize=(9, 4.8), layout="constrained")
    axis.stairs(counts, edges, fill=True, color="#287c8e", alpha=0.82, linewidth=1.0)
    axis.set_xscale("log")
    axis.set_xlim(minimum, maximum)
    ticks = [tick for tick in (1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000) if minimum <= tick <= maximum]
    axis.set_xticks(ticks)
    axis.xaxis.set_major_formatter(ScalarFormatter())
    axis.set_xlabel("Exon length (bp)")
    axis.set_ylabel("Unique exon count")
    axis.grid(axis="y", color="#d8dedf", linewidth=0.7)
    figure.savefig(out / "exon_length_distribution.png", dpi=180)
    plt.close(figure)


def build_exon_mask(contigs, exons):
    masks = {name: np.zeros(info["length"], dtype=bool) for name, info in contigs.items()}
    for name, start, end, _ in exons:
        masks[name][start - 1:end] = True
    return masks


def sample_controls(contigs, exons, seed, attempts=200):
    masks = build_exon_mask(contigs, exons)
    names = list(contigs)
    weights = np.cumsum([contigs[name]["length"] for name in names], dtype=np.float64)
    rng = np.random.default_rng(seed)
    controls = []
    for _, start0, end0, _ in exons:
        need = end0 - start0 + 1
        for _ in range(attempts):
            index = min(int(np.searchsorted(weights, rng.random() * weights[-1], side="right")), len(names) - 1)
            name = names[index]
            span = contigs[name]["length"]
            if span < need:
                continue
            start = int(rng.integers(1, span - need + 2))
            if not masks[name][start - 1:start - 1 + need].any():
                controls.append((name, start, start + need - 1, "+" if rng.random() < 0.5 else "-"))
                break
    if not controls:
        raise ValueError("Could not place any non-exon control at the requested lengths")
    return controls


def read_sequences(path, contigs):
    wanted = set(contigs)
    sequences = {name: sequence.upper() for name, sequence in fasta_records(path, wanted)}
    missing = wanted - sequences.keys()
    if missing:
        raise ValueError(f"FASTA missing sequences for CWT contigs: {', '.join(sorted(missing))}")
    for name, sequence in sequences.items():
        if len(sequence) != contigs[name]["length"]:
            raise ValueError(f"FASTA length for {name} ({len(sequence)}) differs from CWT matrix "
                             f"({contigs[name]['length']})")
    return sequences


def aggt_candidates(sequences):
    """Precompute, per strand, non-exon start positions carrying the transcript-oriented AG acceptor,
    plus the donor dimer mask used to test the exon end. In transcript 5'->3' orientation an internal
    exon reads ...AG | exon | GT... On the plus strand that is genomic AG before the body and GT after;
    on the minus strand the reverse complement makes it genomic AC before and CT after the body."""
    codes = {base: ord(base) for base in "ACGT"}
    starts = {"+": [], "-": []}
    names = {"+": [], "-": []}
    end_masks = {"+": {}, "-": {}}
    for name, sequence in sequences.items():
        arr = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        if len(arr) < 3:
            continue
        left, right = arr[:-1], arr[1:]
        dimer = {pair: (left == codes[pair[0]]) & (right == codes[pair[1]]) for pair in ("AG", "GT", "AC", "CT")}
        end_masks["+"][name] = dimer["GT"]
        end_masks["-"][name] = dimer["CT"]
        for strand, acceptor in (("+", "AG"), ("-", "AC")):
            start = np.nonzero(dimer[acceptor])[0] + 3
            start = start[start <= len(sequence)]
            starts[strand].append(start.astype(np.int64))
            names[strand].append(np.full(start.shape, name, dtype=object))
    candidates = {}
    for strand in ("+", "-"):
        position = np.concatenate(starts[strand]) if starts[strand] else np.empty(0, dtype=np.int64)
        contig = np.concatenate(names[strand]) if names[strand] else np.empty(0, dtype=object)
        candidates[strand] = (position, contig)
    return candidates, end_masks


def sample_aggt_controls(contigs, exons, sequences, seed, length_offset=0, attempts=200):
    masks = build_exon_mask(contigs, exons)
    candidates, end_masks = aggt_candidates(sequences)
    if candidates["+"][0].size == 0 and candidates["-"][0].size == 0:
        raise ValueError("No AG/GT-flanked non-exon sites available in the provided FASTA")
    rng = np.random.default_rng(seed)
    controls = []
    for _, start0, end0, _ in exons:
        need = end0 - start0 + 1 + length_offset
        for _ in range(attempts):
            strand = "+" if rng.random() < 0.5 else "-"
            position, contig = candidates[strand]
            if position.size == 0:
                continue
            index = int(rng.integers(position.size))
            name, start = str(contig[index]), int(position[index])
            end = start + need - 1
            if end > contigs[name]["length"] - 2:
                continue
            if not end_masks[strand][name][end]:
                continue
            if masks[name][start - 1:end].any():
                continue
            controls.append((name, start, end, strand))
            break
    if not controls:
        raise ValueError("Could not place any AG/GT-flanked non-exon control at the requested lengths")
    return controls


def render_group(names, matrix_dir, fasta, args, output_dirs):
    matrix, contigs, widths = load_matrix(matrix_dir)
    exons = read_exons(args.gff, contigs)
    requested = {name.removeprefix("abs-").removesuffix("-mutated") for name in names}
    intervals = {"exons": exons}
    if "control" in requested:
        intervals["control"] = sample_controls(contigs, exons, args.seed)
    if "control-aggt" in requested:
        sequences = read_sequences(fasta, contigs)
        intervals["control-aggt"] = sample_aggt_controls(contigs, exons, sequences, args.seed, length_offset=1)

    for name in names:
        output = output_dirs[name]
        group = name.removeprefix("abs-").removesuffix("-mutated")
        if name.startswith("abs-"):
            render_absolute_intervals(intervals[group], matrix, contigs, widths, output, args.rows, name)
        else:
            render_intervals(intervals[group], matrix, contigs, widths, output,
                             args.bins, args.flank, args.rows, name)
            if group == "exons":
                plot_exon_length_distribution(exons, output)


def draw_plot(args):
    selected = [name for name in PLOT_NAMES if not getattr(args, f"no_{name.replace('-', '_')}")]
    if not selected:
        raise ValueError("At least one plot must be selected")
    args.out.mkdir(parents=True, exist_ok=False)
    output_dirs = {name: args.out / name for name in selected}
    for output in output_dirs.values():
        output.mkdir()

    groups = (([name for name in selected if not name.endswith("-mutated")], args.fasta, "original"),
              ([name for name in selected if name.endswith("-mutated")], None, "mutated"))
    with tempfile.TemporaryDirectory(prefix=".cwt-work-", dir=args.out) as temporary:
        work = pathlib.Path(temporary)
        for group_names, fasta, label in groups:
            if not group_names:
                continue
            group_work = work / label
            group_work.mkdir()
            if fasta is None:
                fasta = make_mutated_fasta(args.fasta, args.gff, args.seed, work / "mutated.fa")
            matrix_dir = run_anno_cwt(args.anno_cwt, fasta, group_work)
            print(f"Plotting {', '.join(group_names)}", file=sys.stderr, flush=True)
            render_group(group_names, matrix_dir, fasta, args, output_dirs)


def read_bed(path, contigs):
    intervals = []
    with open_text(path) as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE), 1):
            if not row or row[0].startswith(("#", "track", "browser")):
                continue
            if len(row) < 3:
                raise ValueError(f"Expected at least 3 BED columns at line {line_number}")
            name, start, end = row[0], int(row[1]) + 1, int(row[2])
            strand = row[5] if len(row) > 5 and row[5] in ("+", "-") else "+"
            if name not in contigs:
                raise ValueError(f"BED contig absent from CWT matrix: {name}")
            oriented_bounds(contigs[name], start, end, strand)
            intervals.append((name, start, end, strand))
    if not intervals:
        raise ValueError("No BED intervals found")
    return intervals


def draw_roi(intervals, matrix, contigs, widths, out, flank, pixels):
    for name, start, end, strand in intervals:
        contig = contigs[name]
        first, last = oriented_bounds(contig, start, end, strand)
        before, after = min(flank, first), min(flank, contig["length"] - last)
        values = matrix[contig[strand] + first - before:contig[strand] + last + after]
        if len(values) > 100000:
            raise ValueError("Use a region of at most 100,000 bp for real/imaginary/phase detail")
        bins = min(pixels, len(values))
        body_end = before + last - first
        panels = ((resample_mean(np.real(values), bins).T, "Real (bin mean)", "RdBu_r"),
                  (resample_mean(np.imag(values), bins).T, "Imaginary (bin mean)", "RdBu_r"),
                  (np.log1p(resample_mean(power(values), bins)).T, "log(1 + mean power)", "viridis"))
        figure, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True, layout="constrained")
        for axis, (panel, title, cmap) in zip(axes, panels):
            maximum = max(float(np.nanmax(np.abs(panel))), 1e-6)
            image = axis.imshow(panel, aspect="auto", origin="lower", interpolation="nearest",
                                extent=(0, len(values), -0.5, len(widths) - 0.5), cmap=cmap,
                                vmin=0 if cmap == "viridis" else -maximum, vmax=maximum)
            axis.set_yticks(range(len(widths)), widths)
            axis.set_ylabel(f"{contig.get('parameter', 'Kernel width')} (bp)")
            axis.set_title(title, loc="left", fontsize=11)
            figure.colorbar(image, ax=axis, pad=0.01)
        left_boundary = ("start", start) if strand == "+" else ("end", end)
        right_boundary = ("end", end) if strand == "+" else ("start", start)
        ticks = []
        if before:
            ticks.append((0, f"5' flank\n{before} bp", "left"))
        ticks.append((before, f"{left_boundary[0]}\n{left_boundary[1]:,}", "right"))
        ticks.append((body_end, f"{right_boundary[0]}\n{right_boundary[1]:,}", "left"))
        if after:
            ticks.append((len(values), f"3' flank\n{after} bp", "right"))
        axes[-1].set_xticks([position for position, _, _ in ticks],
                            [label for _, label, _ in ticks])
        for label, (_, _, alignment) in zip(axes[-1].get_xticklabels(), ticks):
            label.set_horizontalalignment(alignment)
        axes[-1].tick_params(axis="x", labelsize=8, pad=6)
        axes[-1].set_xlabel("5' to 3' position; start/end mark genomic exon boundaries")
        figure.suptitle(f"{name}:{start:,}-{end:,} ({strand}) | CWT")
        tag = "plus" if strand == "+" else "minus"
        figure.savefig(out / f"{name}_{start:06d}_{end:06d}_{tag}.png", dpi=160)
        plt.close(figure)


def region_plot(args):
    matrix, contigs, widths = load_matrix(args.matrix)
    draw_roi([(args.seqid, args.start, args.end, args.strand)], matrix, contigs, widths,
             args.out, args.flank, args.pixels)


def regions_plot(args):
    matrix, contigs, widths = load_matrix(args.matrix)
    draw_roi(read_bed(args.bed, contigs), matrix, contigs, widths, args.out, args.flank, args.pixels)


def positive(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("Expected a positive integer")
    return value


def nonnegative(text):
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("Expected a nonnegative integer")
    return value


def add_pipeline_arguments(command, output_name):
    command.add_argument("fasta", type=pathlib.Path)
    command.add_argument("--out", type=pathlib.Path, default=pathlib.Path(output_name))
    command.add_argument("--anno-cwt", type=pathlib.Path, default=pathlib.Path("./anno_cwt"))


def main():
    parser = argparse.ArgumentParser(description="Generate CWT with anno_cwt, then plot and normalize genomic intervals")
    commands = parser.add_subparsers(dest="command", required=True)
    draw = commands.add_parser("draw", help="run exon/control plots; --no-* skips a plot")
    add_pipeline_arguments(draw, "plots-draw")
    draw.add_argument("--gff", type=pathlib.Path, required=True)
    draw.add_argument("--bins", type=positive, default=200)
    draw.add_argument("--flank", type=nonnegative, default=100)
    draw.add_argument("--rows", type=positive, default=1500)
    draw.add_argument("--seed", type=int, default=42)
    for name in PLOT_NAMES:
        draw.add_argument(f"--no-{name}", action="store_true")
    draw.set_defaults(handler=draw_plot)
    for name, handler in (("region", region_plot), ("regions", regions_plot)):
        command = commands.add_parser(name)
        add_pipeline_arguments(command, name)
        if name == "region":
            command.add_argument("--seqid", required=True)
            command.add_argument("--start", type=positive, required=True)
            command.add_argument("--end", type=positive, required=True)
            command.add_argument("--strand", choices=("+", "-"), default="+")
        else:
            command.add_argument("--bed", type=pathlib.Path, required=True)
        command.add_argument("--flank", type=nonnegative, default=10)
        command.add_argument("--pixels", type=positive, default=1600)
        command.set_defaults(handler=handler)
    args = parser.parse_args()
    try:
        if args.command == "draw":
            args.handler(args)
        else:
            args.out.mkdir(parents=True, exist_ok=False)
            with tempfile.TemporaryDirectory(prefix=".cwt-work-", dir=args.out) as temporary:
                args.matrix = run_anno_cwt(args.anno_cwt, args.fasta, pathlib.Path(temporary))
                args.handler(args)
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()