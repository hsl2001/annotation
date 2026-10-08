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
plt.rcParams["svg.fonttype"] = "none"
import numpy as np

PLOT_NAMES = ("abs-exons", "abs-exons-mutated", "abs-control", "abs-control-aggt",
              "abs-control-aggt-p1", "abs-first-exons", "abs-last-exons",
              "abs-control-no-repeat", "abs-control-aggt-no-repeat", "abs-control-aggt-p1-no-repeat")
REPEAT_FEATURE_TYPES = frozenset(("repeat_region", "mobile_element"))
CONTROL_REPLICATES = 2
ABS_RADIUS = 300

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


def gff_records(path):
    with open_text(path) as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE), 1):
            if not row or row[0].startswith("#"):
                if row and row[0] == "##FASTA":
                    break
                continue
            if len(row) != 9:
                raise ValueError(f"Expected 9 GFF/GTF columns at line {line_number}")
            yield row


def gff_attributes(text):
    return {key.strip(): value for key, value in
            (field.split("=", 1) for field in text.split(";") if "=" in field)}


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


def read_exon_sets(path, contigs):
    intervals = set()
    gene_ids = set()
    parents_by_id = {}
    exon_rows = []
    for row in gff_records(path):
        attributes = gff_attributes(row[8])
        identifier = attributes.get("ID")
        if row[2] == "gene" and identifier:
            gene_ids.add(identifier)
        if identifier:
            parents_by_id[identifier] = attributes.get("Parent", "").split(",")
        if row[2] == "exon":
            exon_rows.append((row, attributes.get("Parent", "").split(",")))

    transcript_exons = defaultdict(lambda: defaultdict(set))
    for row, parents in exon_rows:
        name, start, end, strand = row[0], int(row[3]), int(row[4]), row[6]
        if name not in contigs:
            raise ValueError(f"Exon contig absent from CWT matrix: {name}")
        interval = name, start, end, strand
        oriented_bounds(contigs[name], start, end, strand)
        intervals.add(interval)
        for parent in filter(None, parents):
            ancestors, pending = set(), [parent]
            while pending:
                ancestor = pending.pop()
                if ancestor in ancestors:
                    continue
                ancestors.add(ancestor)
                if ancestor in gene_ids:
                    transcript_exons[ancestor][parent].add(interval)
                else:
                    pending.extend(parents_by_id.get(ancestor, ()))
    if not intervals:
        raise ValueError("No exon features found; CDS is not silently substituted for exon")
    first_exons, last_exons = set(), set()
    for transcripts in transcript_exons.values():
        for exons in transcripts.values():
            strand = next(iter(exons))[3]
            if strand == "+":
                first_position = min(exon[1] for exon in exons)
                last_position = max(exon[2] for exon in exons)
                first_exons.update(exon for exon in exons if exon[1] == first_position)
                last_exons.update(exon for exon in exons if exon[2] == last_position)
            else:
                first_position = max(exon[2] for exon in exons)
                last_position = min(exon[1] for exon in exons)
                first_exons.update(exon for exon in exons if exon[2] == first_position)
                last_exons.update(exon for exon in exons if exon[1] == last_position)

    sort_intervals = lambda selected: sorted(
        selected, key=lambda item: (item[2] - item[1] + 1, item))
    return {
        "exons": sort_intervals(intervals),
        "first-exons": sort_intervals(first_exons),
        "last-exons": sort_intervals(last_exons),
    }


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


def boundary_power(matrix, contig, interval, boundary, radius=ABS_RADIUS):
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
    image_values = np.log1p(display)
    maximum = max(float(np.quantile(image_values[np.isfinite(image_values)], 0.995)), 1e-6)
    parameter = next(iter(contigs.values())).get("parameter", "Kernel width")
    positions = np.arange(-ABS_RADIUS, ABS_RADIUS + 1)
    extent = (-ABS_RADIUS - 0.5, ABS_RADIUS + 0.5, 0, len(intervals))
    for boundary_index, boundary in enumerate(("start", "end")):
        figure, axes = plt.subplots(2, len(widths), figsize=(4.1 * len(widths), 6),
                                    gridspec_kw={"height_ratios": [4, 1.3]}, layout="constrained", squeeze=False)
        for scale, width in enumerate(widths):
            axis, profile = axes[:, scale]
            image = axis.imshow(image_values[:, boundary_index, scale], aspect="auto", origin="lower",
                                interpolation="nearest", extent=extent, cmap="viridis", vmin=0, vmax=maximum)
            axis.set_title(f"{parameter} {width} bp")
            axis.set_ylabel("Interval rank" if scale == 0 else "")
            profile.plot(positions, profiles[boundary_index, scale], color="#256c87", linewidth=1.3)
            profile.set_ylabel("Mean power" if scale == 0 else "")
            profile.set_ylim(0.6, 1.2)
            for panel in (axis, profile):
                panel.set_xticks((-ABS_RADIUS, 0, ABS_RADIUS),
                                (f"-{ABS_RADIUS}", "0", f"+{ABS_RADIUS}"), fontsize=8)
                panel.set_xlim(extent[:2])
            profile.set_xlabel(f"Position from {boundary} boundary (bp)")
        figure.colorbar(image, ax=list(axes[0]), label="log(1 + power), common scale; top 0.5% clipped",
                        shrink=0.65, pad=0.01)
        figure.savefig(out / f"{boundary}.png", dpi=160)
        plt.close(figure)


def splice_positions(path, contigs):
    transcripts = defaultdict(list)
    for row in gff_records(path):
        if row[2] != "exon":
            continue
        name, start, end = row[0], int(row[3]), int(row[4])
        if name not in contigs:
            raise ValueError(f"Exon contig absent from assembly FASTA: {name}")
        for parent in filter(None, gff_attributes(row[8]).get("Parent", "").split(",")):
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


def make_mutated_fasta(fasta, gff, seed, output):
    sequences = {}
    for name, sequence in fasta_records(fasta):
        if name in sequences:
            raise ValueError(f"Duplicate FASTA contig: {name}")
        sequences[name] = bytearray(sequence, "ascii")
    if not sequences:
        raise ValueError("No sequences found in assembly FASTA")
    positions = splice_positions(gff, {name: len(sequence) for name, sequence in sequences.items()})
    rng = np.random.default_rng(seed)
    alphabet = "ACGT"
    for name, sites in positions.items():
        sequence = sequences[name]
        for index in sorted(sites):
            current = chr(sequence[index]).upper()
            options = alphabet.replace(current, "") if current in alphabet else alphabet
            sequence[index] = ord(options[int(rng.integers(len(options)))])
    with open(output, "w") as stream:
        for name, sequence in sequences.items():
            stream.write(f">{name}\n")
            for offset in range(0, len(sequence), 60):
                stream.write(sequence[offset:offset + 60].decode("ascii") + "\n")
    return output


def build_exon_mask(contigs, exons):
    masks = {name: np.zeros(info["length"], dtype=bool) for name, info in contigs.items()}
    for name, start, end, _ in exons:
        masks[name][start - 1:end] = True
    return masks


def read_repeat_masks(path, contigs):
    masks = {name: np.zeros(info["length"], dtype=bool) for name, info in contigs.items()}
    found = False
    for row in gff_records(path):
        if row[2] not in REPEAT_FEATURE_TYPES:
            continue
        name, start, end = row[0], int(row[3]), int(row[4])
        if name not in contigs:
            continue
        if start < 1 or end < start or end > contigs[name]["length"]:
            raise ValueError(f"Invalid repeat interval: {name}:{start}-{end}")
        masks[name][start - 1:end] = True
        found = True
    if not found:
        return None
    return masks


def control_exclusion_masks(contigs, exons, repeat_masks=None):
    masks = build_exon_mask(contigs, exons)
    if repeat_masks is not None:
        for name in masks:
            masks[name] |= repeat_masks[name]
    return masks


def sample_controls(contigs, exons, seed, attempts=200, repeat_masks=None):
    masks = control_exclusion_masks(contigs, exons, repeat_masks)
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


def prepare_aggt_controls(contigs, exons, sequences, repeat_masks=None):
    masks = control_exclusion_masks(contigs, exons, repeat_masks)
    candidates, end_masks = aggt_candidates(sequences)
    if candidates["+"][0].size == 0 and candidates["-"][0].size == 0:
        raise ValueError("No AG/GT-flanked non-exon sites available in the provided FASTA")
    return masks, candidates, end_masks


def sample_aggt_controls(contigs, exons, sequences, seed, length_offset=0, attempts=200, prepared=None):
    masks, candidates, end_masks = prepared or prepare_aggt_controls(contigs, exons, sequences)
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
    intervals = read_exon_sets(args.gff, contigs)
    repeat_names = [name for name in names if name.endswith("-no-repeat")]
    repeat_masks = read_repeat_masks(args.gff, contigs) if repeat_names else None
    if repeat_names and repeat_masks is None:
        print("Skipping repeat-free plots: no repeat annotations on selected contigs",
              file=sys.stderr, flush=True)
        for name in repeat_names:
            output_dirs[name].rmdir()
        names = [name for name in names if name not in repeat_names]
    prepared = {}
    if any("control-aggt" in name for name in names):
        sequences = read_sequences(fasta, contigs)
        prepared[False] = prepare_aggt_controls(contigs, intervals["exons"], sequences)
        if repeat_masks is not None:
            masks, candidates, end_masks = prepared[False]
            prepared[True] = ({name: mask | repeat_masks[name] for name, mask in masks.items()},
                              candidates, end_masks)
    for name in names:
        group = name.removeprefix("abs-").removesuffix("-mutated")
        base_group = group.removesuffix("-no-repeat")
        repeat_free = group.endswith("-no-repeat")
        aggt = name.startswith("abs-control-aggt")
        for replicate in range(1, CONTROL_REPLICATES + 1) if aggt else (0,):
            output = output_dirs[name]
            if aggt:
                offset = int(base_group.endswith("-p1"))
                selected = sample_aggt_controls(
                    contigs, intervals["exons"], sequences,
                    args.seed + offset * CONTROL_REPLICATES + replicate - 1,
                    length_offset=offset, prepared=prepared[repeat_free])
                output = output / f"replicate-{replicate}"
                output.mkdir()
            elif base_group == "control":
                selected = sample_controls(contigs, intervals["exons"], args.seed,
                                           repeat_masks=repeat_masks if repeat_free else None)
            else:
                selected = intervals[group]
            if repeat_free:
                print(f"Sampled {len(selected):,}/{len(intervals['exons']):,} {group} "
                      f"{output.name}; unmatched lengths skipped after 200 attempts", file=sys.stderr, flush=True)
            render_absolute_intervals(selected, matrix, contigs, widths, output, args.rows, name)


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
    intervals = ([(args.seqid, args.start, args.end, args.strand)] if args.command == "region"
                 else read_bed(args.bed, contigs))
    draw_roi(intervals, matrix, contigs, widths, args.out, args.flank, args.pixels)


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
    draw = commands.add_parser(
        "draw", help="run absolute exon/control plots; --no-* skips a plot",
        description="draw generates abs-* patterns only. Repeat-free controls exclude bodies overlapping GFF "
                    "repeat_region or mobile_element features; plotted flanks are not excluded. "
                    f"AG/GT variants use {CONTROL_REPLICATES} replicates; abs plots use "
                    f"+/-{ABS_RADIUS} bp boundary windows.")
    add_pipeline_arguments(draw, "plots-draw")
    draw.add_argument("--gff", type=pathlib.Path, required=True)
    draw.add_argument("--rows", type=positive, default=1500)
    draw.add_argument("--seed", type=int, default=42)
    for name in PLOT_NAMES:
        draw.add_argument(f"--no-{name}", action="store_true")
    draw.set_defaults(handler=draw_plot)
    for name in ("region", "regions"):
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
        command.set_defaults(handler=region_plot)
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