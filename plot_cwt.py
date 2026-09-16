import argparse
import csv
import gzip
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_matrix(directory):
    with open(directory / "contigs.tsv") as stream:
        if stream.readline().strip() != "# anno-cwt-v1":
            raise ValueError("Unsupported CWT index version")
        fields = stream.readline().rstrip().split("\t")
        if fields[0] != "# kernel_widths":
            raise ValueError("Missing kernel widths")
        widths = [int(value) for value in ",".join(fields[1:]).split(",") if value]
        contigs = {}
        expected = 0
        for row in csv.DictReader(stream, delimiter="\t"):
            name = row["seqid"]
            length, plus, minus = (int(row[key]) for key in ("length", "plus_offset", "minus_offset"))
            if name in contigs or length <= 0 or plus != expected or minus != plus + length:
                raise ValueError(f"Invalid contig offsets: {name}")
            contigs[name] = {"length": length, "+": plus, "-": minus}
            expected += 2 * length
    if not widths or not contigs:
        raise ValueError("Missing CWT metadata")
    path = directory / "matrix.bin"
    expected_bytes = expected * len(widths) * np.dtype("<c8").itemsize
    if not path.is_file() or path.stat().st_size != expected_bytes:
        raise ValueError("Expected a little-endian complex64 matrix matching contigs.tsv")
    matrix = np.memmap(path, dtype="<c8", mode="r", shape=(expected, len(widths)))
    return matrix, contigs, widths


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
    raw_count = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE), 1):
            if not row or row[0].startswith("#"):
                if row and row[0] == "##FASTA":
                    break
                continue
            if len(row) != 9:
                raise ValueError(f"Expected 9 GFF/GTF columns at line {line_number}")
            if row[2] != "exon":
                continue
            name, start, end, strand = row[0], int(row[3]), int(row[4]), row[6]
            if name not in contigs:
                raise ValueError(f"Exon contig absent from CWT matrix: {name}")
            oriented_bounds(contigs[name], start, end, strand)
            raw_count += 1
            intervals.add((name, start, end, strand))
    if not intervals:
        raise ValueError("No exon features found; CDS is not silently substituted for exon")
    return sorted(intervals, key=lambda item: (item[2] - item[1] + 1, item)), raw_count


def finite_mean(values, axis=0):
    valid = np.isfinite(values)
    count = valid.sum(axis=axis)
    total = np.where(valid, values, 0.0).sum(axis=axis, dtype=np.float64)
    return np.divide(total, count, out=np.full_like(total, np.nan), where=count > 0)


_KIND = {
    "exons": ("All {n:,} unique annotated exons | {bins} body bins | {rows:,} display rows\n"
              "Every exon included; display rows average adjacent length-sorted exons",
              "Exon rank (short to long)"),
    "control": ("{n:,} length-matched random non-exon controls | {bins} body bins | {rows:,} display rows\n"
                "Arbitrary genomic cut points; body avoids every annotated exon", "Control rank (short to long)"),
}


def render_intervals(intervals, matrix, contigs, widths, out, bins, flank, rows, stem):
    intervals = sorted(intervals, key=lambda item: (item[2] - item[1] + 1, item))
    out.mkdir(parents=True, exist_ok=False)
    columns = bins + 2 * flank
    normalized = np.lib.format.open_memmap(out / f"{stem}.power.npy", mode="w+", dtype="float32",
                                           shape=(len(intervals), len(widths), columns))
    total = np.zeros((len(widths), columns), dtype=np.float64)
    counts = np.zeros_like(total, dtype=np.int64)
    with open(out / f"{stem}.tsv", "w") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("row", "seqid", "start", "end", "strand", "length"))
        for row, (name, start, end, strand) in enumerate(intervals):
            values = exon_power(matrix, contigs[name], start, end, strand, bins, flank)
            normalized[row] = values
            valid = np.isfinite(values)
            total += np.where(valid, values, 0.0)
            counts += valid
            writer.writerow((row, name, start, end, strand, end - start + 1))
            if (row + 1) % 10000 == 0:
                print(f"Normalized {row + 1:,}/{len(intervals):,} {stem}", file=sys.stderr, flush=True)
    normalized.flush()
    del normalized
    normalized = np.load(out / f"{stem}.power.npy", mmap_mode="r", allow_pickle=False)
    display_rows = min(rows, len(intervals))
    edges = np.linspace(0, len(intervals), display_rows + 1).astype(np.int64)
    display = np.empty((display_rows, len(widths), columns), dtype=np.float32)
    for row, (first, last) in enumerate(zip(edges[:-1], edges[1:])):
        display[row] = finite_mean(normalized[first:last])
    mean = np.divide(total, counts, out=np.full_like(total, np.nan), where=counts > 0)
    np.save(out / "mean_power.npy", mean)
    image_values = np.log1p(display)
    finite = image_values[np.isfinite(image_values)]
    maximum = max(float(np.quantile(finite, 0.995)), 1e-6)
    figure, axes = plt.subplots(2, len(widths), figsize=(4.1 * len(widths), 9),
                                gridspec_kw={"height_ratios": [4, 1.3]}, layout="constrained", squeeze=False)
    ticks = [flank, flank + bins / 2, flank + bins]
    labels = ["0%", "50%", "100%"]
    if flank:
        ticks = [0, *ticks, columns]
        labels = [f"-{flank} bp", *labels, f"+{flank} bp"]
    title, rank_label = _KIND[stem]
    for scale, width in enumerate(widths):
        axis, profile = axes[:, scale]
        image = axis.imshow(image_values[:, scale], aspect="auto", origin="lower", interpolation="nearest",
                            extent=(0, columns, 0, len(intervals)), cmap="viridis", vmin=0, vmax=maximum)
        axis.set_title(f"Kernel width {width} bp")
        axis.set_ylabel(rank_label if scale == 0 else "")
        profile.plot(np.arange(columns) + 0.5, mean[scale], color="#256c87", linewidth=1.3)
        profile.set_ylabel("Mean power" if scale == 0 else "")
        for panel in (axis, profile):
            panel.axvline(flank, color="#bf5547", linestyle="--", linewidth=0.8)
            panel.axvline(flank + bins, color="#bf5547", linestyle="--", linewidth=0.8)
            panel.set_xticks(ticks, labels, fontsize=8)
            panel.set_xlim(0, columns)
        profile.set_xlabel("5' to 3': normalized body + genomic flanks")
    figure.colorbar(image, ax=list(axes[0]), label="log(1 + mean power), common scale; top 0.5% clipped",
                    shrink=0.65, pad=0.01)
    figure.suptitle(title.format(n=len(intervals), bins=bins, rows=display_rows), fontsize=14)
    figure.savefig(out / f"{stem}.png", dpi=160)
    plt.close(figure)
    meta = {
        "count": len(intervals), "display_rows": display_rows,
        "min_len": intervals[0][2] - intervals[0][1] + 1, "max_len": intervals[-1][2] - intervals[-1][1] + 1,
        "plus": sum(item[3] == "+" for item in intervals), "minus": sum(item[3] == "-" for item in intervals),
        "contigs": sorted({item[0] for item in intervals}),
        "mean_body_power_by_scale": finite_mean(mean[:, flank:flank + bins], axis=1).tolist(),
    }
    return mean, meta


def exon_plot(args):
    matrix, contigs, widths = load_matrix(args.matrix)
    exons, raw_count = read_exons(args.gff, contigs)
    mean, meta = render_intervals(exons, matrix, contigs, widths, args.out,
                                  args.bins, args.flank, args.rows, "exons")
    summary = {
        "matrix": str(args.matrix.resolve()), "annotation": str(args.gff.resolve()),
        "raw_exon_rows": raw_count, "unique_exons": meta["count"],
        "duplicate_exon_rows": raw_count - meta["count"],
        "plus_exons": meta["plus"], "minus_exons": meta["minus"],
        "contigs": meta["contigs"], "kernel_widths": widths,
        "bins": args.bins, "flank_bp": args.flank, "display_rows": meta["display_rows"],
        "min_exon_bp": meta["min_len"], "max_exon_bp": meta["max_len"],
        "mean_exon_power_by_scale": meta["mean_body_power_by_scale"],
        "orientation": "Transcript 5-prime to 3-prime; reverse-complement CWT for minus strand",
        "normalization": "Area-weighted mean power per bin; length only, no per-exon amplitude normalization",
        "edge_policy": "Contig-external flanks are NaN and excluded from means",
        "reference_policy": "Every exon feature, including UTR and noncoding exons; identical strand/coordinates deduplicated",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


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
    controls, skipped = [], 0
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
        else:
            skipped += 1
    if not controls:
        raise ValueError("Could not place any non-exon control at the requested lengths")
    return controls, skipped


def control_plot(args):
    matrix, contigs, widths = load_matrix(args.matrix)
    exons, _ = read_exons(args.gff, contigs)
    controls, skipped = sample_controls(contigs, exons, args.seed)
    mean, meta = render_intervals(controls, matrix, contigs, widths, args.out,
                                  args.bins, args.flank, args.rows, "control")
    summary = {
        "matrix": str(args.matrix.resolve()), "annotation": str(args.gff.resolve()), "seed": args.seed,
        "controls": meta["count"], "skipped_lengths": skipped,
        "plus": meta["plus"], "minus": meta["minus"], "contigs": meta["contigs"], "kernel_widths": widths,
        "bins": args.bins, "flank_bp": args.flank, "display_rows": meta["display_rows"],
        "min_len": meta["min_len"], "max_len": meta["max_len"],
        "mean_body_power_by_scale": meta["mean_body_power_by_scale"],
        "sampling": "Length-matched to annotated exons; body fully non-exonic; random contig, position and strand",
        "orientation": "5-prime to 3-prime in the sampled strand; reverse-complement CWT for minus strand",
        "normalization": "Area-weighted mean power per bin; identical to exon normalization",
        "edge_policy": "Contig-external flanks are NaN and excluded from means",
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.compare:
        exon_mean = np.load(args.compare / "mean_power.npy")
        if exon_mean.shape != mean.shape:
            raise ValueError("Comparison mean_power.npy shape differs; rerun with matching bins/flank")
        low5, high5 = max(0, args.flank - 5), args.flank + 6
        low3, high3 = args.flank + args.bins - 6, min(mean.shape[1], args.flank + args.bins + 5)
        print("\nLocal boundary peak power, exon vs non-exon control (window +/-5 bp):")
        print(f"{'width':>6} {'boundary':>10} {'exon':>9} {'control':>9} {'exon/ctrl':>10}")
        for scale, width in enumerate(widths):
            for label, (low, high) in (("5'", (low5, high5)), ("3'", (low3, high3))):
                e = float(np.nanmax(exon_mean[scale, low:high]))
                c = float(np.nanmax(mean[scale, low:high]))
                print(f"{width:>6} {label:>10} {e:>9.4f} {c:>9.4f} {e / c if c else float('nan'):>10.2f}")


def region_plot(args):
    matrix, contigs, widths = load_matrix(args.matrix)
    if args.seqid not in contigs:
        raise ValueError(f"Unknown contig: {args.seqid}")
    contig = contigs[args.seqid]
    first, last = oriented_bounds(contig, args.start, args.end, args.strand)
    values = matrix[contig[args.strand] + first:contig[args.strand] + last]
    if len(values) > 100000:
        raise ValueError("Use a region of at most 100,000 bp for real/imaginary/phase detail")
    bins = min(args.pixels, len(values))
    panels = ((resample_mean(np.real(values), bins).T, "Real (bin mean)", "RdBu_r"),
              (resample_mean(np.imag(values), bins).T, "Imaginary (bin mean)", "RdBu_r"),
              (np.log1p(resample_mean(power(values), bins)).T, "log(1 + mean power)", "viridis"))
    figure, axes = plt.subplots(3, 1, figsize=(13, 8), sharex=True, layout="constrained")
    for axis, (image_values, title, colors) in zip(axes, panels):
        maximum = max(float(np.abs(image_values).max()), 1e-6)
        image = axis.imshow(image_values, aspect="auto", origin="lower", interpolation="nearest",
                            extent=(0, len(values), -0.5, len(widths) - 0.5), cmap=colors,
                            vmin=0 if colors == "viridis" else -maximum, vmax=maximum)
        axis.set_yticks(range(len(widths)), widths)
        axis.set_ylabel("Kernel width (bp)")
        axis.set_title(title, loc="left", fontsize=11)
        figure.colorbar(image, ax=axis, pad=0.01)
    axes[-1].set_xlabel("Offset from 5' end (bp); minus strand runs toward smaller genomic coordinates")
    figure.suptitle(f"{args.seqid}:{args.start:,}-{args.end:,} ({args.strand}) | CWT")
    figure.savefig(args.out, dpi=160)
    plt.close(figure)
    print(args.out)


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


def main():
    parser = argparse.ArgumentParser(description="Memory-mapped CWT visualization and all-exon length normalization")
    commands = parser.add_subparsers(dest="command", required=True)
    region = commands.add_parser("region")
    region.add_argument("matrix", type=pathlib.Path)
    region.add_argument("--seqid", required=True)
    region.add_argument("--start", type=positive, required=True)
    region.add_argument("--end", type=positive, required=True)
    region.add_argument("--strand", choices=("+", "-"), default="+")
    region.add_argument("--pixels", type=positive, default=1600)
    region.add_argument("--out", type=pathlib.Path, required=True)
    region.set_defaults(handler=region_plot)
    exons = commands.add_parser("exons")
    exons.add_argument("matrix", type=pathlib.Path)
    exons.add_argument("--gff", type=pathlib.Path, required=True)
    exons.add_argument("--out", type=pathlib.Path, required=True)
    exons.add_argument("--bins", type=positive, default=200)
    exons.add_argument("--flank", type=nonnegative, default=100)
    exons.add_argument("--rows", type=positive, default=1500)
    exons.set_defaults(handler=exon_plot)
    control = commands.add_parser("control")
    control.add_argument("matrix", type=pathlib.Path)
    control.add_argument("--gff", type=pathlib.Path, required=True)
    control.add_argument("--out", type=pathlib.Path, required=True)
    control.add_argument("--bins", type=positive, default=200)
    control.add_argument("--flank", type=nonnegative, default=100)
    control.add_argument("--rows", type=positive, default=1500)
    control.add_argument("--seed", type=int, default=42)
    control.add_argument("--compare", type=pathlib.Path)
    control.set_defaults(handler=control_plot)
    args = parser.parse_args()
    try:
        args.handler(args)
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    main()