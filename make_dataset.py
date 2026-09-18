#!/usr/bin/env python3
"""Combine FASTA/GFF3 dataset pairs in a directory without identifier collisions."""

import argparse
import gzip
import pathlib
import re
import sys
from collections import defaultdict

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
GFF_SUFFIXES = (".gff", ".gff3")
IDENTIFIER_ATTRIBUTES = {
    "ID", "Parent", "Derives_from", "gene_id", "transcript_id",
    "exon_id", "protein_id",
}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def base_suffix(path, suffixes):
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    for suffix in suffixes:
        if name.endswith(suffix):
            return suffix
    return None


def source_files(directory, suffixes):
    return sorted(
        path for path in directory.iterdir()
        if path.is_file() and base_suffix(path, suffixes) is not None
    )


def fasta_sequences(path):
    sequences = {}
    current = None
    length = 0
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = length
                current = line[1:].split(None, 1)[0]
                if not current:
                    raise ValueError(f"{path}:{line_number}: empty FASTA sequence ID")
                if current in sequences:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA sequence ID {current}")
                length = 0
            elif current is None:
                if line.strip():
                    raise ValueError(f"{path}:{line_number}: FASTA sequence before header")
            else:
                length += len(line.strip())
    if current is not None:
        sequences[current] = length
    if not sequences:
        raise ValueError(f"{path}: no FASTA sequences found")
    return sequences


def gff_sequences(path):
    sequence_ids = set()
    maximum = defaultdict(int)
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 GFF3 fields")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid GFF3 coordinates") from error
            if start < 1 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid GFF3 interval")
            sequence_ids.add(fields[0])
            maximum[fields[0]] = max(maximum[fields[0]], end)
    if not sequence_ids:
        raise ValueError(f"{path}: no GFF3 records found")
    return sequence_ids, maximum


def pair_score(fasta, gff):
    fasta_lengths = fasta_sequences(fasta)
    gff_ids, gff_maximum = gff_sequences(gff)
    shared = gff_ids & fasta_lengths.keys()
    valid = {seqid for seqid in shared if gff_maximum[seqid] <= fasta_lengths[seqid]}
    return len(valid), len(gff_ids), fasta_lengths, gff_ids


def pair_inputs(fastas, gffs):
    candidates = []
    for fasta in fastas:
        for gff in gffs:
            valid, total, _, _ = pair_score(fasta, gff)
            if valid == total:
                candidates.append((valid, fasta, gff))
    if not candidates:
        raise ValueError("could not pair any FASTA with a GFF3 by sequence IDs and coordinates")

    pairs = []
    remaining_fasta = set(fastas)
    remaining_gff = set(gffs)
    while remaining_fasta and remaining_gff:
        available = [candidate for candidate in candidates
                     if candidate[1] in remaining_fasta and candidate[2] in remaining_gff]
        if not available:
            break
        best_score = max(candidate[0] for candidate in available)
        best = [candidate for candidate in available if candidate[0] == best_score]
        if len(best) != 1:
            details = ", ".join(f"{f.name} + {g.name}" for _, f, g in best)
            raise ValueError(f"ambiguous FASTA/GFF3 pairing: {details}; use --fasta and --gff")
        _, fasta, gff = best[0]
        pairs.append((fasta, gff))
        remaining_fasta.remove(fasta)
        remaining_gff.remove(gff)
    if remaining_fasta or remaining_gff:
        missing = ", ".join(path.name for path in sorted(remaining_fasta | remaining_gff))
        raise ValueError(f"unpaired input files: {missing}")
    return pairs


def prefix_for(path):
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    for suffix in FASTA_SUFFIXES + GFF_SUFFIXES:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_.-") or "dataset"
    return name + "__"


def prefixed_identifier(value, prefix):
    return ",".join(prefix + item for item in value.split(","))


def rewrite_attributes(text, prefix):
    fields = []
    for field in text.split(";"):
        if "=" not in field:
            fields.append(field)
            continue
        key, value = field.split("=", 1)
        if key in IDENTIFIER_ATTRIBUTES:
            value = prefixed_identifier(value, prefix)
        fields.append(f"{key}={value}")
    return ";".join(fields)


def write_fasta(path, output, prefix):
    with open_text(path) as source:
        for line_number, line in enumerate(source, 1):
            if line.startswith(">"):
                fields = line[1:].rstrip("\r\n").split(None, 1)
                if not fields or not fields[0]:
                    raise ValueError(f"{path}:{line_number}: empty FASTA sequence ID")
                header = ">" + prefix + fields[0]
                if len(fields) == 2:
                    header += " " + fields[1]
                output.write(header + "\n")
            else:
                output.write(line)


def write_gff(path, output, prefix):
    with open_text(path) as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 GFF3 fields")
            fields[0] = prefix + fields[0]
            fields[8] = rewrite_attributes(fields[8], prefix)
            output.write("\t".join(fields) + "\n")


def parse_explicit_pairs(values):
    pairs = []
    for value in values:
        fasta, separator, gff = value.partition(":")
        if not separator or not fasta or not gff:
            raise ValueError(f"invalid --pair value {value!r}; use FASTA:GFF3")
        pairs.append((pathlib.Path(fasta), pathlib.Path(gff)))
    return pairs


def main():
    parser = argparse.ArgumentParser(description="Combine FASTA/GFF3 datasets with prefixed IDs")
    parser.add_argument("-d", "--directory", type=pathlib.Path, required=True,
                        help="directory containing FASTA and GFF/GFF3 files")
    parser.add_argument("--pair", action="append", default=[], metavar="FASTA:GFF3",
                        help="explicit input pair; repeat for multiple datasets")
    parser.add_argument("--fasta", nargs="+", type=pathlib.Path,
                        help="explicit FASTA files, paired by order with --gff")
    parser.add_argument("--gff", nargs="+", type=pathlib.Path,
                        help="explicit GFF/GFF3 files, paired by order with --fasta")
    parser.add_argument("-o", "--output-prefix", default="combined",
                        help="output basename inside -d (default: combined)")
    args = parser.parse_args()

    try:
        if args.pair:
            if args.fasta or args.gff:
                parser.error("use --pair or --fasta/--gff, not both")
            pairs = parse_explicit_pairs(args.pair)
        elif args.fasta or args.gff:
            if not args.fasta or not args.gff or len(args.fasta) != len(args.gff):
                parser.error("--fasta and --gff must be supplied with the same number of files")
            pairs = list(zip(args.fasta, args.gff))
        else:
            fastas = source_files(args.directory, FASTA_SUFFIXES)
            gffs = source_files(args.directory, GFF_SUFFIXES)
            outputs = {args.directory / f"{args.output_prefix}.fa",
                       args.directory / f"{args.output_prefix}.gff3"}
            fastas = [path for path in fastas if path not in outputs]
            gffs = [path for path in gffs if path not in outputs]
            pairs = pair_inputs(fastas, gffs)
        if not pairs:
            raise ValueError("no input pairs found")
        for fasta, gff in pairs:
            if not fasta.is_file() or not gff.is_file():
                raise ValueError(f"input file not found: {fasta} or {gff}")

        args.directory.mkdir(parents=True, exist_ok=True)
        output_fasta = args.directory / f"{args.output_prefix}.fa"
        output_gff = args.directory / f"{args.output_prefix}.gff3"
        prefixes = []
        with open(output_fasta, "w") as fasta_out, open(output_gff, "w") as gff_out:
            gff_out.write("##gff-version 3\n")
            used = set()
            for fasta, gff in pairs:
                prefix = prefix_for(fasta)
                if prefix in used:
                    raise ValueError(f"duplicate dataset prefix {prefix}; use explicit unique input names")
                used.add(prefix)
                write_fasta(fasta, fasta_out, prefix)
                write_gff(gff, gff_out, prefix)
                prefixes.append((fasta.name, gff.name, prefix))
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(f"FASTA: {output_fasta}")
    print(f"GFF3:  {output_gff}")
    for fasta, gff, prefix in prefixes:
        print(f"{prefix}: {fasta} + {gff}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
