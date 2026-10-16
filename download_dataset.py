#!/usr/bin/env python3
"""Download nested plant training sets using the NCBI Datasets v2 CLI.

Usage:
    python3 download_dataset.py --output-dir training-data
    python3 download_dataset.py --output-dir training-data --select-only
"""

import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import TextIO
from zipfile import BadZipFile, ZipFile

from make_dataset import fasta_sequences, open_text, rewrite_attributes, write_fasta


SIZES = (10, 20, 100)
LAND_PLANTS = "3193"
DATASETS_URL = "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets"
PREFERRED_SPECIES = (
    "Oryza sativa", "Zea mays", "Solanum lycopersicum", "Glycine max",
    "Brachypodium distachyon", "Medicago truncatula", "Lotus japonicus",
    "Physcomitrium patens", "Populus trichocarpa", "Setaria viridis",
    "Sorghum bicolor", "Brassica rapa", "Vitis vinifera", "Cucumis sativus",
    "Fragaria vesca", "Nicotiana tabacum", "Hordeum vulgare",
    "Triticum aestivum", "Helianthus annuus", "Phaseolus vulgaris",
)


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def run(command: list[str], capture: bool = False) -> str:
    log("+ " + " ".join(command))
    result = subprocess.run(command, check=True, text=True,
                            stdout=subprocess.PIPE if capture else None)
    return result.stdout if capture else ""


def species_name(record: dict) -> str:
    # NCBI organism names may include subspecies, cultivars, or strain suffixes.
    return " ".join(record["organism"]["organism_name"].split()[:2])


def eligible(record: dict, min_busco: float) -> bool:
    organism = record.get("organism", {})
    name = organism.get("organism_name", "")
    assembly = record.get("assembly_info", {})
    annotation = record.get("annotation_info", {})
    return (
        organism.get("tax_id") != 3702
        and not name.casefold().startswith("arabidopsis thaliana")
        and len(name.split()) >= 2
        and " x " not in name
        and re.fullmatch(r"GCF_\d+\.\d+", record.get("accession", "")) is not None
        and record.get("current_accession", record["accession"]) == record["accession"]
        and assembly.get("assembly_status") == "current"
        and assembly.get("assembly_level") in {"Chromosome", "Complete Genome"}
        and annotation.get("status") in {"Full annotation", "Updated annotation"}
        and annotation.get("stats", {}).get("gene_counts", {}).get("protein_coding", 0) > 0
        and annotation.get("busco", {}).get("complete", -1) >= min_busco
    )


def quality_key(record: dict) -> tuple:
    assembly = record["assembly_info"]
    return (
        assembly.get("refseq_category") == "reference genome",
        record["annotation_info"]["busco"]["complete"],
        assembly["assembly_level"] == "Complete Genome",
        record["annotation_info"].get("release_date", ""),
        assembly.get("release_date", ""),
        record["accession"],
    )


def select_assemblies(records: list[dict], min_busco: float) -> list[dict]:
    best: dict[str, dict] = {}
    for record in records:
        if eligible(record, min_busco):
            species = species_name(record)
            if species not in best or quality_key(record) > quality_key(best[species]):
                best[species] = record
    priority = {name: index for index, name in enumerate(PREFERRED_SPECIES)}
    ranked = sorted(best.values(), key=quality_key, reverse=True)
    ranked.sort(key=lambda record: priority.get(species_name(record), len(priority)))
    if len(ranked) < max(SIZES):
        raise ValueError(
            f"Only {len(ranked)} eligible species; need {max(SIZES)}. "
            "Inspect NCBI metadata or explicitly lower --min-busco; "
            "quality filters will not be relaxed automatically."
        )
    return ranked[:max(SIZES)]


def save_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def ensure_cli(datasets: str) -> str:
    cli = shutil.which(datasets)
    if cli is not None:
        return cli
    path = Path(datasets) if "/" in datasets else Path(__file__).resolve().parent / "bin" / "datasets"
    if not path.is_file():
        run(["curl", "-fL", "-o", str(path), "--create-dirs", DATASETS_URL])
        run(["chmod", "+x", str(path)])
    return str(path.resolve())


def get_selection(cli: str, cache: Path, min_busco: float) -> list[dict]:
    selection = cache / "selection.json"
    if selection.exists():
        records = json.loads(selection.read_text(encoding="utf-8"))
        if not isinstance(records, list) or len(records) != max(SIZES):
            raise ValueError(f"{selection}: expected exactly {max(SIZES)} assemblies")
        if not all(eligible(record, min_busco) for record in records):
            raise ValueError(f"{selection}: an assembly fails the requested quality/exclusion filters")
        if len({species_name(record) for record in records}) != len(records):
            raise ValueError(f"{selection}: duplicate species")
        log(f"Reusing frozen selection: {selection}")
        return records
    text = run([
        cli, "summary", "genome", "taxon", LAND_PLANTS,
        "--annotated", "--assembly-source", "RefSeq",
        "--assembly-level", "chromosome,complete", "--exclude-atypical",
        "--as-json-lines",
    ], capture=True)
    records = select_assemblies(
        [json.loads(line) for line in text.splitlines() if line.strip()], min_busco
    )
    save_json(selection, records)
    return records


def package_inputs(package: Path, records: list[dict]) -> list[tuple[Path, Path]]:
    data = package / "ncbi_dataset" / "data"
    catalog = json.loads((data / "dataset_catalog.json").read_text(encoding="utf-8"))
    entries = {entry["accession"]: entry for entry in catalog["assemblies"]
               if "accession" in entry}
    pairs = []
    for record in records:
        accession = record["accession"]
        if accession not in entries:
            raise ValueError(f"{accession}: missing from the cached package catalog")
        files = entries[accession]["files"]
        paths = []
        for kind in ("GENOMIC_NUCLEOTIDE_FASTA", "GFF3"):
            matches = [item for item in files if item["fileType"] == kind]
            if len(matches) != 1:
                raise ValueError(f"{accession}: expected one {kind} file in catalog")
            item = matches[0]
            path = (data / item["filePath"]).resolve()
            if not path.is_relative_to(data.resolve()):
                raise ValueError(f"Unsafe catalog path: {item['filePath']}")
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing/empty hydrated file: {path}")
            expected = item.get("uncompressedLengthBytes")
            if expected is not None and path.stat().st_size != int(expected):
                raise ValueError(f"Incomplete hydrated file: {path}")
            paths.append(path)
        pairs.append((paths[0], paths[1]))
    return pairs


def download(cli: str, cache: Path, records: list[dict],
             workers: int) -> list[tuple[Path, Path]]:
    accessions = cache / "accessions.txt"
    accessions.write_text("".join(record["accession"] + "\n" for record in records),
                          encoding="utf-8")
    archive = cache / "dataset.zip"
    package = cache / "package"
    if not package.exists():
        if not archive.exists():
            with tempfile.TemporaryDirectory(prefix="download-", dir=cache) as temporary:
                partial = Path(temporary) / "dataset.zip"
                run([cli, "download", "genome", "accession", "--inputfile", str(accessions),
                     "--include", "genome,gff3", "--dehydrated",
                     "--filename", str(partial), "--no-progressbar"])
                partial.replace(archive)
        with tempfile.TemporaryDirectory(prefix="extract-", dir=cache) as temporary:
            unpacked = Path(temporary) / "package"
            unpacked.mkdir()
            with ZipFile(archive) as zipped:
                for item in zipped.infolist():
                    if not (unpacked / item.filename).resolve().is_relative_to(unpacked.resolve()):
                        raise ValueError(f"Unsafe archive path: {item.filename}")
                zipped.extractall(unpacked)
            unpacked.replace(package)
    # Rehydration resumes interrupted downloads; only mark complete after validation.
    complete = package / ".hydrated"
    if not complete.exists():
        run([cli, "rehydrate", "--directory", str(package),
             "--max-workers", str(workers), "--no-progressbar"])
    pairs = package_inputs(package, records)
    complete.touch()
    return pairs


def gff_attributes(value: str) -> dict[str, str]:
    attributes = {}
    for item in value.split(";"):
        key, separator, attribute_value = item.partition("=")
        if separator:
            attributes[key] = attribute_value
    return attributes


def circular_trans_splice_bounds(path: Path, lengths: dict[str, int]) -> tuple[set[str], dict[str, tuple[int, int]]]:
    circular: set[str] = set()
    wrapped: dict[str, str] = {}
    children: dict[str, list[tuple[str, int, int]]] = {}
    with open_text(path) as source:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[0] not in lengths:
                continue
            attributes = gff_attributes(fields[8])
            if fields[2] == "region" and attributes.get("Is_circular", "").casefold() == "true":
                circular.add(fields[0])
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            if 1 <= start <= lengths[fields[0]] < end and attributes.get("ID"):
                wrapped[attributes["ID"]] = fields[0]
            if attributes.get("exception") == "trans-splicing":
                for parent in attributes.get("Parent", "").split(","):
                    if parent:
                        children.setdefault(parent, []).append((fields[0], start, end))

    adjusted = {}
    for identifier, seqid in wrapped.items():
        if seqid not in circular:
            continue
        length = lengths[seqid]
        segments = []
        for child_seqid, start, end in children.get(identifier, []):
            if child_seqid != seqid:
                continue
            if not 1 <= start <= end <= 2 * length:
                raise ValueError(f"{path}: invalid trans-spliced child coordinates for {identifier}")
            if start <= length < end:
                segments.extend(((start, length), (1, end - length)))
            elif start > length:
                segments.append((start - length, end - length))
            else:
                segments.append((start, end))
        if segments:
            adjusted[identifier] = (
                min(start for start, _ in segments),
                max(end for _, end in segments),
            )
    return circular, adjusted


def write_annotation(path: Path, output: TextIO, prefix: str,
                     lengths: dict[str, int]) -> None:
    identifiers: set[str] = set()
    references: set[str] = set()
    features = 0
    circular, wrapped_bounds = circular_trans_splice_bounds(path, lengths)
    with open_text(path) as source:
        for number, line in enumerate(source, 1):
            location = f"{path}:{number}"
            if line.startswith("##FASTA"):
                raise ValueError(f"{location}: embedded FASTA is not supported")
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) != 4 or fields[1] not in lengths:
                    raise ValueError(f"{location}: unknown/invalid sequence-region")
                start, end = int(fields[2]), int(fields[3])
                if not 1 <= start <= end <= lengths[fields[1]]:
                    raise ValueError(f"{location}: sequence-region outside FASTA")
                output.write(f"##sequence-region {prefix}{fields[1]} {start} {end}\n")
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{location}: expected nine GFF3 columns")
            seqid = fields[0]
            if seqid not in lengths:
                raise ValueError(f"{location}: seqid {seqid!r} absent from FASTA")
            start, end = int(fields[3]), int(fields[4])
            sequence_length = lengths[seqid]
            if start > sequence_length or end > sequence_length:
                attributes_by_key = gff_attributes(fields[8])
                if seqid not in circular or end > 2 * sequence_length:
                    raise ValueError(f"{location}: feature outside FASTA sequence")
                if start > sequence_length:
                    if attributes_by_key.get("exception") != "trans-splicing":
                        raise ValueError(f"{location}: unsupported circular feature coordinates")
                    start -= sequence_length
                    end -= sequence_length
                else:
                    identifier = attributes_by_key.get("ID")
                    if identifier not in wrapped_bounds:
                        raise ValueError(f"{location}: unsupported wrapped circular feature")
                    start, end = wrapped_bounds[identifier]
                fields[3], fields[4] = str(start), str(end)
            if not 1 <= start <= end <= sequence_length:
                raise ValueError(f"{location}: feature outside FASTA sequence")
            attributes = []
            for attribute in fields[8].split(";"):
                key, separator, value = attribute.partition("=")
                if key == "ID" and separator:
                    identifiers.add(value)
                elif key in {"Parent", "Derives_from"} and separator:
                    references.update(value.split(","))
                elif key == "Target" and separator:
                    attribute = f"Target={prefix}{value}"
                attributes.append(attribute)
            fields[0] = prefix + seqid
            fields[8] = rewrite_attributes(";".join(attributes), prefix)
            output.write("\t".join(fields) + "\n")
            features += 1
    if features == 0:
        raise ValueError(f"{path}: no annotation features")
    missing = references - identifiers
    if missing:
        raise ValueError(f"{path}: unresolved GFF3 references: {sorted(missing)[:5]}")


def build_outputs(root: Path, records: list[dict],
                  pairs: list[tuple[Path, Path]]) -> None:
    if len(records) != max(SIZES) or len(pairs) != len(records):
        raise ValueError("Expected one FASTA/GFF3 pair for each of 100 assemblies")
    with tempfile.TemporaryDirectory(prefix=".dataset-build-", dir=root) as temporary:
        staging = Path(temporary)
        with ExitStack() as stack:
            outputs = {}
            for size in SIZES:
                directory = staging / f"{size}-genomes"
                directory.mkdir()
                fasta = stack.enter_context((directory / "genome.fasta").open("w"))
                gff = stack.enter_context((directory / "annotation.gff3").open("w"))
                gff.write("##gff-version 3\n")
                outputs[size] = (fasta, gff)
                save_json(directory / "assemblies.json", records[:size])
                (directory / "accessions.txt").write_text(
                    "".join(record["accession"] + "\n" for record in records[:size]),
                    encoding="utf-8",
                )
            for index, (record, (fasta_path, gff_path)) in enumerate(zip(records, pairs), 1):
                prefix = record["accession"] + "__"
                log(f"[{index}/{len(records)}] {prefix} {record['organism']['organism_name']}")
                lengths = fasta_sequences(fasta_path)
                if any(length == 0 for length in lengths.values()):
                    raise ValueError(f"{fasta_path}: empty FASTA sequence")
                for size, (fasta, gff) in outputs.items():
                    if index <= size:
                        write_fasta(fasta_path, fasta, prefix)
                        fasta.write("\n")
                        write_annotation(gff_path, gff, prefix, lengths)
        # Publish only after every assembly has passed validation.
        for size in SIZES:
            destination = root / f"{size}-genomes"
            destination.mkdir(exist_ok=True)
            for name in ("genome.fasta", "annotation.gff3", "accessions.txt", "assemblies.json"):
                (staging / f"{size}-genomes" / name).replace(destination / name)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--cache-dir", type=Path,
                        help="default: OUTPUT_DIR/.ncbi-datasets; keep to resume/reuse selection")
    parser.add_argument("--datasets", default="datasets", help="NCBI Datasets v2 executable")
    parser.add_argument("--workers", type=int, default=10, help="rehydration workers (1-30)")
    parser.add_argument("--min-busco", type=float, default=0.95,
                        help="minimum annotation BUSCO complete fraction (default: 0.95)")
    parser.add_argument("--select-only", action="store_true",
                        help="freeze/print the selection without downloading sequence data")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace existing concatenated files after successful validation")
    args = parser.parse_args()
    if not 1 <= args.workers <= 30 or not 0 <= args.min_busco <= 1:
        parser.error("--workers must be 1-30 and --min-busco must be 0-1")
    root = args.output_dir.resolve()
    cache = (args.cache_dir or root / ".ncbi-datasets").resolve()
    try:
        if not args.select_only and not args.overwrite:
            for size in SIZES:
                for name in ("genome.fasta", "annotation.gff3", "accessions.txt", "assemblies.json"):
                    path = root / f"{size}-genomes" / name
                    if path.exists():
                        raise ValueError(f"Output already exists: {path}; use --overwrite")
        cli = ensure_cli(args.datasets)
        root.mkdir(parents=True, exist_ok=True)
        cache.mkdir(parents=True, exist_ok=True)
        records = get_selection(cli, cache, args.min_busco)
        total = sum(int(record["assembly_stats"]["total_sequence_length"]) for record in records)
        log(f"Selected {len(records)} species; source genome sequence: {total / 1e9:.2f} Gbp")
        if args.select_only:
            for index, record in enumerate(records, 1):
                print(f"{index}\t{record['accession']}\t{record['organism']['organism_name']}")
            return 0
        pairs = download(cli, cache, records, args.workers)
        build_outputs(root, records, pairs)
    except (OSError, ValueError, KeyError, BadZipFile, subprocess.CalledProcessError) as error:
        log(f"ERROR: {error}")
        return 1
    log(f"Created nested 10/20/100-genome sets in {root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
