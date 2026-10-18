#!/usr/bin/env python3
"""Download RefSeq plant genomes via NCBI FTP or combine local FASTA/GFF3 pairs.

Usage:
    python3 make_dataset.py --output-dir training-data
    python3 make_dataset.py --output-dir training-data --select-only
    python3 make_dataset.py -d local_dir -o combined
"""

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import ExitStack
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile
import time
from typing import TextIO
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import zlib

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
GFF_SUFFIXES = (".gff", ".gff3")
IDENTIFIER_ATTRIBUTES = {
    "ID", "Parent", "Derives_from", "gene_id", "transcript_id",
    "exon_id", "protein_id",
}

SIZES = (10, 20, 100)
LAND_PLANTS = "3193"
DATASET_REPORT_URL = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{LAND_PLANTS}/dataset_report"
ASSEMBLY_SUMMARY_URLS = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/assembly_summary.txt",
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt",
)
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


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


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
                fields = line[1:].split(None, 1)
                if not fields:
                    raise ValueError(f"{path}:{line_number}: empty FASTA sequence ID")
                current = fields[0]
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


def parse_explicit_pairs(values):
    pairs = []
    for value in values:
        fasta, separator, gff = value.partition(":")
        if not separator or not fasta or not gff:
            raise ValueError(f"invalid --pair value {value!r}; use FASTA:GFF3")
        pairs.append((Path(fasta), Path(gff)))
    return pairs


def species_name(record: dict) -> str:
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


def select_assemblies(records: list[dict], min_busco: float, count: int = max(SIZES)) -> list[dict]:
    best: dict[str, dict] = {}
    for record in records:
        if eligible(record, min_busco):
            species = species_name(record)
            if species not in best or quality_key(record) > quality_key(best[species]):
                best[species] = record
    priority = {name: index for index, name in enumerate(PREFERRED_SPECIES)}
    ranked = sorted(best.values(), key=quality_key, reverse=True)
    ranked.sort(key=lambda record: priority.get(species_name(record), len(priority)))
    if len(ranked) < count:
        raise ValueError(
            f"Only {len(ranked)} eligible species; need {count}. "
            "Inspect NCBI metadata or explicitly lower --min-busco; "
            "quality filters will not be relaxed automatically."
        )
    return ranked[:count]


def save_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def fetch_url_bytes(url: str, timeout: int = 120) -> bytes:
    last_error: Exception | None = None
    for attempt in range(1, 5):
        try:
            request = Request(url, headers={"User-Agent": "anno-make-dataset/1.0"})
            with urlopen(request, timeout=timeout) as response:
                return response.read()
        except Exception as error:
            last_error = error
            if attempt < 4:
                time.sleep(min(2 ** attempt, 10))
    raise ValueError(f"Failed to fetch {url}: {last_error}")


def fetch_ncbi_plant_reports() -> list[dict]:
    reports: list[dict] = []
    page_token: str | None = None
    while True:
        params = [
            ("filters.reference_only", "false"),
            ("filters.assembly_source", "refseq"),
            ("filters.has_annotation", "true"),
            ("filters.exclude_atypical", "true"),
            ("filters.assembly_level", "chromosome"),
            ("filters.assembly_level", "complete_genome"),
            ("page_size", "1000"),
        ]
        if page_token:
            params.append(("page_token", page_token))
        url = f"{DATASET_REPORT_URL}?{urlencode(params)}"
        payload = json.loads(fetch_url_bytes(url).decode("utf-8"))
        reports.extend(payload.get("reports", []))
        page_token = payload.get("next_page_token")
        if not page_token:
            break
    return reports


def get_selection(cache: Path, min_busco: float, count: int = max(SIZES)) -> list[dict]:
    selection = cache / "selection.json"
    if selection.exists():
        records = json.loads(selection.read_text(encoding="utf-8"))
        if isinstance(records, list) and len(records) >= count:
            subset = records[:count]
            if (all(eligible(record, min_busco) for record in subset)
                    and len({species_name(record) for record in subset}) == len(subset)):
                log(f"Reusing frozen selection: {selection}")
                return subset
    target_count = max(count, max(SIZES))
    records = select_assemblies(fetch_ncbi_plant_reports(), min_busco, target_count)
    save_json(selection, records)
    return records[:count]


def resolve_ftp_dirs(cache: Path, accessions: list[str]) -> dict[str, str]:
    wanted = set(accessions)
    ftp_dirs: dict[str, str] = {}
    summary_cache = cache / "assembly_summary_plant.txt"
    if summary_cache.is_file() and summary_cache.stat().st_size > 0:
        lines = summary_cache.read_text(encoding="utf-8", errors="replace").splitlines()
    else:
        lines = []
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] in wanted and len(fields) > 19 and fields[19].startswith("http"):
            ftp_dirs[fields[0]] = fields[19].rstrip("/")
    if len(ftp_dirs) == len(wanted):
        return ftp_dirs

    for summary_url in ASSEMBLY_SUMMARY_URLS:
        try:
            text = fetch_url_bytes(summary_url, timeout=180).decode("utf-8", errors="replace")
        except Exception:
            continue
        if summary_url == ASSEMBLY_SUMMARY_URLS[0]:
            temp = summary_cache.with_suffix(".tmp")
            temp.write_text(text, encoding="utf-8")
            temp.replace(summary_cache)
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] in wanted and len(fields) > 19 and fields[19].startswith("http"):
                ftp_dirs[fields[0]] = fields[19].rstrip("/")
        if len(ftp_dirs) == len(wanted):
            return ftp_dirs

    for accession in wanted - ftp_dirs.keys():
        digits = accession.split("_")[1].split(".")[0]
        parent = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
        html = fetch_url_bytes(parent, timeout=60).decode("utf-8", errors="replace")
        matches = re.findall(rf'href="({re.escape(accession)}_[^"]+)/"', html)
        if not matches:
            raise ValueError(f"Could not resolve NCBI FTP directory for {accession}")
        ftp_dirs[accession] = parent + matches[0]
    return ftp_dirs


def parse_md5checksums(text: str) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line in text.splitlines():
        parts = line.strip().split()
        if len(parts) == 2:
            filename = parts[1].lstrip("./")
            checksums[filename] = parts[0].lower()
    return checksums


def stream_ftp_gzip_to_target(url: str, target: Path, expected_gz_md5: str | None) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={"User-Agent": "anno-make-dataset/1.0"})
    with tempfile.NamedTemporaryFile(
        prefix=target.name + ".part-", dir=target.parent, delete=False
    ) as temporary:
        temp_path = Path(temporary.name)
    try:
        gz_hasher = hashlib.md5() if expected_gz_md5 else None
        decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        total_out = 0
        with urlopen(request, timeout=180) as response, temp_path.open("wb") as out:
            while chunk := response.read(1 << 20):
                if gz_hasher is not None:
                    gz_hasher.update(chunk)
                data = decompressor.decompress(chunk)
                if data:
                    out.write(data)
                    total_out += len(data)
            tail = decompressor.flush()
            if tail:
                out.write(tail)
                total_out += len(tail)
        if total_out <= 0:
            raise ValueError(f"Empty decompressed output from {url}")
        if gz_hasher is not None and gz_hasher.hexdigest().lower() != expected_gz_md5:
            raise ValueError(
                f"MD5 mismatch for {url}: got {gz_hasher.hexdigest().lower()}, expected {expected_gz_md5}"
            )
        temp_path.replace(target)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise


def download_assembly_ftp(accession: str, ftp_dir: str, data_root: Path) -> tuple[Path, Path]:
    asm_dir = data_root / accession
    fasta_target = asm_dir / "genomic.fna"
    gff_target = asm_dir / "genomic.gff"
    done_marker = asm_dir / ".complete"

    if not fasta_target.is_file() or fasta_target.stat().st_size == 0:
        candidates = sorted(asm_dir.glob("*_genomic.fna")) if asm_dir.is_dir() else []
        if candidates and candidates[0].stat().st_size > 0:
            fasta_target = candidates[0]
    if not gff_target.is_file() or gff_target.stat().st_size == 0:
        candidates = sorted(asm_dir.glob("*_genomic.gff")) if asm_dir.is_dir() else []
        if candidates and candidates[0].stat().st_size > 0:
            gff_target = candidates[0]

    if (fasta_target.is_file() and fasta_target.stat().st_size > 0
            and gff_target.is_file() and gff_target.stat().st_size > 0):
        return fasta_target, gff_target

    base = ftp_dir.rsplit("/", 1)[-1]
    md5_url = f"{ftp_dir}/md5checksums.txt"
    checksums = parse_md5checksums(fetch_url_bytes(md5_url, timeout=60).decode("utf-8", errors="replace"))

    fna_gz_name = f"{base}_genomic.fna.gz"
    gff_gz_name = f"{base}_genomic.gff.gz"
    targets = [
        (f"{ftp_dir}/{fna_gz_name}", fasta_target, checksums.get(fna_gz_name)),
        (f"{ftp_dir}/{gff_gz_name}", gff_target, checksums.get(gff_gz_name)),
    ]
    for url, target, expected_md5 in targets:
        if target.is_file() and target.stat().st_size > 0 and done_marker.is_file():
            continue
        last_error: Exception | None = None
        for attempt in range(1, 5):
            try:
                log(f"FTP download ({attempt}/4): {url} -> {accession}/{target.name}")
                stream_ftp_gzip_to_target(url, target, expected_md5)
                break
            except Exception as error:
                last_error = error
                if attempt < 4:
                    time.sleep(min(2 ** attempt, 10))
        else:
            raise ValueError(f"Failed to download {url}: {last_error}")
    done_marker.touch()
    return fasta_target, gff_target


def download(cache: Path, records: list[dict], workers: int) -> list[tuple[Path, Path]]:
    accessions = cache / "accessions.txt"
    accessions.write_text("".join(record["accession"] + "\n" for record in records), encoding="utf-8")
    data_root = cache / "package" / "ncbi_dataset" / "data"
    data_root.mkdir(parents=True, exist_ok=True)

    ftp_dirs = resolve_ftp_dirs(cache, [record["accession"] for record in records])
    pool_workers = max(1, min(workers, 8, len(records)))
    results: dict[str, tuple[Path, Path]] = {}
    with ThreadPoolExecutor(max_workers=pool_workers) as executor:
        futures = {
            executor.submit(
                download_assembly_ftp,
                record["accession"],
                ftp_dirs[record["accession"]],
                data_root,
            ): record["accession"]
            for record in records
        }
        for future in as_completed(futures):
            acc = futures[future]
            results[acc] = future.result()
    return [results[record["accession"]] for record in records]


def gff_attributes(value: str) -> dict[str, str]:
    attributes = {}
    for item in value.split(";"):
        key, separator, attribute_value = item.partition("=")
        if separator:
            attributes[key] = attribute_value
    return attributes


def normalize_circular_coord(coord: int, length: int) -> int:
    return ((coord - 1) % length) + 1


def normalize_circular_interval(start: int, end: int, length: int) -> tuple[int, int] | None:
    if length <= 0 or start < 1 or end < 1 or start > 2 * length or end > 2 * length:
        return None
    if 1 <= start <= end <= length:
        return start, end
    if start > length and end > length and start <= end:
        return start - length, end - length
    s = normalize_circular_coord(start, length)
    e = normalize_circular_coord(end, length)
    return (s, e) if s <= e else (e, s)


def circular_trans_splice_bounds(
    path: Path, lengths: dict[str, int]
) -> tuple[set[str], dict[str, tuple[int, int]],
           dict[tuple[str, int, int, str], tuple[int, int]]]:
    circular: set[str] = set()
    feature_seqid: dict[str, str] = {}
    children: dict[str, list[tuple[str, int, int]]] = {}
    with open_text(path) as source:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[0] not in lengths:
                continue
            seqid = fields[0]
            attributes = gff_attributes(fields[8])
            if fields[2] == "region" and attributes.get("Is_circular", "").casefold() == "true":
                circular.add(seqid)
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            if start < 1 or end < 1 or start > end or end > lengths[seqid]:
                if start <= 2 * lengths[seqid] and end <= 2 * lengths[seqid]:
                    circular.add(seqid)
            identifier = attributes.get("ID")
            if identifier:
                feature_seqid[identifier] = seqid
            parents = [parent for parent in attributes.get("Parent", "").split(",") if parent]
            for parent in parents:
                children.setdefault(parent, []).append((seqid, start, end))

    adjusted: dict[str, tuple[int, int]] = {}
    adjusted_children: dict[tuple[str, int, int, str], tuple[int, int]] = {}
    for parent_id, child_rows in children.items():
        seqid = feature_seqid.get(parent_id) or (child_rows[0][0] if child_rows else "")
        if seqid not in circular or seqid not in lengths:
            continue
        length = lengths[seqid]
        segments: list[tuple[int, int]] = []
        for child_seqid, start, end in child_rows:
            if child_seqid != seqid:
                continue
            normalized = normalize_circular_interval(start, end, length)
            if normalized is None:
                continue
            segments.append(normalized)
            adjusted_children[(seqid, start, end, parent_id)] = normalized
        if segments:
            adjusted[parent_id] = (
                min(s for s, _ in segments),
                max(e for _, e in segments),
            )
    return circular, adjusted, adjusted_children


def write_annotation(path: Path, output: TextIO, prefix: str,
                     lengths: dict[str, int]) -> None:
    identifiers: set[str] = set()
    references: set[str] = set()
    features = 0
    circular, wrapped_bounds, wrapped_children = circular_trans_splice_bounds(path, lengths)
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
            if start < 1 or start > end or end > sequence_length:
                attributes_by_key = gff_attributes(fields[8])
                if seqid not in circular or start < 1 or end < 1 or start > 2 * sequence_length or end > 2 * sequence_length:
                    raise ValueError(f"{location}: feature outside FASTA sequence")
                identifier = attributes_by_key.get("ID")
                normalized = None
                if identifier and identifier in wrapped_bounds:
                    normalized = wrapped_bounds[identifier]
                if normalized is None:
                    normalized = next(
                        (wrapped_children[(seqid, start, end, parent)]
                         for parent in attributes_by_key.get("Parent", "").split(",")
                         if (seqid, start, end, parent) in wrapped_children),
                        None,
                    )
                if normalized is None:
                    normalized = normalize_circular_interval(start, end, sequence_length)
                if normalized is None:
                    raise ValueError(f"{location}: unsupported circular feature coordinates")
                start, end = normalized
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
                  pairs: list[tuple[Path, Path]],
                  sizes: tuple[int, ...] = SIZES) -> None:
    if len(pairs) != len(records):
        raise ValueError("Expected one FASTA/GFF3 pair for each selected assembly")
    active_sizes = tuple(sorted({s for s in sizes if s <= len(records)}))
    if not active_sizes:
        active_sizes = (len(records),)
    with tempfile.TemporaryDirectory(prefix=".dataset-build-", dir=root) as temporary:
        staging = Path(temporary)
        with ExitStack() as stack:
            outputs = {}
            for size in active_sizes:
                directory = staging / f"{size}-genomes"
                directory.mkdir()
                fasta = stack.enter_context((directory / "genome.fasta").open("w", encoding="utf-8"))
                gff = stack.enter_context((directory / "annotation.gff3").open("w", encoding="utf-8"))
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
        for size in active_sizes:
            destination = root / f"{size}-genomes"
            destination.mkdir(exist_ok=True)
            for name in ("genome.fasta", "annotation.gff3", "accessions.txt", "assemblies.json"):
                (staging / f"{size}-genomes" / name).replace(destination / name)


def run_local_combine(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
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
        with open(output_fasta, "w", encoding="utf-8") as fasta_out, open(output_gff, "w", encoding="utf-8") as gff_out:
            gff_out.write("##gff-version 3\n")
            used = set()
            for fasta, gff in pairs:
                prefix = prefix_for(fasta)
                if prefix in used:
                    raise ValueError(f"duplicate dataset prefix {prefix}; use explicit unique input names")
                used.add(prefix)
                lengths = fasta_sequences(fasta)
                write_fasta(fasta, fasta_out, prefix)
                write_annotation(gff, gff_out, prefix, lengths)
                prefixes.append((fasta.name, gff.name, prefix))
    except (OSError, ValueError) as error:
        parser.error(str(error))

    print(f"FASTA: {output_fasta}")
    print(f"GFF3:  {output_gff}")
    for fasta, gff, prefix in prefixes:
        print(f"{prefix}: {fasta} + {gff}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-d", "--directory", type=Path,
                        help="local directory containing FASTA and GFF/GFF3 files to combine")
    parser.add_argument("--pair", action="append", default=[], metavar="FASTA:GFF3",
                        help="explicit local input pair; repeat for multiple datasets")
    parser.add_argument("--fasta", nargs="+", type=Path,
                        help="explicit local FASTA files, paired by order with --gff")
    parser.add_argument("--gff", nargs="+", type=Path,
                        help="explicit local GFF/GFF3 files, paired by order with --fasta")
    parser.add_argument("-o", "--output-prefix", default="combined",
                        help="output basename inside -d (default: combined)")
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="output directory for NCBI FTP plant datasets (default: .)")
    parser.add_argument("--cache-dir", type=Path,
                        help="default: OUTPUT_DIR/.ncbi-datasets; keep to resume/reuse selection")
    parser.add_argument("-n", "--num-genomes", type=int, default=max(SIZES),
                        help=f"number of plant genomes to download (default: {max(SIZES)})")
    parser.add_argument("--workers", type=int, default=10, help="download workers (1-30)")
    parser.add_argument("--min-busco", type=float, default=0.95,
                        help="minimum annotation BUSCO complete fraction (default: 0.95)")
    parser.add_argument("--select-only", action="store_true",
                        help="freeze/print the selection without downloading sequence data")
    parser.add_argument("--download-only", action="store_true",
                        help="download per-genome FASTA/GFF3 via FTP without building concatenated sets")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace existing concatenated files after successful validation")
    args = parser.parse_args()

    if args.directory is not None or args.pair or args.fasta or args.gff:
        if args.directory is None:
            args.directory = Path(".")
        return run_local_combine(args, parser)

    if not 1 <= args.workers <= 30 or not 0 <= args.min_busco <= 1 or args.num_genomes < 1:
        parser.error("--workers must be 1-30, --min-busco must be 0-1, and -n must be >= 1")
    root = args.output_dir.resolve()
    cache = (args.cache_dir or root / ".ncbi-datasets").resolve()
    active_sizes = tuple(s for s in SIZES if s <= args.num_genomes) or (args.num_genomes,)
    try:
        if not args.select_only and not args.download_only and not args.overwrite:
            for size in active_sizes:
                for name in ("genome.fasta", "annotation.gff3", "accessions.txt", "assemblies.json"):
                    path = root / f"{size}-genomes" / name
                    if path.exists():
                        raise ValueError(f"Output already exists: {path}; use --overwrite")
        root.mkdir(parents=True, exist_ok=True)
        cache.mkdir(parents=True, exist_ok=True)
        records = get_selection(cache, args.min_busco, args.num_genomes)
        total = sum(int(record["assembly_stats"]["total_sequence_length"]) for record in records)
        log(f"Selected {len(records)} species; source genome sequence: {total / 1e9:.2f} Gbp")
        if args.select_only:
            for index, record in enumerate(records, 1):
                print(f"{index}\t{record['accession']}\t{record['organism']['organism_name']}")
            return 0
        pairs = download(cache, records, args.workers)
        if not args.download_only:
            build_outputs(root, records, pairs, active_sizes)
    except (OSError, ValueError) as error:
        log(f"ERROR: {error}")
        return 1
    log("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
