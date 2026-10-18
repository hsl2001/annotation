#!/usr/bin/env python3
"""Protein gene-level and gffcompare CDS exon-level accuracy."""

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import gzip
import math
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np


BASE_TO_INT = np.full(256, 4, dtype=np.uint8)
for _ch, _val in zip(b"ACGTacgt", (0, 1, 2, 3, 0, 1, 2, 3)):
    BASE_TO_INT[_ch] = _val
COMPLEMENT = np.array([3, 2, 1, 0, 4], dtype=np.uint8)

STOP_AA = ord("*")
CODON_TABLE = np.full(125, ord("X"), dtype=np.uint8)
_STANDARD_CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_BMAP = {"A": 0, "C": 1, "G": 2, "T": 3}
for _codon, _aa in _STANDARD_CODONS.items():
    _idx = _BMAP[_codon[0]] * 25 + _BMAP[_codon[1]] * 5 + _BMAP[_codon[2]]
    CODON_TABLE[_idx] = ord(_aa)
ATG_CODE = _BMAP["A"] * 25 + _BMAP["T"] * 5 + _BMAP["G"]
STOP_CODES = np.array([
    _BMAP[c[0]] * 25 + _BMAP[c[1]] * 5 + _BMAP[c[2]]
    for c in ("TAA", "TAG", "TGA")
], dtype=np.uint8)

BLOSUM62_DIAG = np.zeros(256, dtype=np.float64)
for _aa_ch, _diag in {
    "A": 4, "R": 5, "N": 6, "D": 6, "C": 9, "Q": 5, "E": 5, "G": 6, "H": 8, "I": 4,
    "L": 4, "K": 5, "M": 5, "F": 6, "P": 7, "S": 4, "T": 5, "W": 11, "Y": 7, "V": 4,
}.items():
    BLOSUM62_DIAG[ord(_aa_ch)] = _diag
_BLAST_LAMBDA = 0.267
_BLAST_K = 0.041
_LN2 = math.log(2.0)
_LN_K = math.log(_BLAST_K)


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def attributes(text):
    result = {}
    for field in text.split(";"):
        field = field.strip()
        if not field:
            continue
        key, separator, value = field.partition("=")
        if not separator:
            parts = field.split(None, 1)
            if len(parts) == 2:
                key, value = parts
            else:
                continue
        result[key] = value.strip('"')
    return result


def gff_records(path, wanted_features=None):
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.rstrip() == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected GFF3/GTF with 9 columns")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 1 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            kind = fields[2].lower()
            if kind in ("exon", "cds") and fields[6] not in ("+", "-"):
                raise ValueError(f"{path}:{line_number}: invalid strand")
            if wanted_features is not None and kind not in wanted_features:
                continue
            yield fields, attributes(fields[8]), start, end


def feature_parents(info):
    return [parent for parent in (info.get("Parent") or info.get("transcript_id") or info.get("gene_id") or "").split(",")
            if parent]


def transcript_intervals(path, feature):
    metadata, intervals = {}, defaultdict(list)
    wanted_features = {"mrna", "transcript", feature}
    for fields, info, start, end in gff_records(path, wanted_features):
        kind = fields[2].lower()
        identifier = info.get("ID") or info.get("transcript_id")
        if kind in ("mrna", "transcript") and identifier:
            gene = (info.get("Parent") or info.get("gene_id") or identifier).split(",")[0]
            metadata[identifier] = (gene or identifier, fields[0], fields[6])
        elif kind == feature:
            interval = (fields[0], fields[6], start, end)
            for parent in feature_parents(info):
                intervals[parent].append(interval)
                if info.get("gene_id"):
                    metadata.setdefault(parent, (info["gene_id"], fields[0], fields[6]))
    models = {}
    for transcript, segments in intervals.items():
        seqid, strand = segments[0][:2]
        gene, seqid, strand = metadata.get(transcript, (transcript, seqid, strand))
        if all(segment[:2] == (seqid, strand) for segment in segments):
            deduped = sorted(set(segments), key=lambda item: (item[2], item[3]))
            models[transcript] = gene, seqid, strand, deduped
    return models


def encode_genome(genome):
    """Convert genome dict (seqid -> str or np.ndarray) into encoded uint8 NumPy arrays."""
    encoded = {}
    for name, seq in genome.items():
        if isinstance(seq, np.ndarray) and seq.dtype == np.uint8:
            encoded[name] = seq
        else:
            raw = seq.encode("ascii") if isinstance(seq, str) else bytes(seq)
            encoded[name] = BASE_TO_INT[np.frombuffer(raw, dtype=np.uint8)]
    return encoded


def select_longest_transcripts(transcripts):
    """Select one representative transcript per (seqid, strand, gene) with consistent tie-breaking."""
    selected = {}
    for transcript, record in transcripts.items():
        key = (record["seqid"], record["strand"], record["gene"])
        previous = selected.get(key)
        if previous is None or (record["length"], transcript) > (previous[1]["length"], previous[0]):
            selected[key] = (transcript, record)
    return {transcript: record for transcript, record in selected.values()}


def read_genes(path, genome=None, validate=False):
    """Map each coding gene to its longest-CDS transcript and CDS-spanning region.

    Returns gene_id -> (seqid, start, end, strand, transcript_id).
    """
    if genome is not None:
        transcripts = select_longest_transcripts(coding_transcripts(path, genome, validate=validate))
        return {
            record["gene"]: (
                record["seqid"],
                min(interval[0] for interval in record["intervals"]),
                max(interval[1] for interval in record["intervals"]),
                record["strand"],
                transcript,
            )
            for transcript, record in transcripts.items()
        }
    transcripts = transcript_intervals(path, "cds")
    genes = {}
    for transcript, (gene, seqid, strand, intervals) in transcripts.items():
        coords = np.asarray([(s, e) for _, _, s, e in intervals], dtype=np.int64)
        length = int(np.sum(coords[:, 1] - coords[:, 0] + 1))
        start = int(np.min(coords[:, 0]))
        end = int(np.max(coords[:, 1]))
        previous = genes.get(gene)
        if previous is None or (length, transcript) > (previous[0], previous[5]):
            genes[gene] = (length, seqid, start, end, strand, transcript)
    return {gene: value[1:] for gene, value in genes.items()}


def read_genome(path):
    sequences = {name: sequence.upper() for name, sequence in read_fasta(path)}
    if not sequences:
        raise ValueError(f"{path}: no FASTA sequences found")
    return sequences


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def coding_transcripts(path, genome, validate=True):
    """Read coding transcripts and optionally apply the paper's invalid-model filters using NumPy."""
    transcripts = transcript_intervals(path, "cds")
    encoded_genome = encode_genome(genome)
    valid = {}
    for transcript, (gene, seqid, strand, intervals) in transcripts.items():
        if seqid not in encoded_genome:
            raise ValueError(f"{path}: sequence {seqid} is missing from the genome FASTA")
        seq_arr = encoded_genome[seqid]
        coords = np.asarray([(s, e) for _, _, s, e in intervals], dtype=np.int64)
        if coords.shape[0] > 1:
            coords = coords[np.argsort(coords[:, 0], kind="mergesort")]
        if int(np.max(coords[:, 1])) > seq_arr.size:
            raise ValueError(f"{path}: CDS extends beyond sequence {seqid}")
        if validate and coords.shape[0] > 1:
            introns = coords[1:, 0] - coords[:-1, 1] - 1
            if np.any(introns <= 1):
                continue
        cds = np.concatenate([seq_arr[s - 1:e] for s, e in coords])
        if strand == "-":
            cds = COMPLEMENT[cds[::-1]]
        coding_length = int(cds.size)
        n_codons = coding_length // 3
        if validate:
            if coding_length < 6 or coding_length % 3 != 0:
                continue
            codons = (
                cds[0::3] * np.uint8(25)
                + cds[1::3] * np.uint8(5)
                + cds[2::3]
            )
            if codons[0] != ATG_CODE or int(codons[-1]) not in (int(STOP_CODES[0]), int(STOP_CODES[1]), int(STOP_CODES[2])):
                continue
            aa = CODON_TABLE[codons]
            if n_codons > 2 and np.any(aa[1:-1] == STOP_AA):
                continue
            protein = aa[:-1].tobytes().decode("ascii")
        else:
            if n_codons > 0:
                trimmed = cds[:n_codons * 3]
                codons = (
                    trimmed[0::3] * np.uint8(25)
                    + trimmed[1::3] * np.uint8(5)
                    + trimmed[2::3]
                )
                aa = CODON_TABLE[codons]
                protein = aa[aa != STOP_AA].tobytes().decode("ascii")
            else:
                protein = ""
        valid[transcript] = {
            "gene": gene,
            "seqid": seqid,
            "strand": strand,
            "intervals": [(int(s), int(e)) for s, e in coords],
            "length": coding_length,
            "protein": protein,
        }
    return valid


def filtered_gff(path, genome, output, longest_only=False, validate=None, transcripts=None):
    if validate is None:
        validate = longest_only
    if transcripts is None:
        transcripts = coding_transcripts(path, genome, validate=validate)
    if longest_only:
        transcripts = select_longest_transcripts(transcripts)

    genes = {}
    for record in transcripts.values():
        key = record["seqid"], record["strand"], record["gene"]
        start = min(interval[0] for interval in record["intervals"])
        end = max(interval[1] for interval in record["intervals"])
        previous = genes.get(key, (start, end))
        genes[key] = min(start, previous[0]), max(end, previous[1])
    written = set()
    with open(output, "w") as stream:
        stream.write("##gff-version 3\n")
        for transcript, record in sorted(transcripts.items()):
            start = min(interval[0] for interval in record["intervals"])
            end = max(interval[1] for interval in record["intervals"])
            gene = record["gene"]
            gene_key = (record["seqid"], record["strand"], gene)
            if gene_key not in written:
                written.add(gene_key)
                gene_start, gene_end = genes[gene_key]
                stream.write(
                    f"{record['seqid']}\taccuracy\tgene\t{gene_start}\t{gene_end}\t.\t"
                    f"{record['strand']}\t.\tID={gene}\n"
                )
            stream.write(
                f"{record['seqid']}\taccuracy\ttranscript\t{start}\t{end}\t.\t"
                f"{record['strand']}\t.\tID={transcript};Parent={gene}\n"
            )
            for number, (exon_start, exon_end) in enumerate(sorted(record["intervals"]), 1):
                stream.write(
                    f"{record['seqid']}\taccuracy\texon\t{exon_start}\t{exon_end}\t.\t"
                    f"{record['strand']}\t.\tID={transcript}.exon{number};Parent={transcript}\n"
                )
    return len(transcripts)


def parse_gffcompare_stats(path):
    pattern = re.compile(
        r"^\s*Exon level:\s*"
        r"([0-9.]+)\s*\|\s*([0-9.]+)"
    )
    with open(path) as stream:
        for line in stream:
            match = pattern.match(line)
            if not match:
                continue
            recall, precision = float(match.group(1)) / 100, float(match.group(2)) / 100
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
            return {"exon": (precision, recall, f1)}
    raise ValueError(f"could not parse gffcompare exon statistics: {path}")


def run_gffcompare(reference, query, workdir, label):
    prefix = workdir / label
    subprocess.run(
        ["gffcompare", "--no-exon-merge", "--strict-match", "-e", "0",
         "-r", str(reference), "-o", str(prefix), str(query)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    stats_path = pathlib.Path(f"{prefix}.stats")
    if not stats_path.exists() and prefix.is_file():
        stats_path = prefix
    return parse_gffcompare_stats(stats_path)


def scores(tp, fp, fn):
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, recall, f1


def genome_fasta(genome, workdir):
    """Return a plain, gffread-indexable FASTA path (decompressing .gz if needed)."""
    if str(genome).endswith(".gz"):
        with gzip.GzipFile(filename=genome, mode="rb") as source, \
                tempfile.NamedTemporaryFile(mode="wb", prefix="genome-", suffix=".fa", dir=workdir, delete=False) as target:
            shutil.copyfileobj(source, target)
            return pathlib.Path(target.name)
    link = workdir / "genome.fa"
    source = pathlib.Path(genome).resolve()
    if pathlib.Path(genome).absolute() == link.absolute():
        return link
    if link.is_symlink():
        (workdir / "genome.fa.fai").unlink(missing_ok=True)
    if link.resolve() == source:
        return link
    if link.is_symlink():
        link.unlink()
    if link.exists():
        raise ValueError(f"Refusing to overwrite existing FASTA: {link}")
    link.symlink_to(source)
    return link


def write_protein_gff(path, genes, output):
    selected_genes = set(genes)
    selected_transcripts = {value[4] for value in genes.values()}
    with open(output, "w") as target:
        target.write("##gff-version 3\n")
        for fields, info, start, end in gff_records(path):
            feature = fields[2].lower()
            identifier = info.get("ID") or info.get("transcript_id")
            parents = [parent for parent in feature_parents(info) if parent in selected_transcripts]
            keep = ((feature == "gene" and identifier in selected_genes) or
                    (feature in ("mrna", "transcript") and identifier in selected_transcripts) or
                    (feature == "cds" and bool(parents)))
            if not keep:
                continue
            if fields[6] not in ("+", "-"):
                raise ValueError(f"{path}: selected coding feature has invalid strand")
            if feature == "cds" and "Parent" in info:
                info["Parent"] = ",".join(parents)
                fields[8] = ";".join(f"{key}={value}" for key, value in info.items())
            target.write("\t".join(fields) + "\n")


def read_fasta(path):
    name, chunks = None, []
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                header = line[1:].split()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name, chunks = header[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"{path}:{line_number}: FASTA sequence before header")
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def extract_proteins(gff, genome, genes, out_fasta, workdir, tag, transcripts=None):
    """Translate the longest-CDS transcript of each gene using NumPy, header renamed to gene id."""
    if not genes:
        out_fasta.write_text("")
        return 0
    if transcripts is None:
        if isinstance(genome, (str, pathlib.Path)):
            genome = read_genome(genome)
        transcripts = coding_transcripts(gff, genome, validate=False)
    written = 0
    with open(out_fasta, "w") as target:
        for gene, (_, _, _, _, transcript) in genes.items():
            record = transcripts.get(transcript)
            if record is None:
                continue
            sequence = record["protein"].replace(".", "").replace("*", "")
            if not sequence:
                continue
            target.write(f">{gene}\n{sequence}\n")
            written += 1
    return written


def _union_length(intervals):
    """Compute merged length of 1-based inclusive intervals [(start, end), ...] using NumPy."""
    if len(intervals) == 1:
        s, e = intervals[0]
        return e - s + 1
    arr = np.asarray(intervals, dtype=np.int64)
    arr = arr[np.argsort(arr[:, 0], kind="mergesort")]
    starts = arr[:, 0]
    ends = arr[:, 1]
    max_prev_end = np.maximum.accumulate(ends[:-1])
    new_seg = starts[1:] > (max_prev_end + 1)
    seg_starts = np.empty(int(np.sum(new_seg)) + 1, dtype=np.int64)
    seg_ends = np.empty_like(seg_starts)
    seg_starts[0] = starts[0]
    seg_starts[1:] = starts[1:][new_seg]
    seg_ends[:-1] = max_prev_end[new_seg]
    seg_ends[-1] = max(int(max_prev_end[-1]), int(ends[-1]))
    return int(np.sum(seg_ends - seg_starts + 1))


def _parse_blast_tabular(stdout):
    """Aggregate HSPs per (query, subject) with merged query/subject coverage."""
    grouped = {}
    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 10:
            continue
        query, subject = parts[0], parts[1]
        hit_score, hit_evalue = float(parts[2]), float(parts[3])
        qstart, qend, qlen = int(parts[4]), int(parts[5]), int(parts[6])
        sstart, send, slen = int(parts[7]), int(parts[8]), int(parts[9])
        if qstart > qend:
            qstart, qend = qend, qstart
        if sstart > send:
            sstart, send = send, sstart
        key = (query, subject)
        entry = grouped.get(key)
        if entry is None:
            grouped[key] = [hit_score, hit_evalue, qlen, slen, [(qstart, qend)], [(sstart, send)]]
        else:
            if hit_score > entry[0]:
                entry[0] = hit_score
            if hit_evalue < entry[1]:
                entry[1] = hit_evalue
            entry[4].append((qstart, qend))
            entry[5].append((sstart, send))
    hits = []
    for (query, subject), (best_score, min_evalue, qlen, slen, q_ivs, s_ivs) in grouped.items():
        qcov = _union_length(q_ivs) / qlen if qlen else 0.0
        scov = _union_length(s_ivs) / slen if slen else 0.0
        hits.append((query, subject, best_score, min_evalue, qcov, scov))
    return hits


def find_overlapping_pairs(query_genes, ref_genes):
    """Find all (query, subject) pairs on the same (seqid, strand) with overlapping spans using NumPy."""
    q_by_key = defaultdict(list)
    r_by_key = defaultdict(list)
    for q_gene, (seqid, start, end, strand, *_) in query_genes.items():
        q_by_key[(seqid, strand)].append((start, end, q_gene))
    for r_gene, (seqid, start, end, strand, *_) in ref_genes.items():
        r_by_key[(seqid, strand)].append((start, end, r_gene))

    pairs = set()
    components = []
    for key, q_list in q_by_key.items():
        r_list = r_by_key.get(key)
        if not r_list:
            continue
        q_list.sort(key=lambda item: (item[0], item[1]))
        r_list.sort(key=lambda item: (item[0], item[1]))
        q_coords = np.asarray([(s, e) for s, e, _ in q_list], dtype=np.int64)
        r_coords = np.asarray([(s, e) for s, e, _ in r_list], dtype=np.int64)
        overlap_mask = (q_coords[:, 0:1] <= r_coords[None, :, 1]) & (r_coords[None, :, 0] <= q_coords[:, 1:2])
        q_idx, r_idx = np.nonzero(overlap_mask)
        if q_idx.size == 0:
            continue
        q_adj = defaultdict(list)
        r_adj = defaultdict(list)
        for qi, ri in zip(q_idx.tolist(), r_idx.tolist()):
            q_name = q_list[qi][2]
            r_name = r_list[ri][2]
            pairs.add((q_name, r_name))
            q_adj[qi].append(ri)
            r_adj[ri].append(qi)
        visited_q = set()
        visited_r = set()
        for start_qi in q_adj:
            if start_qi in visited_q:
                continue
            comp_q, comp_r = [], []
            stack_q = [start_qi]
            visited_q.add(start_qi)
            while stack_q:
                curr_q = stack_q.pop()
                comp_q.append(q_list[curr_q][2])
                for nxt_r in q_adj[curr_q]:
                    if nxt_r not in visited_r:
                        visited_r.add(nxt_r)
                        comp_r.append(r_list[nxt_r][2])
                        for nxt_q in r_adj[nxt_r]:
                            if nxt_q not in visited_q:
                                visited_q.add(nxt_q)
                                stack_q.append(nxt_q)
            components.append((comp_q, comp_r))
    return pairs, components


def _self_blast_stats(seq, total_ref_aa):
    """Return (bitscore, evalue) for an exact protein self-match using NumPy BLOSUM62."""
    arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    raw = float(np.sum(BLOSUM62_DIAG[arr]))
    bitscore = (_BLAST_LAMBDA * raw - _LN_K) / _LN2
    evalue = _BLAST_K * arr.size * total_ref_aa * math.exp(-_BLAST_LAMBDA * raw)
    return bitscore, evalue


def blast_hits(query_fasta, ref_fasta, workdir, evalue, threads, components=None, bitscore=50.0):
    columns = "6 qseqid sseqid bitscore evalue qstart qend qlen sstart send slen"
    if components is not None:
        if not components:
            return []
        q_seqs = dict(read_fasta(query_fasta))
        r_seqs = dict(read_fasta(ref_fasta))
        total_ref_aa = sum(len(seq) for seq in r_seqs.values())
        if total_ref_aa == 0:
            return []
        exact_hits = []
        batches = []
        cur_q, cur_r = [], []
        for comp_q, comp_r in components:
            q_valid = [q for q in comp_q if q in q_seqs]
            r_valid = [r for r in comp_r if r in r_seqs]
            if not q_valid or not r_valid:
                continue
            if len(q_valid) == 1 and len(r_valid) == 1:
                q_id, r_id = q_valid[0], r_valid[0]
                seq = q_seqs[q_id]
                if seq == r_seqs[r_id] and len(seq) >= 35 and "X" not in seq:
                    b_score, e_val = _self_blast_stats(seq, total_ref_aa)
                    if b_score >= bitscore + 10.0 and e_val <= evalue * 1e-3:
                        exact_hits.append((q_id, r_id, b_score, e_val, 1.0, 1.0))
                        continue
            if cur_q and (len(cur_q) + len(q_valid) > 35 or len(cur_r) + len(r_valid) > 35):
                batches.append((cur_q, cur_r))
                cur_q, cur_r = [], []
            cur_q.extend(q_valid)
            cur_r.extend(r_valid)
        if cur_q:
            batches.append((cur_q, cur_r))
        if not batches:
            return exact_hits

        def run_batch(item):
            idx, (batch_q, batch_r) = item
            q_path = workdir / f"batch_{idx}.q.faa"
            r_path = workdir / f"batch_{idx}.r.faa"
            q_path.write_text("".join(f">{q}\n{q_seqs[q]}\n" for q in batch_q))
            r_path.write_text("".join(f">{r}\n{r_seqs[r]}\n" for r in batch_r))
            res = subprocess.run(
                ["blastp", "-query", str(q_path), "-subject", str(r_path),
                 "-dbsize", str(total_ref_aa), "-evalue", str(evalue),
                 "-outfmt", columns, "-max_target_seqs", str(max(500, len(batch_r)))],
                check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
            )
            q_path.unlink(missing_ok=True)
            r_path.unlink(missing_ok=True)
            return res.stdout

        with ThreadPoolExecutor(max_workers=max(1, threads)) as executor:
            outputs = list(executor.map(run_batch, enumerate(batches)))
        return exact_hits + _parse_blast_tabular("\n".join(outputs))

    database = workdir / "refdb"
    subprocess.run(["makeblastdb", "-in", str(ref_fasta), "-dbtype", "prot", "-out", str(database)],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run(
        ["blastp", "-query", str(query_fasta), "-db", str(database), "-evalue", str(evalue),
         "-outfmt", columns, "-max_target_seqs", "500", "-num_threads", str(threads)],
        check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    return _parse_blast_tabular(result.stdout)


def region_overlap(a, b):
    return a[0] == b[0] and a[1] <= b[2] and b[1] <= a[2]


def gene_level(reference_gff, query_gff, genome, workdir, *, bitscore, evalue, coverage, threads,
               ref_transcripts=None, query_transcripts=None):
    encoded = encode_genome(genome) if isinstance(genome, dict) else encode_genome(read_genome(genome))
    if ref_transcripts is None:
        ref_transcripts = coding_transcripts(reference_gff, encoded, validate=True)
    if query_transcripts is None:
        query_transcripts = coding_transcripts(query_gff, encoded, validate=False)
    ref_tx = select_longest_transcripts(ref_transcripts)
    query_tx = select_longest_transcripts(query_transcripts)
    ref_genes = {
        rec["gene"]: (
            rec["seqid"],
            min(iv[0] for iv in rec["intervals"]),
            max(iv[1] for iv in rec["intervals"]),
            rec["strand"],
            tid,
        )
        for tid, rec in ref_tx.items()
    }
    query_genes = {
        rec["gene"]: (
            rec["seqid"],
            min(iv[0] for iv in rec["intervals"]),
            max(iv[1] for iv in rec["intervals"]),
            rec["strand"],
            tid,
        )
        for tid, rec in query_tx.items()
    }
    ref_fasta, query_fasta = workdir / "reference.faa", workdir / "query.faa"
    positives = extract_proteins(reference_gff, encoded, ref_genes, ref_fasta, workdir, "reference", transcripts=ref_tx)
    predicted = extract_proteins(query_gff, encoded, query_genes, query_fasta, workdir, "query", transcripts=query_tx)
    if positives == 0 or predicted == 0:
        return scores(0, predicted, positives), (0, positives, predicted)

    valid_pairs, components = find_overlapping_pairs(query_genes, ref_genes)
    if not valid_pairs:
        return scores(0, predicted, positives), (0, positives, predicted)

    candidates = []
    for query, subject, hit_score, hit_evalue, qcov, scov in blast_hits(
        query_fasta, ref_fasta, workdir, evalue, threads, components=components, bitscore=bitscore
    ):
        if (query, subject) not in valid_pairs:
            continue
        if hit_score < bitscore or hit_evalue > evalue:
            continue
        if qcov < coverage or scov < coverage:
            continue
        candidates.append((hit_score, query, subject))

    candidates.sort(reverse=True)
    used_query, used_ref, tp = set(), set(), 0
    for _, query, subject in candidates:
        if query in used_query or subject in used_ref:
            continue
        used_query.add(query)
        used_ref.add(subject)
        tp += 1
    return scores(tp, predicted - tp, positives - tp), (tp, positives, predicted)


def main():
    parser = argparse.ArgumentParser(description="Protein gene-level and gffcompare CDS exon-level accuracy")
    parser.add_argument("-r", "--reference", required=True, help="reference GFF3/GTF(.gz)")
    parser.add_argument("-q", "--query", required=True, help="query GFF3/GTF(.gz)")
    parser.add_argument("-g", "--genome", required=True, help="genome FASTA(.gz)")
    parser.add_argument("--bitscore", type=float, default=50.0, help="minimum blastp bit score")
    parser.add_argument("--evalue", type=float, default=1e-5, help="maximum blastp e-value")
    parser.add_argument("--coverage", type=float, default=0.70,
                        help="minimum best-HSP coverage of both proteins (fraction)")
    parser.add_argument("--threads", type=int, default=4, help="blastp threads")
    parser.add_argument("--workdir", type=pathlib.Path, help="keep protein/blast intermediates here")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    if not math.isfinite(args.bitscore) or args.bitscore < 0:
        parser.error("--bitscore must be finite and nonnegative")
    if not math.isfinite(args.evalue) or args.evalue <= 0:
        parser.error("--evalue must be finite and positive")
    if not 0 <= args.coverage < 1:
        parser.error("--coverage must be in [0, 1)")

    keep = args.workdir is not None
    if keep:
        workdir = args.workdir
    else:
        base = "tmp" if pathlib.Path("tmp").is_dir() else "."
        workdir = pathlib.Path(tempfile.mkdtemp(prefix="anno-gene-", dir=base))
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        genome = encode_genome(read_genome(args.genome))
        ref_transcripts = coding_transcripts(args.reference, genome, validate=True)
        query_transcripts = coding_transcripts(args.query, genome, validate=False)
        gene, counts = gene_level(args.reference, args.query, genome, workdir,
                                  bitscore=args.bitscore, evalue=args.evalue,
                                  coverage=args.coverage, threads=args.threads,
                                  ref_transcripts=ref_transcripts,
                                  query_transcripts=query_transcripts)
        reference_longest = workdir / "reference.longest.cds.gff3"
        query_cds = workdir / "query.cds.gff3"
        reference_count = filtered_gff(args.reference, genome, reference_longest,
                                       longest_only=True, validate=True, transcripts=ref_transcripts)
        query_count = filtered_gff(args.query, genome, query_cds,
                                   validate=False, transcripts=query_transcripts)
        exon_gff = scores(0, query_count, reference_count)
        if reference_count and query_count:
            exon_gff = run_gffcompare(reference_longest, query_cds, workdir, "gffcompare_exon")["exon"]
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)
    print(f"gene-level (protein): F1={gene[2]:.4%} precision={gene[0]:.4%} recall={gene[1]:.4%}")
    print(f"  genes TP={counts[0]:,} P(reference)={counts[1]:,} PP(predicted)={counts[2]:,}")
    print(f"gffcompare CDS exon-level: F1={exon_gff[2]:.2%} precision={exon_gff[0]:.1%} recall={exon_gff[1]:.1%}")
    print(f"  filtered transcripts: reference={reference_count:,} query={query_count:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
