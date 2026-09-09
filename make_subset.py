import gzip
import os

fasta_path = "Col-CEN_v1.2.fasta.gz"
gff_path = "Col-CEN_v1.2_genes.araport11.gff3"
if not os.path.exists(gff_path) and os.path.exists(gff_path + ".gz"):
    gff_path = gff_path + ".gz"

out_fasta = "subset_genome.fasta"
out_gff = "subset_annotation.gff3"

max_total_len = 500000 # 13000000 # 13Mb

print(f"Starting subset extraction from {fasta_path}...")

# 각 염색체별로 추출된 길이를 기록할 딕셔너리
kept_chroms = {}
total_extracted_bp = 0

# FASTA에서 서열을 염색체별로 읽으며 누적 100Mb까지 저장
current_chr = None
current_seq_parts = []

def write_chr_seq(name, seq_data, limit):
    truncated_seq = "".join(seq_data)[:limit]
    with open(out_fasta, "a" if os.path.exists(out_fasta) else "w") as out:
        out.write(f">{name}\n")
        for i in range(0, len(truncated_seq), 80):
            out.write(truncated_seq[i:i+80] + "\n")
    return len(truncated_seq)

# 기존 출력 파일 삭제
if os.path.exists(out_fasta):
    os.remove(out_fasta)

with gzip.open(fasta_path, "rt") as f:
    for line in f:
        if line.startswith(">"):
            if current_chr is not None:
                bp_needed = max_total_len - total_extracted_bp
                if bp_needed <= 0:
                    break
                kept_len = write_chr_seq(current_chr, current_seq_parts, bp_needed)
                kept_chroms[current_chr] = kept_len
                total_extracted_bp += kept_len
                print(f"Extracted {kept_len} bp from {current_chr} (Total: {total_extracted_bp} bp)")
                if total_extracted_bp >= max_total_len:
                    current_chr = None
                    break
            
            current_chr = line.strip().split()[0][1:]
            current_seq_parts = []
        else:
            if current_chr is not None:
                current_seq_parts.append(line.strip())

# 마지막 염색체 처리
if current_chr is not None and total_extracted_bp < max_total_len:
    bp_needed = max_total_len - total_extracted_bp
    kept_len = write_chr_seq(current_chr, current_seq_parts, bp_needed)
    kept_chroms[current_chr] = kept_len
    total_extracted_bp += kept_len
    print(f"Extracted {kept_len} bp from {current_chr} (Total: {total_extracted_bp} bp)")

print(f"FASTA extraction complete. Total extracted: {total_extracted_bp} bp across {len(kept_chroms)} chromosomes.")
print(f"Chromosomes kept: {kept_chroms}")

# GFF3에서 주석 추출
print(f"Extracting matching annotations from {gff_path}...")
annotation_count = 0
gff_opener = gzip.open if gff_path.endswith(".gz") else open
with gff_opener(gff_path, "rt") as f, open(out_gff, "w") as out:
    for line in f:
        if line.startswith("#"):
            out.write(line)
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 5:
            seqid = parts[0]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue
            
            if seqid in kept_chroms:
                limit = kept_chroms[seqid]
                if start >= 1 and end <= limit:
                    out.write(line)
                    annotation_count += 1

print(f"GFF3 extraction complete. Total annotation records written: {annotation_count}")
