#define _POSIX_C_SOURCE 200809L
#include "anno_data.h"
#include "anno_hmm.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
  int start, end, phase;
} CDS;

typedef struct {
  char *id, *gene, *alias;
  int sequence, strand, copy, count, capacity, selected;
  long total;
  CDS *blocks;
} Transcript;

typedef struct {
  char *parent, *gene;
  int sequence, strand, copy;
  CDS block;
} PendingCDS;

static char *copy_string(const char *text) {
  char *copy = anno_alloc(strlen(text) + 1, 1);
  strcpy(copy, text);
  return copy;
}

static void *grow(void *memory, size_t count, size_t size) {
  if (size && count > (size_t)-1 / size)
    anno_fail("Annotation capacity overflow");
  void *result = realloc(memory, count * size);
  if (!result)
    anno_fail("Out of memory loading annotation");
  return result;
}

Dataset *dataset_read(const char *path) {
  gzFile stream = gzopen(path, "rb");
  if (!stream)
    anno_fail("Cannot open FASTA: %s", path);
  kseq_t *record = kseq_init(stream);
  Dataset *dataset = anno_alloc(1, sizeof(*dataset));
  int result;
  while ((result = kseq_read(record)) >= 0) {
    if (!record->seq.l)
      continue;
    if (record->seq.l > INT_MAX - MAX_WAVE_SIZE)
      anno_fail("FASTA record exceeds supported length: %s", record->name.s);
    for (int sequence = 0; sequence < dataset->count; sequence++)
      if (!strcmp(dataset->sequences[sequence].name, record->name.s))
        anno_fail("Duplicate FASTA identifier: %s", record->name.s);
    dataset->sequences = grow(dataset->sequences, dataset->count + 1, sizeof(Sequence));
    Sequence *sequence = &dataset->sequences[dataset->count++];
    memset(sequence, 0, sizeof(*sequence));
    sequence->name = copy_string(record->name.s);
    sequence->bases = copy_string(record->seq.s);
    sequence->length = (int)record->seq.l;
  }
  int gz_error;
  gzerror(stream, &gz_error);
  if (result < -1 || (gz_error != Z_OK && gz_error != Z_STREAM_END))
    anno_fail("Malformed or unreadable FASTA: %s", path);
  kseq_destroy(record);
  gzclose(stream);
  if (!dataset->count)
    anno_fail("No nonempty sequences in FASTA: %s", path);
  return dataset;
}

void dataset_free(Dataset *dataset) {
  for (int index = 0; index < dataset->count; index++) {
    Sequence *sequence = &dataset->sequences[index];
    free(sequence->name);
    free(sequence->bases);
    for (int strand = 0; strand < 2; strand++) {
      free(sequence->coding[strand]);
      free(sequence->states[strand]);
    }
  }
  free(dataset->sequences);
  free(dataset);
}

static int hex_digit(char character) {
  if (character >= '0' && character <= '9') return character - '0';
  if (character >= 'a' && character <= 'f') return character - 'a' + 10;
  if (character >= 'A' && character <= 'F') return character - 'A' + 10;
  return -1;
}

static char *attribute(const char *attributes, const char *key) {
  const char *cursor = attributes;
  while (*cursor) {
    while (*cursor == ';' || isspace((unsigned char)*cursor)) cursor++;
    const char *name = cursor;
    while (*cursor && *cursor != '=' && *cursor != ';' &&
           !isspace((unsigned char)*cursor)) cursor++;
    size_t name_length = (size_t)(cursor - name);
    while (isspace((unsigned char)*cursor)) cursor++;
    int encoded = *cursor == '=';
    if (encoded) cursor++;
    while (isspace((unsigned char)*cursor)) cursor++;
    int quoted = *cursor == '"';
    if (quoted) cursor++;
    const char *value = cursor;
    while (*cursor && (quoted ? *cursor != '"' : *cursor != ';')) cursor++;
    const char *end = cursor;
    if (quoted && *cursor == '"') cursor++;
    while (end > value && isspace((unsigned char)end[-1])) end--;
    if (strlen(key) == name_length && !strncmp(name, key, name_length)) {
      char *result = anno_alloc((size_t)(end - value) + 1, 1);
      size_t used = 0;
      while (value < end) {
        if (encoded && *value == '%' && end - value >= 3 &&
            hex_digit(value[1]) >= 0 && hex_digit(value[2]) >= 0) {
          int decoded = 16 * hex_digit(value[1]) + hex_digit(value[2]);
          if (!decoded) anno_fail("NUL in GFF3 attribute");
          result[used++] = (char)decoded;
          value += 3;
        } else {
          result[used++] = *value++;
        }
      }
      return result;
    }
    while (*cursor && *cursor != ';') cursor++;
  }
  return NULL;
}

static int coordinate(const char *text, int length, long line_number) {
  char *end;
  errno = 0;
  long value = strtol(text, &end, 10);
  if (errno || *end || value < 1 || value > length)
    anno_fail("Out-of-range GFF coordinate at line %ld", line_number);
  return (int)value;
}

static int transcript_find(Transcript *transcripts, int count, int sequence,
                            int strand, int copy, const char *id) {
  for (int index = 0; index < count; index++)
    if (transcripts[index].sequence == sequence && transcripts[index].strand == strand &&
        transcripts[index].copy == copy &&
        (!strcmp(transcripts[index].id, id) ||
         (transcripts[index].alias && !strcmp(transcripts[index].alias, id))))
      return index;
  return -1;
}

static int block_compare(const void *first, const void *second) {
  const CDS *left = first, *right = second;
  if (left->start != right->start) return left->start < right->start ? -1 : 1;
  return (left->end > right->end) - (left->end < right->end);
}

void annotation_read(Dataset *dataset, FILE *stream) {
  Transcript *transcripts = NULL;
  PendingCDS *pending = NULL;
  int transcript_count = 0, pending_count = 0;
  char *line = NULL;
  size_t capacity = 0;
  long line_number = 0;
  while (getline(&line, &capacity, stream) >= 0) {
    line_number++;
    if (!strncmp(line, "##FASTA", 7)) break;
    if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
    line[strcspn(line, "\r\n")] = 0;
    char *fields[9], *cursor = line;
    for (int field = 0; field < 9; field++) {
      fields[field] = cursor;
      if (field < 8) {
        char *tab = strchr(cursor, '\t');
        if (!tab) anno_fail("Expected 9 tab-separated GFF fields at line %ld", line_number);
        *tab = 0;
        cursor = tab + 1;
      }
    }
    int sequence = -1;
    for (int index = 0; index < dataset->count; index++)
      if (!strcmp(dataset->sequences[index].name, fields[0])) sequence = index;
    if (sequence < 0) continue;
    int is_cds = !strcmp(fields[2], "CDS");
    int is_transcript = !strcmp(fields[2], "mRNA") || !strcmp(fields[2], "transcript");
    if (!is_cds && !is_transcript) continue;
    if (strcmp(fields[6], "+") && strcmp(fields[6], "-"))
      anno_fail("CDS/transcript requires strand at line %ld", line_number);
    int strand = fields[6][0] == '-';
    int start = coordinate(fields[3], dataset->sequences[sequence].length, line_number) - 1;
    int end = coordinate(fields[4], dataset->sequences[sequence].length, line_number);
    if (start >= end) anno_fail("Reversed GFF interval at line %ld", line_number);
    char *id = attribute(fields[8], "ID");
    char *parent = attribute(fields[8], "Parent");
    char *gene = attribute(fields[8], "gene_id");
    char *transcript_id = attribute(fields[8], "transcript_id");
    char *copy_text = attribute(fields[8], "extra_copy_number");
    int copy = 0;
    if (copy_text) {
      char *copy_end;
      errno = 0;
      long parsed = strtol(copy_text, &copy_end, 10);
      if (errno || !*copy_text || *copy_end || parsed < 0 || parsed > INT_MAX)
        anno_fail("Invalid Liftoff copy number at line %ld", line_number);
      copy = (int)parsed;
    }
    if (is_transcript) {
      const char *name = transcript_id ? transcript_id : id;
      const char *gene_name = gene ? gene : parent;
      if (!name || !*name || !gene_name || !*gene_name)
        anno_fail("Transcript needs ID and gene Parent (or GTF IDs), line %ld", line_number);
      if (transcript_find(transcripts, transcript_count, sequence, strand, copy, name) >= 0)
        anno_fail("Duplicate transcript ID: %s", name);
      transcripts = grow(transcripts, transcript_count + 1, sizeof(Transcript));
      Transcript *transcript = &transcripts[transcript_count++];
      memset(transcript, 0, sizeof(*transcript));
      transcript->id = copy_string(name);
      transcript->alias = id ? copy_string(id) : NULL;
      transcript->gene = copy_string(gene_name);
      transcript->sequence = sequence;
      transcript->strand = strand;
      transcript->copy = copy;
    } else {
      const char *names = transcript_id ? transcript_id : parent;
      if (!names || !*names)
        anno_fail("CDS needs transcript Parent or transcript_id, line %ld", line_number);
      int phase = -1;
      if (strcmp(fields[7], ".")) {
        if (strlen(fields[7]) != 1 || fields[7][0] < '0' || fields[7][0] > '2')
          anno_fail("Invalid CDS phase at line %ld", line_number);
        phase = fields[7][0] - '0';
      }
      char *parents = copy_string(names), *save = NULL;
      for (char *name = strtok_r(parents, ",", &save); name; name = strtok_r(NULL, ",", &save)) {
        pending = grow(pending, pending_count + 1, sizeof(PendingCDS));
        PendingCDS *record = &pending[pending_count++];
        record->parent = copy_string(name);
        record->gene = gene ? copy_string(gene) : NULL;
        record->sequence = sequence;
        record->strand = strand;
        record->copy = copy;
        record->block = (CDS){start, end, phase};
      }
      free(parents);
    }
    free(id); free(parent); free(gene); free(transcript_id); free(copy_text);
  }
  free(line);
  if (ferror(stream)) anno_fail("Error reading annotation");
  for (int index = 0; index < pending_count; index++) {
    PendingCDS *record = &pending[index];
    int found = transcript_find(transcripts, transcript_count, record->sequence,
                                record->strand, record->copy, record->parent);
    if (found < 0) {
      if (!record->gene)
        anno_fail("Cannot resolve transcript %s to a gene; supply mRNA Parent or gene_id", record->parent);
      transcripts = grow(transcripts, transcript_count + 1, sizeof(Transcript));
      found = transcript_count++;
      memset(&transcripts[found], 0, sizeof(Transcript));
      transcripts[found].id = copy_string(record->parent);
      transcripts[found].gene = copy_string(record->gene);
      transcripts[found].sequence = record->sequence;
      transcripts[found].strand = record->strand;
      transcripts[found].copy = record->copy;
    }
    Transcript *transcript = &transcripts[found];
    if (record->gene && strcmp(record->gene, transcript->gene))
      anno_fail("Conflicting gene for transcript %s", transcript->id);
    if (transcript->count == transcript->capacity) {
      transcript->capacity = transcript->capacity ? transcript->capacity * 2 : 8;
      transcript->blocks = grow(transcript->blocks, transcript->capacity, sizeof(CDS));
    }
    transcript->blocks[transcript->count++] = record->block;
    free(record->parent); free(record->gene);
  }
  free(pending);
  for (int index = 0; index < transcript_count; index++) {
    Transcript *transcript = &transcripts[index];
    if (!transcript->count) continue;
    qsort(transcript->blocks, transcript->count, sizeof(CDS), block_compare);
    int unique = 0;
    for (int block = 0; block < transcript->count; block++) {
      CDS current = transcript->blocks[block];
      if (unique && current.start == transcript->blocks[unique - 1].start &&
          current.end == transcript->blocks[unique - 1].end) continue;
      if (unique && current.start < transcript->blocks[unique - 1].end)
        anno_fail("Overlapping CDS within transcript %s", transcript->id);
      transcript->blocks[unique++] = current;
      transcript->total += current.end - current.start;
    }
    transcript->count = unique;
  }
  for (int index = 0; index < transcript_count; index++) {
    Transcript *candidate = &transcripts[index];
    if (!candidate->count) continue;
    candidate->selected = 1;
    for (int other = 0; other < transcript_count; other++) {
      Transcript *competitor = &transcripts[other];
      if (competitor->sequence == candidate->sequence && competitor->strand == candidate->strand &&
          competitor->copy == candidate->copy &&
          !strcmp(competitor->gene, candidate->gene) &&
          (competitor->total > candidate->total ||
           (competitor->total == candidate->total && strcmp(competitor->id, candidate->id) < 0)))
        candidate->selected = 0;
    }
  }
  for (int sequence = 0; sequence < dataset->count; sequence++)
    for (int strand = 0; strand < 2; strand++) {
      dataset->sequences[sequence].coding[strand] = anno_alloc(dataset->sequences[sequence].length, 1);
      dataset->sequences[sequence].states[strand] = anno_alloc(dataset->sequences[sequence].length, 1);
    }
  for (int index = 0; index < transcript_count; index++) {
    Transcript *transcript = &transcripts[index];
    if (!transcript->selected) continue;
    dataset->selected_transcripts++;
    Sequence *sequence = &dataset->sequences[transcript->sequence];
    unsigned char *coding = sequence->coding[transcript->strand];
    unsigned char *states = sequence->states[transcript->strand];
    int phase = 0, previous_end = -1;
    for (int block = 0; block < transcript->count; block++) {
      int block_index = transcript->strand ? transcript->count - block - 1 : block;
      CDS cds = transcript->blocks[block_index];
      int start = transcript->strand ? sequence->length - cds.end : cds.start;
      int end = transcript->strand ? sequence->length - cds.start : cds.end;
      if (!block && cds.phase >= 0) phase = (3 - cds.phase) % 3;
      if (cds.phase >= 0 && cds.phase != (3 - phase) % 3)
        anno_fail("Inconsistent CDS phase in transcript %s", transcript->id);
      if (previous_end >= 0 && start > previous_end) {
        int gap = start - previous_end;
        for (int position = previous_end; position < start; position++) {
          int offset = position - previous_end;
          int state = gap < 5 ? HMM_UNKNOWN : HMM_DONOR(phase) +
                      (offset < 2 ? offset : (position >= start - 2 ? 3 + position - (start - 2) : 2));
          states[position] = states[position] ? HMM_UNKNOWN : state;
        }
      }
      for (int position = start; position < end; position++) {
        coding[position] = 1;
        states[position] = states[position] ? HMM_UNKNOWN : HMM_C0 + phase;
        phase = (phase + 1) % 3;
      }
      previous_end = end;
    }
  }
  for (int index = 0; index < transcript_count; index++) {
    free(transcripts[index].id); free(transcripts[index].gene);
    free(transcripts[index].alias); free(transcripts[index].blocks);
  }
  free(transcripts);
  if (!dataset->selected_transcripts)
    anno_fail("No usable CDS transcripts match FASTA identifiers");
  fprintf(stderr, "Selected %d transcripts by summed CDS length (one per gene)\n",
          dataset->selected_transcripts);
}