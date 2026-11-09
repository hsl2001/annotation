/**
 * @file test.c
 * @brief Small, hand-checkable unit tests for every public surf function.
 */

#include "../src/include/crf.h"
#include "../src/include/lbfgs.h"
#include "../src/include/surf.h"
#include "../src/include/utils.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int failures;
static int checks;

#define CHECK(condition)                                                       \
  do {                                                                         \
    checks++;                                                                  \
    if (!(condition)) {                                                        \
      failures++;                                                              \
      fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__,         \
              #condition);                                                     \
    }                                                                          \
  } while (0)

static FILE *make_temp(char path[], const char *contents) {
  strcpy(path, "/tmp/surf-test-XXXXXX");
  int descriptor = mkstemp(path);
  CHECK(descriptor >= 0);
  if (descriptor < 0)
    return NULL;
  FILE *file = fdopen(descriptor, "w+");
  CHECK(file != NULL);
  if (!file) {
    close(descriptor);
    unlink(path);
    return NULL;
  }
  fputs(contents, file);
  fflush(file);
  rewind(file);
  return file;
}

static int close_enough(double actual, double expected) {
  return fabs(actual - expected) < 1e-10;
}

static double quadratic(const double *x, double *gradient, void *context) {
  (void)context;
  gradient[0] = 2.0 * (x[0] - 3.0);
  return (x[0] - 3.0) * (x[0] - 3.0);
}

/* Each character is one base: B is background; S/E/I are plus, s/e/i minus. */
static void label_path(const char *symbols, uint8_t *path) {
  const char *alphabet = "BSEIsei";
  const uint8_t classes[] = {INTERGENIC, PLUS_START_SINGLE_EXON, PLUS_END_EXON,
                            PLUS_INTRON, MINUS_START_SINGLE_EXON,
                            MINUS_END_EXON, MINUS_INTRON};
  for (size_t p = 0; symbols[p]; ++p) {
    const char *state = strchr(alphabet, symbols[p]);
    CHECK(state != NULL);
    path[p] = state ? classes[state - alphabet] : INTERGENIC;
  }
}

/* Noncoding fixtures emit only gene, mRNA, and genomic-order exon records. */
static void check_gff_transcript(FILE *output, const Interval *exons, int count,
                                char strand, unsigned long gene) {
  for (int row = 0; row < count + 2; ++row) {
    char line[256] = {0}, type[32] = {0}, attributes[128] = {0};
    char actual_strand = 0, expected[128];
    long long start = 0, end = 0;
    CHECK(fgets(line, sizeof(line), output) != NULL);
    CHECK(sscanf(line, "%*s %*s %31s %lld %lld %*s %c %*s %127s", type,
                 &start, &end, &actual_strand, attributes) == 5);
    CHECK(!strcmp(type, row == 0 ? "gene" : row == 1 ? "mRNA" : "exon"));
    CHECK(start == (row < 2 ? exons[0].start : exons[row - 2].start) + 1 &&
          end == (row < 2 ? exons[count - 1].end : exons[row - 2].end) &&
          actual_strand == strand);
    if (row == 0)
      snprintf(expected, sizeof(expected), "ID=gene%lu", gene);
    else if (row == 1)
      snprintf(expected, sizeof(expected), "ID=gene%lu.t1;Parent=gene%lu", gene,
               gene);
    else
      snprintf(expected, sizeof(expected),
               "ID=gene%lu.t1.exon%d;Parent=gene%lu.t1", gene, row - 1, gene);
    CHECK(!strcmp(attributes, expected));
  }
}

static void test_utils(void) {
  CHECK(base_index('A') == 0);
  CHECK(base_index('C') == 1);
  CHECK(base_index('G') == 2);
  CHECK(base_index('T') == 3);
  CHECK(base_index('N') == -1);
  CHECK(complement('A') == 'T');
  CHECK(complement('C') == 'G');
  CHECK(complement('G') == 'C');
  CHECK(complement('T') == 'A');
  CHECK(complement('N') == 'N');
  CHECK(smul(6, 7) == 42);
  int *zeroes = scalloc(3, sizeof(*zeroes));
  CHECK(zeroes[0] == 0 && zeroes[1] == 0 && zeroes[2] == 0);
  zeroes = srealloc(zeroes, 4, sizeof(*zeroes));
  zeroes[3] = 9;
  CHECK(zeroes[3] == 9);
  free(zeroes);
  char *copy = sstrdup("surf");
  CHECK(!strcmp(copy, "surf"));
  free(copy);
  double before = whattime();
  double after = whattime();
  CHECK(isfinite(before) && after >= before);
  cry("utility test message");
}

static void test_fasta_and_gff_input(void) {
  char fasta_path[32], gff_path[32];
  FILE *fasta = make_temp(fasta_path, ">chr1 description\natg aaa taa\n");
  FILE *gff = make_temp(gff_path,
                        "##gff-version 3\n"
                        "chr1\tanno\texon\t1\t9\t.\t+\t.\tID=e1;Parent=tx1\n");
  CHECK(fasta && gff);
  if (!fasta || !gff)
    return;
  fclose(fasta);
  fclose(gff);

  int count = 0;
  Sequence *sequences = NULL;
  read_fasta(fasta_path, &count, &sequences);
  CHECK(count == 1);
  CHECK(!strcmp(sequences[0].name, "chr1"));
  CHECK(!memcmp(sequences[0].bases, "ATGAAATAA", 9));
  CHECK(sequences[0].length == 9);
  read_gff3(gff_path, sequences, count);
  CHECK(sequences[0].labels != NULL);
  for (int i = 0; i < 9; i++)
    CHECK(sequences[0].labels[i] == PLUS_START_SINGLE_EXON);
  free_sequences(sequences, count);
  unlink(fasta_path);
  unlink(gff_path);
}

static void test_cwt_window_independence(void) {
  enum { LENGTH = 2400 };
  char bases[LENGTH];
  float whole[LENGTH * CWT_CHANNELS], split[LENGTH * CWT_CHANNELS];
  for (int i = 0; i < LENGTH; i++)
    bases[i] = "ACGTN"[(i * 13 + i / 7) % 5];
  cwt_power(bases, LENGTH, 0, LENGTH, whole);
  cwt_power(bases, LENGTH, 0, 331, split);
  cwt_power(bases, LENGTH, 331, 1400, split + 331 * CWT_CHANNELS);
  cwt_power(bases, LENGTH, 1400, LENGTH, split + 1400 * CWT_CHANNELS);
  for (int i = 0; i < LENGTH * CWT_CHANNELS; i++)
    CHECK(fabsf(whole[i] - split[i]) < 1e-5f);
}

static void test_intervals_and_gff_output(void) {
  IntervalList list = {0};
  push_interval(&list, 4, 8);
  push_interval(&list, 1, 3);
  CHECK(list.count == 2 && list.capacity >= 2);
  qsort(list.items, (size_t)list.count, sizeof(*list.items), compare_interval);
  CHECK(list.items[0].start == 1 && list.items[1].start == 4);
  CHECK(compare_interval(&list.items[0], &list.items[1]) < 0);
  free(list.items);

  char bases[] = "CCCATGAAATAAGGG";
  uint8_t path[15];
  for (int i = 0; i < 15; i++)
    path[i] = PLUS_START_SINGLE_EXON;
  Sequence sequence = {(char *)"chr1", bases, 15, NULL};
  FILE *output = tmpfile();
  unsigned long gene = 0;
  CHECK(output != NULL);
  CHECK(write_gff(output, &sequence, path, &gene) == 1);
  CHECK(gene == 1);
  rewind(output);
  char text[2048] = {0};
  CHECK(fread(text, 1, sizeof(text) - 1, output) > 0);
  CHECK(strstr(text, "gene1") != NULL);
  CHECK(strstr(text, "\tCDS\t") != NULL);
  fclose(output);

  sequence.bases = "CCCCCCCCCCCCCCC";
  for (int strand = 0; strand < 2; strand++) {
    for (int i = 0; i < 15; i++)
      path[i] = strand ? MINUS_START_SINGLE_EXON : PLUS_START_SINGLE_EXON;
    output = tmpfile();
    CHECK(output != NULL);
    if (!output)
      continue;
    gene = 0;
    CHECK(write_gff(output, &sequence, path, &gene) == 1 && gene == 1);
    rewind(output);
    memset(text, 0, sizeof(text));
    CHECK(fread(text, 1, sizeof(text) - 1, output) > 0);
    CHECK(strstr(text, "\tgene\t") && strstr(text, "\tmRNA\t") &&
          strstr(text, "\texon\t") && !strstr(text, "\tCDS\t") &&
          !strstr(text, "UTR\t") && strstr(text, strand ? "\t-\t" : "\t+\t"));
    fclose(output);
  }
}

static void test_multiexon_gff(void) {
  const char *bases[] = {"NCCATGANNAAACCNNGGGTAACTN",
                         "NAGTTACCCNNGGTTTNNTCATGGN"};
  const int starts[][3] = {{1, 9, 16}, {1, 11, 18}};
  const int ends[][3] = {{7, 14, 24}, {9, 16, 24}};
  const char *expected[] = {
      "chr1\tanno\tgene\t2\t24\t.\t+\t.\tID=gene1\n"
      "chr1\tanno\tmRNA\t2\t24\t.\t+\t.\tID=gene1.t1;Parent=gene1\n"
      "chr1\tanno\texon\t2\t7\t.\t+\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
      "chr1\tanno\texon\t10\t14\t.\t+\t.\tID=gene1.t1.exon2;Parent=gene1.t1\n"
      "chr1\tanno\texon\t17\t24\t.\t+\t.\tID=gene1.t1.exon3;Parent=gene1.t1\n"
      "chr1\tanno\tfive_prime_UTR\t2\t3\t.\t+\t.\tID=gene1.t1.utr5p1;Parent="
      "gene1.t1\n"
      "chr1\tanno\tCDS\t4\t7\t.\t+\t0\tID=gene1.t1.cds1;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t10\t14\t.\t+\t2\tID=gene1.t1.cds2;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t17\t22\t.\t+\t0\tID=gene1.t1.cds3;Parent=gene1.t1\n"
      "chr1\tanno\tthree_prime_UTR\t23\t24\t.\t+\t.\tID=gene1.t1.utr3p1;Parent="
      "gene1.t1\n",
      "chr1\tanno\tgene\t2\t24\t.\t-\t.\tID=gene1\n"
      "chr1\tanno\tmRNA\t2\t24\t.\t-\t.\tID=gene1.t1;Parent=gene1\n"
      "chr1\tanno\texon\t2\t9\t.\t-\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
      "chr1\tanno\texon\t12\t16\t.\t-\t.\tID=gene1.t1.exon2;Parent=gene1.t1\n"
      "chr1\tanno\texon\t19\t24\t.\t-\t.\tID=gene1.t1.exon3;Parent=gene1.t1\n"
      "chr1\tanno\tfive_prime_UTR\t23\t24\t.\t-\t.\tID=gene1.t1.utr5p1;Parent="
      "gene1.t1\n"
      "chr1\tanno\tCDS\t4\t9\t.\t-\t0\tID=gene1.t1.cds1;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t12\t16\t.\t-\t2\tID=gene1.t1.cds2;Parent=gene1.t1\n"
      "chr1\tanno\tCDS\t19\t22\t.\t-\t0\tID=gene1.t1.cds3;Parent=gene1.t1\n"
      "chr1\tanno\tthree_prime_UTR\t2\t3\t.\t-\t.\tID=gene1.t1.utr3p1;Parent="
      "gene1.t1\n"};
  const uint8_t classes[][3] = {
      {PLUS_START_SINGLE_EXON, INTERGENIC, PLUS_END_EXON},
      {MINUS_END_EXON, INTERGENIC, MINUS_START_SINGLE_EXON}};
  for (int strand = 0; strand < 2; strand++) {
    uint8_t path[25] = {0};
    for (int e = 0; e < 3; e++) {
      for (int p = starts[strand][e]; p < ends[strand][e]; p++)
        path[p] = classes[strand][e];
      if (e < 2)
        for (int p = ends[strand][e]; p < starts[strand][e + 1]; p++)
          path[p] = strand ? MINUS_INTRON : PLUS_INTRON;
    }
    Sequence sequence = {(char *)"chr1", (char *)bases[strand], 25, NULL};
    FILE *output = tmpfile();
    CHECK(output != NULL);
    if (!output)
      continue;
    unsigned long gene = 0;
    CHECK(write_gff(output, &sequence, path, &gene) == 1 && gene == 1);
    rewind(output);
    char text[2048] = {0};
    CHECK(fread(text, 1, sizeof(text) - 1, output) > 0);
    CHECK(!strcmp(text, expected[strand]));
    fclose(output);

    char file_path[32];
    FILE *file = make_temp(file_path, "");
    CHECK(file != NULL);
    if (file) {
      fclose(file);
      file = fopen(file_path, "r");
      CHECK(file != NULL);
      if (file) {
        gene = 0;
        CHECK(write_gff(file, &sequence, path, &gene) == -1 && gene == 0);
        fclose(file);
      }
      unlink(file_path);
    }
  }
}

static void test_multiexon_gff_input(void) {
  char path[32];
  FILE *gff =
      make_temp(path, "plus\tanno\texon\t2\t7\t.\t+\t.\tParent=tx1\n"
                      "plus\tanno\texon\t10\t14\t.\t+\t.\tParent=tx1\n"
                      "plus\tanno\texon\t17\t24\t.\t+\t.\tParent=tx1\n"
                      "minus\tanno\texon\t2\t9\t.\t-\t.\tParent=tx2\n"
                      "minus\tanno\texon\t12\t16\t.\t-\t.\tParent=tx2\n"
                      "minus\tanno\texon\t19\t24\t.\t-\t.\tParent=tx2\n"
                      "plus\tanno\texon\t11\t13\t.\t-\t.\tParent=nested1\n"
                      "plus\tanno\texon\t11\t11\t.\t+\t.\tParent=nested2\n"
                      "plus\tanno\texon\t13\t13\t.\t+\t.\tParent=nested2\n"
                      "minus\tanno\texon\t13\t15\t.\t+\t.\tParent=nested3\n"
                      "minus\tanno\texon\t13\t13\t.\t-\t.\tParent=nested4\n"
                      "minus\tanno\texon\t15\t15\t.\t-\t.\tParent=nested4\n");
  CHECK(gff != NULL);
  if (!gff)
    return;
  fclose(gff);
  Sequence sequences[] = {
      {(char *)"plus", (char *)"NCCATGANNAAACCNNGGGTAACTN", 25, NULL},
      {(char *)"minus", (char *)"NAGTTACCCNNGGTTTNNTCATGGN", 25, NULL}};
  read_gff3(path, sequences, 2);
  CHECK(sequences[0].labels[1] == PLUS_START_SINGLE_EXON);
  CHECK(sequences[0].labels[7] == PLUS_INTRON);
  /* Nested transcripts must not overwrite the unmarked internal exons. */
  for (int p = 9; p < 14; ++p)
    CHECK(sequences[0].labels[p] == INTERGENIC);
  CHECK(sequences[0].labels[16] == PLUS_END_EXON);
  CHECK(sequences[1].labels[1] == MINUS_END_EXON);
  CHECK(sequences[1].labels[9] == MINUS_INTRON);
  for (int p = 11; p < 16; ++p)
    CHECK(sequences[1].labels[p] == INTERGENIC);
  CHECK(sequences[1].labels[18] == MINUS_START_SINGLE_EXON);
  free(sequences[0].labels);
  free(sequences[1].labels);
  unlink(path);
}

static void test_gff_roundtrip(void) {
  const struct {
    const char *symbols;
    int exons;
    char strand;
  } cases[] = {{"BSSB", 1, '+'},
               {"BssB", 1, '-'},
               {"BSSIIEEB", 2, '+'},
               {"BeeiissB", 2, '-'},
               {"BSSIIBBIIBBIIEEB", 4, '+'},
               {"BeeiiBBiiBBiissB", 4, '-'}};
  const Interval exons[] = {{1, 3}, {5, 7}, {9, 11}, {13, 15}};
  char bases[] = "CCCCCCCCCCCCCCCC";
  for (size_t c = 0; c < sizeof(cases) / sizeof(*cases); ++c) {
    char file_path[32];
    FILE *input = make_temp(file_path, "");
    CHECK(input != NULL);
    if (!input)
      continue;
    /* Reverse input order so reconstruction cannot rely on GFF row order. */
    for (int e = cases[c].exons - 1; e >= 0; --e)
      fprintf(input, "chr1\tanno\texon\t%lld\t%lld\t.\t%c\t.\tParent=tx\n",
              (long long)exons[e].start + 1, (long long)exons[e].end,
              cases[c].strand);
    CHECK(fclose(input) == 0);
    Sequence sequence = {(char *)"chr1", bases,
                         (int64_t)strlen(cases[c].symbols), NULL};
    read_gff3(file_path, &sequence, 1);
    unlink(file_path);
    uint8_t expected[16], decoded[16];
    label_path(cases[c].symbols, expected);
    CHECK(!memcmp(sequence.labels, expected, (size_t)sequence.length));
    Decoder *decoder = decoder_create(sequence.length);
    CHECK(decoder != NULL);
    if (!decoder) {
      free(sequence.labels);
      continue;
    }
    memset(decoded, 255, sizeof(decoded));
    for (int64_t p = 0; p < sequence.length; ++p) {
      double emissions[CLASSES];
      for (int state = 0; state < CLASSES; ++state)
        emissions[state] = state == sequence.labels[p] ? 100.0 : 0.0;
      decoder_step(decoder, emissions, NULL);
    }
    CHECK(decoder_finish(decoder, decoded) == 0);
    CHECK(!memcmp(decoded, expected, (size_t)sequence.length));
    decoder_destroy(decoder);
    FILE *output = tmpfile();
    CHECK(output != NULL);
    if (output) {
      unsigned long gene = 0;
      CHECK(write_gff(output, &sequence, decoded, &gene) == 1 && gene == 1);
      rewind(output);
      check_gff_transcript(output, exons, cases[c].exons, cases[c].strand, 1);
      CHECK(fgetc(output) == EOF && !ferror(output));
      fclose(output);
    }
    free(sequence.labels);
  }
}

static void test_gff_neighboring_genes(void) {
  const char *symbols[] = {"BSSIIEEBSSIIEEBssB", "BeeiissBeeiissBSSB"};
  const Interval first[] = {{1, 3}, {5, 7}};
  const Interval second[] = {{8, 10}, {12, 14}};
  const Interval single[] = {{15, 17}};
  char bases[32];
  memset(bases, 'C', sizeof(bases));
  for (int strand = 0; strand < 2; ++strand) {
    uint8_t path[32];
    label_path(symbols[strand], path);
    Sequence sequence = {(char *)"chr1", bases,
                         (int64_t)strlen(symbols[strand]), NULL};
    FILE *output = tmpfile();
    CHECK(output != NULL);
    if (!output)
      continue;
    unsigned long gene = 10;
    CHECK(write_gff(output, &sequence, path, &gene) == 3 && gene == 13);
    rewind(output);
    check_gff_transcript(output, first, 2, strand ? '-' : '+', 11);
    check_gff_transcript(output, second, 2, strand ? '-' : '+', 12);
    check_gff_transcript(output, single, 1, strand ? '+' : '-', 13);
    CHECK(fgetc(output) == EOF && !ferror(output));
    fclose(output);
  }
}

static void test_gff_invalid_paths(void) {
  const char *symbols[] = {
      "BBBB", "BIIB", "BiiB", "BEEB", "BeeB", /* Orphan runs. */
      "BIIEEB", "BiissB",                    /* Missing opening exon. */
      "BSSIIB", "BeeiiB", "BssiiB",          /* Incomplete or reversed. */
      "BSSII", "Beeii",                      /* Truncated terminal intron. */
      "BSSIIBBIIB", "BeeiiBBiiB",            /* Unclosed internal exons. */
      "BSSiiEEB", "BeeIIssB",                /* Wrong-strand intron. */
      "BSSIIBBiiEEB", "BeeiiBBIIssB",        /* Mixed internal chain. */
      "BSSIIssB", "BeeiiEEB",                /* Wrong-strand closing exon. */
      "BSSIIEEiiB", "BeeiissIIB"};           /* Missing terminal separator. */
  char bases[32];
  memset(bases, 'C', sizeof(bases));
  for (size_t c = 0; c < sizeof(symbols) / sizeof(*symbols); ++c) {
    uint8_t path[32];
    label_path(symbols[c], path);
    Sequence sequence = {(char *)"chr1", bases, (int64_t)strlen(symbols[c]),
                         NULL};
    FILE *output = tmpfile();
    CHECK(output != NULL);
    if (!output)
      continue;
    unsigned long gene = 9;
    CHECK(write_gff(output, &sequence, path, &gene) == 0 && gene == 9);
    rewind(output);
    CHECK(fgetc(output) == EOF && !ferror(output));
    fclose(output);
  }
}

static void test_gff_label_conflicts(void) {
  char path[32], bases[65] = {0};
  memset(bases, 'C', 64);
  memcpy(bases + 9, "ATGTAA", 6);
  Sequence sequence = {(char *)"chr", bases, 64, NULL};
  /* Reject malformed candidates and preserve the longest noncoding transcript.
   */
  FILE *gff =
      make_temp(path, "chr\tanno\texon\t2\t4\t.\t+\t.\tParent=valid\n"
                      "chr\tanno\texon\t33\t38\t.\t+\t.\tParent=valid\n"
                      "chr\tanno\texon\t2\t25\t.\t+\t.\tParent=overlap\n"
                      "chr\tanno\texon\t20\t45\t.\t+\t.\tParent=overlap\n"
                      "chr\tanno\texon\t2\t25\t.\t+\t.\tParent=adjacent\n"
                      "chr\tanno\texon\t26\t45\t.\t+\t.\tParent=adjacent\n"
                      "chr\tanno\texon\t10\t15\t.\t+\t.\tParent=orf\n"
                      "chr\tanno\texon\t1\t1\t.\t-\t.\tParent=enclosing\n"
                      "chr\tanno\texon\t40\t42\t.\t-\t.\tParent=enclosing\n"
                      "chr\tanno\texon\t39\t41\t.\t-\t.\tParent=touching\n");
  CHECK(gff != NULL);
  if (!gff)
    return;
  fclose(gff);
  read_gff3(path, &sequence, 1);
  uint8_t expected[64] = {0};
  memset(expected + 1, PLUS_START_SINGLE_EXON, 3);
  memset(expected + 4, PLUS_INTRON, 28);
  memset(expected + 32, PLUS_END_EXON, 6);
  CHECK(!memcmp(sequence.labels, expected, sizeof(expected)));
  uint8_t *labels = sequence.labels;
  memset(labels, 255, sizeof(expected));
  read_gff3(path, &sequence, 1);
  CHECK(sequence.labels == labels &&
        !memcmp(labels, expected, sizeof(expected)));
  free(sequence.labels);
  unlink(path);
}

static void test_decoder_and_grammar(void) {
  CHECK(CLASSES == 7 && GRAMMAR_MAX_BACKPOINTERS == 6);
  CHECK(sizeof(gene_grammar.backpointers) / sizeof(gene_grammar.backpointers[0]) ==
        CLASSES);
  CHECK(sizeof(gene_grammar.from) / sizeof(gene_grammar.from[0]) == CLASSES);
  CHECK(sizeof(gene_grammar.from[0]) / sizeof(gene_grammar.from[0][0]) ==
        GRAMMAR_MAX_BACKPOINTERS);
  const int backpointers[CLASSES] = {6, 2, 2, 3, 3, 2, 3};
  for (int state = 0; state < CLASSES; ++state) {
    int count = gene_grammar.backpointers[state];
    CHECK(count == backpointers[state]);
    CHECK(count > 0 && count <= GRAMMAR_MAX_BACKPOINTERS);
    if (count <= 0 || count > GRAMMAR_MAX_BACKPOINTERS)
      continue;
    for (int edge = 0; edge < count; ++edge) {
      int from = gene_grammar.from[state][edge];
      CHECK(from >= 0 && from < CLASSES);
      for (int previous = 0; previous < edge; ++previous)
        CHECK(from != gene_grammar.from[state][previous]);
    }
  }

  Decoder *decoder = decoder_create(1);
  CHECK(decoder != NULL);
  double emissions[CLASSES] = {0};
  uint8_t path[1] = {255};
  decoder_step(decoder, emissions, NULL);
  CHECK(decoder_finish(decoder, path) == 0 && path[0] == INTERGENIC);
  decoder_destroy(decoder);
  CHECK(decoder_create(0) == NULL);

  decoder = decoder_create(8);
  CHECK(decoder != NULL);
  uint8_t gene_path[8] = {255};
  for (int position = 0; position < 8; position++) {
    for (int class = 0; class < CLASSES; class++)
      emissions[class] =
          class == (position == 0 || position == 7 ? INTERGENIC
                                                   : PLUS_START_SINGLE_EXON)
              ? 100.0
              : 0.0;
    decoder_step(decoder, emissions, NULL);
  }
  CHECK(decoder_finish(decoder, gene_path) == 0);
  CHECK(gene_path[0] == INTERGENIC && gene_path[7] == INTERGENIC);
  for (int position = 1; position < 7; position++)
    CHECK(gene_path[position] == PLUS_START_SINGLE_EXON);
  decoder_destroy(decoder);
}

static void test_crf_and_lbfgs(void) {
  Model model = {0};
  float features[INPUT_CHANNELS] = {0};
  double scores[CLASSES];
  model.emission_weights[2 * (INPUT_CHANNELS + 1) + 0] = 2.0;
  model.emission_weights[2 * (INPUT_CHANNELS + 1) + INPUT_CHANNELS] = 1.0;
  features[0] = 3.0f;
  compute_emissions(&model, features, scores);
  CHECK(close_enough(scores[2], 7.0));
  CHECK(close_enough(scores[0], 0.0));
  CHECK(is_model_finite(&model));
  model.transition_weights[0] = -2.0;
  CHECK(is_model_finite(&model));
  model.transition_weights[0] = 0.0;

  uint8_t bases[4] = {0}, labels[4] = {INTERGENIC};
  float cwt[4 * CWT_CHANNELS] = {0};
  Sample sample = {.length = 1, .bases = bases, .cwt = cwt, .classes = labels};
  Sample invalid_sample = sample;
  invalid_sample.length = INT64_MAX;
  CHECK(train(&model, &invalid_sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  invalid_sample.length = -1;
  CHECK(train(&model, &invalid_sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  invalid_sample.length = 0;
  CHECK(train(&model, &invalid_sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  ModelObjectiveContext context = {&sample, 1, 1, 1, 0};
  double parameters[PARAMETERS] = {0};
  double gradient[PARAMETERS];
  double objective = model_objective(parameters, gradient, &context);
  CHECK(isfinite(objective) && close_enough(objective, 0.0));
  CHECK(fabs(gradient[0]) < 1e-12);
  /* A one-base path has zero NLL, isolating all three parameter groups' L2. */
  parameters[0] = 2;
  parameters[EMISSIONS] = -3;
  parameters[EMISSIONS + TRANSITIONS] = 4;
  objective = model_objective(parameters, gradient, &context);
  CHECK(close_enough(objective, 0.5 * L2 * 29));
  for (int i = 0; i < PARAMETERS; ++i)
    CHECK(close_enough(gradient[i], L2 * parameters[i]));
  Sample repeated[] = {sample, sample};
  ModelObjectiveContext duplicated = {repeated, 2, 2, 2, 0};
  CHECK(close_enough(model_objective(parameters, gradient, &duplicated), objective));
  for (int i = 0; i < PARAMETERS; ++i)
    CHECK(close_enough(gradient[i], L2 * parameters[i]));
  model.emission_weights[INPUT_CHANNELS] = NAN;
  CHECK(!is_model_finite(&model));
  model.emission_weights[INPUT_CHANNELS] = 0.0;
  CHECK(train(&model, &sample, 1, 1, 1) >= 0);
  invalid_sample = sample;
  invalid_sample.length = 4;
  labels[1] = CLASSES;
  CHECK(train(&model, &invalid_sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  labels[1] = INTERGENIC;
  bases[1] = NUCLEOTIDES + 1;
  CHECK(train(&model, &invalid_sample, 1, 1, 1) == LBFGS_INVALID_ARGUMENT);
  bases[1] = 0;
  labels[1] = PLUS_START_SINGLE_EXON;
  labels[2] = MINUS_START_SINGLE_EXON;
  CHECK(train(&model, &invalid_sample, 1, 1, 1) < 0);

  double x[1] = {0};
  int completed = -1;
  LbfgsStatus status = lbfgs_minimize(1, x, quadratic, NULL, 20, &completed);
  CHECK(status == LBFGS_CONVERGED || status == LBFGS_STOPPING_CRITERION ||
        status == LBFGS_MAX_ITERATIONS);
  CHECK(fabs(x[0] - 3.0) < 1e-4);
  CHECK(completed >= 0);
  CHECK(lbfgs_minimize(0, x, quadratic, NULL, 1, &completed) ==
        LBFGS_INVALID_ARGUMENT);

  char model_path[32];
  FILE *model_file = make_temp(model_path, "");
  CHECK(model_file != NULL);
  if (model_file)
    fclose(model_file);
  CHECK(save_model(model_path, &model) == 0);
  Model loaded = {0};
  CHECK(load_model(model_path, &loaded) == 0);
  CHECK(!memcmp(&model, &loaded, sizeof(model)));
  model_file = fopen(model_path, "r+b");
  CHECK(model_file != NULL);
  if (model_file) {
    char magic[4] = {0};
    CHECK(fread(magic, 1, sizeof(magic), model_file) == sizeof(magic));
    CHECK(!memcmp(magic, "SURF", sizeof(magic)));
    rewind(model_file);
    CHECK(fwrite("SURF", 1, 4, model_file) == 4);
    CHECK(fclose(model_file) == 0);
    loaded.weights[0] = 123.0;
    Model preserved = loaded;
    CHECK(!memcmp(&preserved, &loaded, sizeof(loaded)));
  }
  unlink(model_path);
}

static void test_crf_gradient(void) {
  uint8_t bases[] = {0, 1, NUCLEOTIDES, 3, 0, 1, 2, NUCLEOTIDES, 3};
  uint8_t labels[] = {0, 1, 1, 0, 0};
  float cwt[5 * CWT_CHANNELS];
  for (int i = 0; i < 5 * CWT_CHANNELS; i++)
    cwt[i] = (float)(i - 7) / 10;
  Sample sample = {.length = 5,
                   .bases = bases + 2,
                   .cwt = cwt,
                   .classes = labels,
                   .offset = 2,
                   .sequence_length = sizeof(bases)};
  ModelObjectiveContext context = {&sample, 1, 5, 1, 0};
  double parameters[PARAMETERS], gradient[PARAMETERS], scratch[PARAMETERS];
  for (int i = 0; i < PARAMETERS; i++)
    parameters[i] = (double)(i % 7 - 3) / 10;
  CHECK(isfinite(model_objective(parameters, gradient, &context)));
  Sample full = {.length = sizeof(bases), .bases = bases};
  double whole[TRANSITIONS], cut[TRANSITIONS];
  const double *bias = parameters + EMISSIONS, *pwm = bias + TRANSITIONS;
  for (int offset = 0; offset <= 4; offset += 2) {
    Sample window = sample;
    window.bases = bases + offset;
    window.offset = offset;
    for (int p = 0; p < 5; p += 4) {
      CHECK(compute_transitions(bias, pwm, &full, offset + p, whole));
      CHECK(compute_transitions(bias, pwm, &window, p, cut));
      CHECK(!memcmp(whole, cut, sizeof(whole)));
    }
  }
  for (int i = 0; i < PARAMETERS; i++) {
    double original = parameters[i];
    parameters[i] = original + 1e-6;
    double plus = model_objective(parameters, scratch, &context);
    parameters[i] = original - 1e-6;
    double minus = model_objective(parameters, scratch, &context);
    parameters[i] = original;
    CHECK(isfinite(plus) && isfinite(minus) &&
          fabs(gradient[i] - (plus - minus) / 2e-6) < 1e-7);
  }
  Sample samples[] = {sample,
                      {0},
                      {4, bases, cwt, labels, 0, 0},
                      {1, bases, cwt, labels, 0, 0}};
  ModelObjectiveContext parallel = {samples, 4, 10, 1, 0};
  double serial_loss = model_objective(parameters, gradient, &parallel);
  CHECK(parallel.failure_status == 0 && isfinite(serial_loss));
  parallel.threads = 4;
  CHECK(close_enough(serial_loss,
                     model_objective(parameters, scratch, &parallel)));
  CHECK(parallel.failure_status == 0);
  for (int i = 0; i < PARAMETERS; ++i)
    CHECK(close_enough(gradient[i], scratch[i]));
  Model serial, threaded;
  memcpy(serial.weights, parameters, sizeof(parameters));
  threaded = serial;
  CHECK(train(&serial, samples, 4, 2, 1) >= 0);
  CHECK(train(&threaded, samples, 4, 2, 4) >= 0);
  for (int i = 0; i < PARAMETERS; ++i)
    CHECK(fabs(serial.weights[i] - threaded.weights[i]) < 1e-9);
}

static void test_multiexon_crf_gradient(void) {
  enum { LENGTH = 16, CUT = 6 };
  const char *symbols[] = {"BSSIIBBIIBBIIEEB", "BeeiiBBiiBBiissB"};
  uint8_t bases[LENGTH], labels[2][LENGTH];
  float cwt[LENGTH * CWT_CHANNELS];
  for (int p = 0; p < LENGTH; ++p)
    bases[p] = (uint8_t)(p % NUCLEOTIDES);
  for (int i = 0; i < LENGTH * CWT_CHANNELS; ++i)
    cwt[i] = (float)(i % 7 - 3) / 10;
  Sample samples[4];
  for (int strand = 0; strand < 2; ++strand) {
    label_path(symbols[strand], labels[strand]);
    CHECK(labels[strand][CUT - 1] == INTERGENIC &&
          labels[strand][CUT] == INTERGENIC);
    /* Cut inside an internal exon, retaining DNA outside each window. */
    samples[2 * strand] =
        (Sample){CUT, bases, cwt, labels[strand], 0, LENGTH};
    samples[2 * strand + 1] =
        (Sample){LENGTH - CUT, bases + CUT, cwt + CUT * CWT_CHANNELS,
                 labels[strand] + CUT, CUT, LENGTH};
  }
  ModelObjectiveContext context = {samples, 4, 2 * LENGTH, 1, 0};
  double parameters[PARAMETERS], gradient[PARAMETERS], scratch[PARAMETERS];
  for (int i = 0; i < PARAMETERS; ++i)
    parameters[i] = (double)(i % 7 - 3) / 100;
  double objective = model_objective(parameters, gradient, &context);
  CHECK(context.failure_status == 0 && isfinite(objective) && objective > 0);
  for (int i = 0; i < PARAMETERS; ++i)
    CHECK(isfinite(gradient[i]));

  Sample full = {.length = LENGTH, .bases = bases};
  double whole[TRANSITIONS], cut[TRANSITIONS];
  const double *bias = parameters + EMISSIONS, *pwm = bias + TRANSITIONS;
  CHECK(compute_transitions(bias, pwm, &full, CUT + 1, whole));
  CHECK(compute_transitions(bias, pwm, &samples[1], 1, cut));
  CHECK(!memcmp(whole, cut, sizeof(whole)));

  /* Check only parameters tied to background and its intron boundaries. */
  const int selected[] = {
      INTERGENIC * (INPUT_CHANNELS + 1) + INPUT_CHANNELS,
      EMISSIONS + PLUS_INTRON * CLASSES + INTERGENIC,
      EMISSIONS + INTERGENIC * CLASSES + PLUS_INTRON,
      EMISSIONS + MINUS_INTRON * CLASSES + INTERGENIC,
      EMISSIONS + INTERGENIC * CLASSES + MINUS_INTRON,
      EMISSIONS + TRANSITIONS + PLUS_INTRON * BOUNDARY_FEATURES +
          BOUNDARY_RADIUS * NUCLEOTIDES + 1,
      EMISSIONS + TRANSITIONS + (CLASSES + INTERGENIC) * BOUNDARY_FEATURES +
          BOUNDARY_RADIUS * NUCLEOTIDES + 1};
  for (size_t c = 0; c < sizeof(selected) / sizeof(*selected); ++c) {
    int i = selected[c];
    double original = parameters[i];
    parameters[i] = original + 1e-6;
    double plus = model_objective(parameters, scratch, &context);
    CHECK(context.failure_status == 0);
    parameters[i] = original - 1e-6;
    double minus = model_objective(parameters, scratch, &context);
    CHECK(context.failure_status == 0);
    parameters[i] = original;
    CHECK(isfinite(plus) && isfinite(minus) &&
          fabs(gradient[i] - (plus - minus) / 2e-6) < 1e-7);
  }
}

static void test_prediction(void) {
  Model model = {0};
  char bases[] = "A";
  Sequence sequence = {(char *)"tiny", bases, 1, NULL};
  FILE *output = tmpfile();
  CHECK(output != NULL);
  if (!output)
    return;
  CHECK(predict(&model, &sequence, 1, output, NULL) == 0);
  fclose(output);
}

static void test_prediction_parallel_seam(void) {
  enum { LENGTH = 16384 + 37 };
  char bases[LENGTH + 1] = {0}, text[2][1024] = {{0}};
  memset(bases, 'C', LENGTH);
  Sequence sequences[] = {{(char *)"seam", bases, LENGTH, NULL},
                          {(char *)"short", bases, 29, NULL}};
  Model model = {0};
  model.emission_weights[PLUS_START_SINGLE_EXON * (INPUT_CHANNELS + 1) +
                         INPUT_CHANNELS] = 2;
  model.transition_weights[PLUS_START_SINGLE_EXON * CLASSES] = -1;
  model
      .boundary_weights[(CLASSES + PLUS_START_SINGLE_EXON) * BOUNDARY_FEATURES +
                        BOUNDARY_RADIUS * NUCLEOTIDES + 1] = 0.1;
  for (int pass = 0; pass < 2; ++pass) {
    omp_set_num_threads(pass ? 4 : 1);
    FILE *output = tmpfile();
    CHECK(output != NULL);
    if (!output)
      continue;
    CHECK(predict(&model, sequences, 2, output, NULL) == 2);
    rewind(output);
    CHECK(fread(text[pass], 1, sizeof(text[pass]) - 1, output) > 0);
    CHECK(!ferror(output) && feof(output));
    fclose(output);
  }
  omp_set_num_threads(1);
  CHECK(!strcmp(text[0], text[1]));
  CHECK(!strcmp(
      text[0],
      "seam\tanno\tgene\t2\t16420\t.\t+\t.\tID=gene1\n"
      "seam\tanno\tmRNA\t2\t16420\t.\t+\t.\tID=gene1.t1;Parent=gene1\n"
      "seam\tanno\texon\t2\t16420\t.\t+\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
      "short\tanno\tgene\t2\t28\t.\t+\t.\tID=gene2\n"
      "short\tanno\tmRNA\t2\t28\t.\t+\t.\tID=gene2.t1;Parent=gene2\n"
      "short\tanno\texon\t2\t28\t.\t+\t.\tID=gene2.t1.exon1;Parent=gene2."
      "t1\n"));
}

int main(void) {
  omp_set_dynamic(0);
  omp_set_num_threads(1);
  test_utils();
  test_fasta_and_gff_input();
  test_cwt_window_independence();
  test_intervals_and_gff_output();
  test_multiexon_gff();
  test_multiexon_gff_input();
  test_gff_roundtrip();
  test_gff_neighboring_genes();
  test_gff_invalid_paths();
  test_gff_label_conflicts();
  test_decoder_and_grammar();
  test_crf_and_lbfgs();

  test_crf_gradient();
  test_multiexon_crf_gradient();
  test_prediction();
  test_prediction_parallel_seam();
  fprintf(stderr, "%d checks, %d failures\n", checks, failures);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
