/**
 * @file test.c
 * @brief Unit and pipeline tests for anno. 
 */
#define ANNO_NO_MAIN
#include "anno.c"
#include "unet.c" /* static internals for the gradient test */

#include <sys/wait.h>
#include <unistd.h>

static int failures = 0, checks = 0;

/** @brief Return a deterministic random integer below the given bound. */
static unsigned random_below(uint64_t *state, unsigned bound) {
  return (unsigned)(anno_random_next(state) % bound);
}


static float *test_sequence_power(const AnnoSequence *sequence) {
  float *power =
      xcalloc((size_t)sequence->length * ANNO_CWT_CHANNELS, sizeof(float));
  for (int64_t start = 0; start < sequence->length; start += ANNO_CWT_WINDOW) {
    int64_t end = start + ANNO_CWT_WINDOW;
    if (end > sequence->length)
      end = sequence->length;
    float values[ANNO_CWT_WINDOW * ANNO_CWT_CHANNELS * ANNO_COMPLEX_COMPONENTS];
    anno_cwt(sequence->bases, sequence->length, start, end, values);
    for (int64_t i = 0; i < (end - start) * ANNO_CWT_CHANNELS; i++)
      power[(size_t)start * ANNO_CWT_CHANNELS + (size_t)i] =
          values[ANNO_COMPLEX_COMPONENTS * i] *
              values[ANNO_COMPLEX_COMPONENTS * i] +
          values[ANNO_COMPLEX_COMPONENTS * i + 1] *
              values[ANNO_COMPLEX_COMPONENTS * i + 1];
  }
  return power;
}

/** @brief Record one check. */
#define CHECK(condition)                                                       \
  do {                                                                         \
    checks++;                                                                  \
    if (!(condition)) {                                                        \
      failures++;                                                              \
      fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__,         \
              #condition);                                                     \
    }                                                                          \
  } while (0)

/** @brief Write text to a file. */
static void write_text(const char *path, const char *text) {
  FILE *file = fopen(path, "w");
  if (!file || fputs(text, file) < 0 || fclose(file)) {
    perror(path);
    exit(EXIT_FAILURE);
  }
}

/** @brief Read a whole file into a NUL-terminated buffer (caller frees). */
static char *read_text(FILE *file) {
  long size = ftell(file);
  char *text = xcalloc((size_t)size + 1, 1);
  rewind(file);
  if (fread(text, 1, (size_t)size, file) != (size_t)size)
    text[0] = 0;
  return text;
}

/* ---- hand-made gene pair ------------------------------------------------ */

/*
 * Block B (500 bp, background 'C'), plus-strand transcript, 0-based [a,b):
 *   exons   [40,70) [100,150) [210,300) [360,420)  (first intron is in the
 * 5'UTR) 5'UTR   [40,70) [100,121) CDS     [121,150) [210,300) [360,391)  =
 * 29+90+31 = 150 bp = ATG + 48 GCC + TAA 3'UTR   [391,420) Introns start with
 * GT and end with AG. Genome = B + 1000 'C' + revcomp(B), so the same gene
 * appears on the minus strand at 1999-p.
 */
#define TEST_LENGTH 2000
static const char *const test_gff =
    "##gff-version 3\n"
    "t\tanno\tgene\t41\t420\t.\t+\t.\tID=gene1\n"
    "t\tanno\tmRNA\t41\t420\t.\t+\t.\tID=gene1.t1;Parent=gene1\n"
    "t\tanno\texon\t41\t70\t.\t+\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
    "t\tanno\texon\t101\t150\t.\t+\t.\tID=gene1.t1.exon2;Parent=gene1.t1\n"
    "t\tanno\texon\t211\t300\t.\t+\t.\tID=gene1.t1.exon3;Parent=gene1.t1\n"
    "t\tanno\texon\t361\t420\t.\t+\t.\tID=gene1.t1.exon4;Parent=gene1.t1\n"
    "t\tanno\tfive_prime_UTR\t41\t70\t.\t+\t.\tID=gene1.t1.utr5p1;Parent=gene1."
    "t1\n"
    "t\tanno\tfive_prime_UTR\t101\t121\t.\t+\t.\tID=gene1.t1.utr5p2;Parent="
    "gene1.t1\n"
    "t\tanno\tCDS\t122\t150\t.\t+\t0\tID=gene1.t1.cds1;Parent=gene1.t1\n"
    "t\tanno\tCDS\t211\t300\t.\t+\t1\tID=gene1.t1.cds2;Parent=gene1.t1\n"
    "t\tanno\tCDS\t361\t391\t.\t+\t1\tID=gene1.t1.cds3;Parent=gene1.t1\n"
    "t\tanno\tthree_prime_UTR\t392\t420\t.\t+\t.\tID=gene1.t1.utr3p1;Parent="
    "gene1.t1\n"
    "t\tanno\tgene\t1581\t1960\t.\t-\t.\tID=gene2\n"
    "t\tanno\tmRNA\t1581\t1960\t.\t-\t.\tID=gene2.t1;Parent=gene2\n"
    "t\tanno\texon\t1581\t1640\t.\t-\t.\tID=gene2.t1.exon1;Parent=gene2.t1\n"
    "t\tanno\texon\t1701\t1790\t.\t-\t.\tID=gene2.t1.exon2;Parent=gene2.t1\n"
    "t\tanno\texon\t1851\t1900\t.\t-\t.\tID=gene2.t1.exon3;Parent=gene2.t1\n"
    "t\tanno\texon\t1931\t1960\t.\t-\t.\tID=gene2.t1.exon4;Parent=gene2.t1\n"
    "t\tanno\tfive_prime_UTR\t1880\t1900\t.\t-\t.\tID=gene2.t1.utr5p1;Parent="
    "gene2.t1\n"
    "t\tanno\tfive_prime_UTR\t1931\t1960\t.\t-\t.\tID=gene2.t1.utr5p2;Parent="
    "gene2.t1\n"
    "t\tanno\tCDS\t1610\t1640\t.\t-\t1\tID=gene2.t1.cds1;Parent=gene2.t1\n"
    "t\tanno\tCDS\t1701\t1790\t.\t-\t1\tID=gene2.t1.cds2;Parent=gene2.t1\n"
    "t\tanno\tCDS\t1851\t1879\t.\t-\t0\tID=gene2.t1.cds3;Parent=gene2.t1\n"
    "t\tanno\tthree_prime_UTR\t1581\t1609\t.\t-\t.\tID=gene2.t1.utr3p1;Parent="
    "gene2.t1\n";

/** @brief Build the 2000-bp test genome described above (caller frees). */
static char *test_genome(void) {
  static const int exons[4][2] = {{40, 70}, {100, 150}, {210, 300}, {360, 420}};
  char block[501], *genome = xcalloc(TEST_LENGTH + 1, 1);
  memset(block, 'C', 500);
  block[500] = 0;
  char cds[151];
  memcpy(cds, "ATG", 3);
  for (int c = 1; c < 49; c++)
    memcpy(cds + 3 * c, "GCC", 3);
  memcpy(cds + 147, "TAA", 3);
  int at = 0; /* spliced coordinate */
  for (int e = 0; e < 4; e++)
    for (int p = exons[e][0]; p < exons[e][1]; p++, at++)
      if (p >= 121 && p < 391)
        block[p] = cds[at - 51];
  for (int e = 0; e < 3; e++) {
    block[exons[e][1]] = 'G';
    block[exons[e][1] + 1] = 'T';
    block[exons[e + 1][0] - 2] = 'A';
    block[exons[e + 1][0] - 1] = 'G';
  }
  memset(genome, 'C', TEST_LENGTH);
  memcpy(genome, block, 500);
  for (int p = 0; p < 500; p++)
    genome[TEST_LENGTH - 1 - p] = complement(block[p]);
  return genome;
}

/** @brief FASTA + GFF of the test genome on disk; returns the sequence. */
static AnnoSequence *test_sequences(const char *fasta, const char *gff) {
  char *genome = test_genome();
  size_t size = strlen(genome) + 16;
  char *text = xcalloc(size, 1);
  snprintf(text, size, ">t\n%s\n", genome);
  write_text(fasta, text);
  write_text(gff, test_gff);
  free(text);
  free(genome);
  int count = 0;
  AnnoSequence *sequences = anno_read_fasta(fasta, &count);
  CHECK(count == 1 && sequences[0].length == TEST_LENGTH);
  anno_read_labels(gff, sequences, count);
  return sequences;
}

/** @brief GFF -> label vector checks at hand-picked positions. */
static void test_labels(const AnnoSequence *s) {
  const uint8_t *y = s->labels;
  CHECK(ANNO_CLASSES == 7);
  CHECK(ANNO_INTERGENIC == 0 && ANNO_PLUS_START_SINGLE_EXON == 1 &&
        ANNO_PLUS_END_EXON == 2 && ANNO_PLUS_INTRON == 3 &&
        ANNO_MINUS_START_SINGLE_EXON == 4 && ANNO_MINUS_END_EXON == 5 &&
        ANNO_MINUS_INTRON == 6);
  CHECK(ANNO_PLUS_INTERNAL_EXON == 7 && ANNO_MINUS_INTERNAL_EXON == 8);
  CHECK(y[39] == ANNO_INTERGENIC && y[40] == ANNO_PLUS_START_SINGLE_EXON);
  CHECK(y[69] == ANNO_PLUS_START_SINGLE_EXON);
  CHECK(y[70] == ANNO_PLUS_INTRON && y[99] == ANNO_PLUS_INTRON);
  CHECK(y[120] == ANNO_INTERGENIC && y[121] == ANNO_INTERGENIC);
  CHECK(y[390] == ANNO_PLUS_END_EXON && y[391] == ANNO_PLUS_END_EXON);
  CHECK(y[419] == ANNO_PLUS_END_EXON && y[420] == ANNO_INTERGENIC);
  CHECK(y[1579] == ANNO_INTERGENIC && y[1580] == ANNO_MINUS_END_EXON);
  CHECK(y[1609] == ANNO_MINUS_END_EXON && y[1640] == ANNO_MINUS_INTRON);
  CHECK(y[1879] == ANNO_INTERGENIC &&
        y[1959] == ANNO_MINUS_START_SINGLE_EXON);
  CHECK(y[1960] == ANNO_INTERGENIC);
  for (int i = 0; i < TEST_LENGTH; i++) {
    int expected = i < 40     ? ANNO_INTERGENIC
                   : i < 70   ? ANNO_PLUS_START_SINGLE_EXON
                   : i < 100  ? ANNO_PLUS_INTRON
                   : i < 150  ? ANNO_INTERGENIC
                   : i < 210  ? ANNO_PLUS_INTRON
                   : i < 300  ? ANNO_INTERGENIC
                   : i < 360  ? ANNO_PLUS_INTRON
                   : i < 420  ? ANNO_PLUS_END_EXON
                   : i < 1580 ? ANNO_INTERGENIC
                   : i < 1640 ? ANNO_MINUS_END_EXON
                   : i < 1700 ? ANNO_MINUS_INTRON
                   : i < 1790 ? ANNO_INTERGENIC
                   : i < 1850 ? ANNO_MINUS_INTRON
                   : i < 1900 ? ANNO_INTERGENIC
                   : i < 1930 ? ANNO_MINUS_INTRON
                   : i < 1960 ? ANNO_MINUS_START_SINGLE_EXON
                              : ANNO_INTERGENIC;
    CHECK(y[i] == expected && y[i] < ANNO_CLASSES);
  }
}

/** @brief Reverse-complement augmentation preserves transcriptional roles. */
static void test_strand_swap(void) {
  static const uint8_t swapped[] = {
      ANNO_INTERGENIC, ANNO_MINUS_START_SINGLE_EXON, ANNO_MINUS_END_EXON,
      ANNO_MINUS_INTRON, ANNO_PLUS_START_SINGLE_EXON, ANNO_PLUS_END_EXON,
      ANNO_PLUS_INTRON};
  for (int label = 0; label < ANNO_CLASSES; label++) {
    CHECK(ann_swap_strand(label) == swapped[label]);
    CHECK(ann_swap_strand(ann_swap_strand(label)) == label);
  }
}

/** @brief Print the GFF records of a path into a string (caller frees). */
static char *gff_of_path(const AnnoSequence *s, const uint8_t *path) {
  FILE *file = tmpfile();
  unsigned long genes = 0;
  fputs("##gff-version 3\n", file);
  CHECK(anno_write_gff(file, s, path, &genes) == 2 && genes == 2);
  char *text = read_text(file);
  fclose(file);
  return text;
}

/** @brief Restore known intron-bounded exons before writing decoded GFF. */
static void test_round_trip(const AnnoSequence *s) {
  uint8_t *expected = xcalloc((size_t)s->length, 1);
  memcpy(expected, s->labels, (size_t)s->length);
  static const int internal[4][2] = {
      {100, 150}, {210, 300}, {1700, 1790}, {1850, 1900}};
  for (int e = 0; e < 4; e++)
    for (int i = internal[e][0]; i < internal[e][1]; i++) {
      CHECK(s->labels[i] == ANNO_INTERGENIC);
      expected[i] = e < 2 ? ANNO_PLUS_INTERNAL_EXON : ANNO_MINUS_INTERNAL_EXON;
    }
  char *direct = gff_of_path(s, expected);
  CHECK(!strcmp(direct, test_gff));
  if (strcmp(direct, test_gff))
    fprintf(stderr, "got:\n%s", direct);
  float *probability = xcalloc((size_t)s->length * ANNO_CLASSES, sizeof(float));
  for (int64_t i = 0; i < s->length; i++)
    for (int j = 0; j < ANNO_CLASSES; j++)
      probability[i * ANNO_CLASSES + j] =
          j == s->labels[i] ? 0.9f : 0.1f / (ANNO_CLASSES - 1);
  uint8_t *path = xcalloc((size_t)s->length, 1);
  CHECK(anno_decode(s->length, probability, path) == 0);
  CHECK(!memcmp(path, expected, (size_t)s->length));
  char *decoded = gff_of_path(s, path);
  CHECK(!strcmp(decoded, test_gff));
  free(direct);
  free(decoded);
  free(path);
  free(probability);
  free(expected);
}

/** @brief Grammar constraints on a probability matrix. */
static void test_grammar(void) {
  Grammar grammar;
  grammar_build(&grammar);
  CHECK(grammar.count ==
        1 + 2 * (3 * ANNO_MIN_EXON_LENGTH + ANNO_MIN_INTRON_LENGTH));
  uint8_t reachable[GRAMMAR_MAX_STATES] = {0};
  reachable[0] = 1;
  for (int pass = 0; pass < grammar.count; ++pass)
    for (int state = 0; state < grammar.count; ++state)
      for (int k = 0; k < grammar.predecessors[state]; ++k)
        if (reachable[grammar.from[state][k]])
          reachable[state] = 1;
  int count = 0;
  for (int state = 0; state < grammar.count; ++state)
    count += reachable[state] != 0;
  CHECK(count == grammar.count);
}

/** @brief Recognize a decoded path using the duration-expanded grammar. */
static int grammar_accepts(const uint8_t *path, int n) {
  Grammar g;
  grammar_build(&g);
  uint8_t active[GRAMMAR_MAX_STATES] = {0};
  if (!n || path[0] != ANNO_INTERGENIC)
    return 0;
  active[0] = 1;
  for (int i = 1; i < n; i++) {
    uint8_t next[GRAMMAR_MAX_STATES] = {0};
    for (int s = 0; s < g.count; s++)
      if (g.cls[s] == path[i])
        for (int k = 0; k < g.predecessors[s]; k++)
          next[s] |= active[g.from[s][k]];
    memcpy(active, next, sizeof(active));
  }
  return active[0] != 0;
}

/** @brief Recognize a short list of constant-label runs. */
static int grammar_accepts_runs(const uint8_t *labels, const int *lengths,
                                int count) {
  uint8_t path[512];
  int n = 0;
  for (int r = 0; r < count; r++) {
    if (lengths[r] < 0 || lengths[r] > (int)sizeof(path) - n)
      return 0;
    memset(path + n, labels[r], (size_t)lengths[r]);
    n += lengths[r];
  }
  return grammar_accepts(path, n);
}

/** @brief Terminal roles, strand isolation and all duration boundaries. */
static void test_grammar_paths(void) {
  CHECK(ANNO_MIN_EXON_LENGTH == 6 && ANNO_MIN_INTRON_LENGTH == 30);
  for (int strand = 0; strand < 2; strand++) {
    uint8_t start = strand ? ANNO_MINUS_START_SINGLE_EXON
                           : ANNO_PLUS_START_SINGLE_EXON;
    uint8_t end = strand ? ANNO_MINUS_END_EXON : ANNO_PLUS_END_EXON;
    uint8_t internal = strand ? ANNO_MINUS_INTERNAL_EXON
                              : ANNO_PLUS_INTERNAL_EXON;
    uint8_t intron = strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON;
    uint8_t single[] = {ANNO_INTERGENIC, start, ANNO_INTERGENIC};
    int single_lengths[] = {10, 6, 10};
    CHECK(grammar_accepts_runs(single, single_lengths, 3));
    single_lengths[1] = 5;
    CHECK(!grammar_accepts_runs(single, single_lengths, 3));
    single_lengths[1] = 6;
    single[1] = end;
    CHECK(!grammar_accepts_runs(single, single_lengths, 3));
    single[1] = internal;
    CHECK(!grammar_accepts_runs(single, single_lengths, 3));
    single[1] = intron;
    single_lengths[1] = 30;
    CHECK(!grammar_accepts_runs(single, single_lengths, 3));

    uint8_t labels[] = {ANNO_INTERGENIC, strand ? end : start, intron,
                        internal, intron, strand ? start : end,
                        ANNO_INTERGENIC};
    int lengths[] = {10, 6, 30, 6, 30, 6, 10};
    CHECK(grammar_accepts_runs(labels, lengths, 7));
    for (int r = 1; r <= 5; r++) {
      lengths[r]--;
      CHECK(!grammar_accepts_runs(labels, lengths, 7));
      lengths[r] += 2;
      CHECK(grammar_accepts_runs(labels, lengths, 7));
      lengths[r]--;
    }
    labels[3] = ANNO_INTERGENIC; /* training targets are not decoded paths */
    CHECK(!grammar_accepts_runs(labels, lengths, 7));
    labels[3] = internal;
    labels[4] = strand ? ANNO_PLUS_INTRON : ANNO_MINUS_INTRON;
    CHECK(!grammar_accepts_runs(labels, lengths, 7));
    labels[4] = intron;
    labels[5] = strand ? ANNO_PLUS_START_SINGLE_EXON : ANNO_MINUS_END_EXON;
    CHECK(!grammar_accepts_runs(labels, lengths, 7));

    uint8_t pair[] = {ANNO_INTERGENIC, strand ? end : start, intron,
                      strand ? start : end, ANNO_INTERGENIC};
    int pair_lengths[] = {10, 6, 30, 6, 10};
    CHECK(grammar_accepts_runs(pair, pair_lengths, 5));
    pair[3] = pair[1]; /* no second start on plus, or second end on minus */
    CHECK(!grammar_accepts_runs(pair, pair_lengths, 5));
    pair[3] = strand ? start : end;
    pair[1] = pair[3]; /* reversed terminal order */
    CHECK(!grammar_accepts_runs(pair, pair_lengths, 5));
    pair[1] = strand ? end : start;
    pair_lengths[2] = 0; /* terminal exons cannot touch without an intron */
    CHECK(!grammar_accepts_runs(pair, pair_lengths, 5));
    pair_lengths[2] = 30;
    pair_lengths[3] = 0; /* a transcript cannot finish in an intron */
    CHECK(!grammar_accepts_runs(pair, pair_lengths, 5));

    uint8_t mixed[] = {ANNO_INTERGENIC, start,
                       strand ? ANNO_PLUS_START_SINGLE_EXON
                              : ANNO_MINUS_START_SINGLE_EXON,
                       ANNO_INTERGENIC};
    int mixed_lengths[] = {10, 6, 6, 10};
    CHECK(!grammar_accepts_runs(mixed, mixed_lengths, 4));
  }
}

static void test_viterbi(void) {
  enum { N = 60 };
  float p[N * ANNO_CLASSES];
  uint8_t path[N];
  for (int strand = 0; strand < 2; strand++) {
    int start = strand ? ANNO_MINUS_START_SINGLE_EXON
                       : ANNO_PLUS_START_SINGLE_EXON;
    /* A 4-bp single-exon signal is widened to the 6-bp minimum. */
    for (int i = 0; i < N; i++)
      for (int j = 0; j < ANNO_CLASSES; j++)
        p[i * ANNO_CLASSES + j] =
            j == (i >= 20 && i < 24 ? start : ANNO_INTERGENIC)
                ? 0.99f : 0.01f / (ANNO_CLASSES - 1);
    CHECK(anno_decode(N, p, path) == 0);
    CHECK(grammar_accepts(path, N));
    int exon = 0;
    for (int i = 0; i < N; i++)
      exon += path[i] == start;
    CHECK(exon == ANNO_MIN_EXON_LENGTH);
    /* A strong 30-bp single exon is decoded exactly on either strand. */
    for (int i = 0; i < N; i++)
      for (int j = 0; j < ANNO_CLASSES; j++)
        p[i * ANNO_CLASSES + j] =
            j == (i >= 20 && i < 50 ? start : ANNO_INTERGENIC)
                ? 0.99f : 0.01f / (ANNO_CLASSES - 1);
    CHECK(anno_decode(N, p, path) == 0);
    for (int i = 0; i < N; i++)
      CHECK(path[i] == (i >= 20 && i < 50 ? start : ANNO_INTERGENIC));
  }
  /* NaN emissions are treated as zero probability, not propagated. */
  p[5 * ANNO_CLASSES] = NAN;
  CHECK(anno_decode(N, p, path) == 0);
  CHECK(grammar_accepts(path, N));
}

/** @brief Favour a decoded path using only the seven model targets. */
static void favour_path(const uint8_t *path, int n, float *p) {
  for (int i = 0; i < n; i++) {
    int target = path[i];
    if (target == ANNO_PLUS_INTERNAL_EXON || target == ANNO_MINUS_INTERNAL_EXON)
      target = ANNO_INTERGENIC;
    for (int j = 0; j < ANNO_CLASSES; j++)
      p[i * ANNO_CLASSES + j] =
          j == target ? 0.9f : 0.1f / (ANNO_CLASSES - 1);
  }
}

/** @brief Single-exon labels and exact GFF output on both strands. */
static void test_single_exons(void) {
  enum { N = 70 };
  static const char gff[] =
      "##gff-version 3\n"
      "s\tanno\tgene\t11\t19\t.\t+\t.\tID=gene1\n"
      "s\tanno\tmRNA\t11\t19\t.\t+\t.\tID=gene1.t1;Parent=gene1\n"
      "s\tanno\texon\t11\t19\t.\t+\t.\tID=gene1.t1.exon1;Parent=gene1.t1\n"
      "s\tanno\tCDS\t11\t19\t.\t+\t0\tID=gene1.t1.cds1;Parent=gene1.t1\n"
      "s\tanno\tgene\t41\t49\t.\t-\t.\tID=gene2\n"
      "s\tanno\tmRNA\t41\t49\t.\t-\t.\tID=gene2.t1;Parent=gene2\n"
      "s\tanno\texon\t41\t49\t.\t-\t.\tID=gene2.t1.exon1;Parent=gene2.t1\n"
      "s\tanno\tCDS\t41\t49\t.\t-\t0\tID=gene2.t1.cds1;Parent=gene2.t1\n";
  char bases[N + 1];
  memset(bases, 'C', N);
  bases[N] = 0;
  memcpy(bases + 10, "ATGAAATAA", 9);
  memcpy(bases + 40, "TTATTTCAT", 9);
  AnnoSequence s = {"s", bases, N, NULL};
  write_text("test-single.gff3", gff);
  anno_read_labels("test-single.gff3", &s, 1);
  uint8_t expected[N] = {0}, path[N];
  memset(expected + 10, ANNO_PLUS_START_SINGLE_EXON, 9);
  memset(expected + 40, ANNO_MINUS_START_SINGLE_EXON, 9);
  CHECK(!memcmp(s.labels, expected, N));
  char *direct = gff_of_path(&s, expected);
  CHECK(!strcmp(direct, gff));
  float p[N * ANNO_CLASSES];
  favour_path(expected, N, p);
  CHECK(anno_decode(N, p, path) == 0);
  CHECK(!memcmp(path, expected, N));
  char *decoded = gff_of_path(&s, path);
  CHECK(!strcmp(decoded, gff));
  free(direct);
  free(decoded);
  free(s.labels);
  remove("test-single.gff3");
}

/** @brief Intergenic internal targets still occupy the whole transcript span. */
static void test_nested_transcripts(const AnnoSequence *original) {
  size_t size = strlen(test_gff) + 256;
  char *gff = xcalloc(size, 1);
  for (int strand = 0; strand < 2; strand++) {
    snprintf(gff, size,
             "%st\ta\texon\t111\t130\t.\t%c\t.\tParent=nested1\n"
             "t\ta\texon\t1861\t1880\t.\t%c\t.\tParent=nested2\n",
             test_gff, strand ? '-' : '+', strand ? '+' : '-');
    write_text("test-nested.gff3", gff);
    AnnoSequence s = {original->name, original->bases, original->length, NULL};
    anno_read_labels("test-nested.gff3", &s, 1);
    CHECK(!memcmp(s.labels, original->labels, (size_t)s.length));
    for (int i = 110; i < 130; i++)
      CHECK(s.labels[i] == ANNO_INTERGENIC);
    for (int i = 1860; i < 1880; i++)
      CHECK(s.labels[i] == ANNO_INTERGENIC);
    free(s.labels);
  }
  free(gff);
  remove("test-nested.gff3");
}

/** @brief Invalid signals must decode into complete, single-strand genes. */
static void test_invalid_signals(void) {
  enum { N = 180 };
  for (int strand = 0; strand < 2; strand++) {
    int end = strand ? ANNO_MINUS_END_EXON : ANNO_PLUS_END_EXON;
    int first = strand ? ANNO_MINUS_END_EXON : ANNO_PLUS_START_SINGLE_EXON;
    int intron = strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON;
    int wrong_last = strand ? ANNO_PLUS_START_SINGLE_EXON : ANNO_MINUS_END_EXON;
    for (int signal = 0; signal < 4; signal++) {
      uint8_t truth[N] = {0}, path[N];
      if (signal == 0)
        memset(truth + 20, end, 30); /* an end exon cannot stand alone */
      else if (signal == 1)
        memset(truth + 50, intron, 60); /* an intron needs both terminal exons */
      else if (signal == 2) {
        memset(truth + 20, first, 30);
        memset(truth + 50, intron, 60);
        memset(truth + 110, wrong_last, 30);
      } else {
        memset(truth + 20, strand ? ANNO_MINUS_START_SINGLE_EXON
                                 : ANNO_PLUS_START_SINGLE_EXON, 30);
        memset(truth + 50, strand ? ANNO_PLUS_START_SINGLE_EXON
                                 : ANNO_MINUS_START_SINGLE_EXON, 30);
      }
      CHECK(!grammar_accepts(truth, N));
      float p[N * ANNO_CLASSES];
      favour_path(truth, N, p);
      CHECK(anno_decode(N, p, path) == 0);
      CHECK(grammar_accepts(path, N));
      CHECK(memcmp(path, truth, N) != 0);
    }
  }
}

/** @brief Decode all exon types and introns exactly at their minimum lengths. */
static void test_minimum_decode(void) {
  enum { N = 98 };
  for (int strand = 0; strand < 2; strand++) {
    uint8_t truth[N] = {0}, path[N];
    memset(truth + 10, strand ? ANNO_MINUS_END_EXON
                             : ANNO_PLUS_START_SINGLE_EXON, 6);
    memset(truth + 16, strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON, 30);
    memset(truth + 46, strand ? ANNO_MINUS_INTERNAL_EXON
                             : ANNO_PLUS_INTERNAL_EXON, 6);
    memset(truth + 52, strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON, 30);
    memset(truth + 82, strand ? ANNO_MINUS_START_SINGLE_EXON
                             : ANNO_PLUS_END_EXON, 6);
    float p[N * ANNO_CLASSES];
    favour_path(truth, N, p);
    CHECK(grammar_accepts(truth, N));
    CHECK(anno_decode(N, p, path) == 0);
    CHECK(!memcmp(path, truth, N));
  }
}

/** @brief Start/stop codons exactly at exon edges (UTR intron next to CDS). */
static void test_utr_intron_edges(void) {
  enum { N = 280 };
  /* genomic order: inter | exon [20,40) | intron [40,100) |
   * exon [100,160) | intron [160,220) | exon [220,250) | inter. */
  for (int strand = 0; strand < 2; strand++) {
    uint8_t truth[N], path[N];
    for (int i = 0; i < N; i++)
      truth[i] = i < 20    ? ANNO_INTERGENIC
                 : i < 40  ? (strand ? ANNO_MINUS_END_EXON
                                     : ANNO_PLUS_START_SINGLE_EXON)
                 : i < 100 ? (strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON)
                 : i < 160 ? (strand ? ANNO_MINUS_INTERNAL_EXON
                                     : ANNO_PLUS_INTERNAL_EXON)
                 : i < 220 ? (strand ? ANNO_MINUS_INTRON : ANNO_PLUS_INTRON)
                 : i < 250 ? (strand ? ANNO_MINUS_START_SINGLE_EXON
                                     : ANNO_PLUS_END_EXON)
                           : ANNO_INTERGENIC;
    float p[N * ANNO_CLASSES];
    favour_path(truth, N, p);
    CHECK(anno_decode(N, p, path) == 0);
    CHECK(!memcmp(path, truth, N));
  }
}

/** @brief Backpropagation equals central finite differences of the loss. */
static void test_gradient(void) {
  enum { B = 2 };
  size_t xs = (size_t)B * ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW,
         ys = (size_t)B * ANNO_CLASSES * ANNO_UNET_WINDOW;
  float *x = xcalloc(xs, sizeof(float)), *y = xcalloc(ys, sizeof(float)),
        *m = xcalloc((size_t)B * ANNO_UNET_WINDOW, sizeof(float));
  uint64_t seed = 11;
  for (size_t i = 0; i < xs; i++)
    x[i] = anno_random_symmetric(&seed);
  /* window 1 has no minus-strand classes: their gradient must still flow */
  for (int b = 0; b < B; b++)
    for (int i = 0; i < ANNO_UNET_WINDOW; i++) {
      int c = (i / 300) % ANNO_CLASSES;
      if (b == 1 && c >= ANNO_MINUS_START_SINGLE_EXON)
        c = ANNO_INTERGENIC;
      y[(size_t)b * ANNO_CLASSES * ANNO_UNET_WINDOW +
        (size_t)c * ANNO_UNET_WINDOW + i] = 1.0f;
      m[b * ANNO_UNET_WINDOW + i] = i < ANNO_UNET_WINDOW - 97;
    }
  UnetConfig config = unet_default_config();
  Unet *u = unet_create(&config);
  float *gw[UNET_LAYER_COUNT], *gb[UNET_LAYER_COUNT];
  for (int l = 0; l < u->nlayer; l++) {
    gw[l] = xcalloc((size_t)u->layer[l].in * u->layer[l].out * u->layer[l].k,
                    sizeof(float));
    gb[l] = xcalloc((size_t)u->layer[l].out, sizeof(float));
  }
  double loss = 0, plus = 0, minus = 0;
  CHECK(batch_gradient(u, x, y, m, B, gw, gb, &loss) == 0);
  /* encoder input, skip source, decoder concat and output layers */
  static const int layers[] = {0, 3, 6, 8, 10, 14};
  double worst = 0;
  int tested = 0;
  for (int t = 0; t < 6; t++) {
    Layer *l = &u->layer[layers[t]];
    size_t count = (size_t)l->in * l->out * l->k, best = 0;
    for (size_t i = 1; i < count; i++) /* largest gradient: well above noise */
      if (fabsf(gw[layers[t]][i]) > fabsf(gw[layers[t]][best]))
        best = i;
    float w = l->w[best], h = 1e-2f * fmaxf(fabsf(w), 1e-2f);
    float *scratch[UNET_LAYER_COUNT], *sb[UNET_LAYER_COUNT];
    for (int k = 0; k < u->nlayer; k++) {
      scratch[k] = xcalloc((size_t)u->layer[k].in * u->layer[k].out *
                               u->layer[k].k,
                           sizeof(float));
      sb[k] = xcalloc((size_t)u->layer[k].out, sizeof(float));
    }
    l->w[best] = w + h;
    batch_gradient(u, x, y, m, B, scratch, sb, &plus);
    l->w[best] = w - h;
    batch_gradient(u, x, y, m, B, scratch, sb, &minus);
    l->w[best] = w;
    for (int k = 0; k < u->nlayer; k++) {
      free(scratch[k]);
      free(sb[k]);
    }
    double numeric = (plus - minus) / (2.0 * (double)h);
    double analytic = gw[layers[t]][best];
    if (fabs(numeric) > 1e-7) {
      tested++;
      worst = fmax(worst, fabs(numeric - analytic) / fabs(numeric));
    }
  }
  CHECK(tested == 6);
  CHECK(worst < 5e-2);
  if (worst >= 5e-2)
    fprintf(stderr, "gradient relative error %g\n", worst);
  for (int l = 0; l < u->nlayer; l++) {
    free(gw[l]);
    free(gb[l]);
  }
  unet_destroy(u);
  free(x);
  free(y);
  free(m);
}

/** @brief Longest ORF, stop codons and CDS phase. */
static void test_orf(void) {
  int64_t start = INT64_MAX, end = 0;
  /* frame 2: ATG AAA TAG -> [2,11) */
  CHECK(longest_orf("CCATGAAATAGCC", 13, &start, &end) && start == 2 &&
        end == 11);
  start = INT64_MAX;
  CHECK(longest_orf("ATGCCCCCCTGACC", 14, &start, &end) && start == 0 &&
        end == 12);
  start = INT64_MAX;
  CHECK(!longest_orf("ATGCCCCCC", 9, &start, &end)); /* no stop -> no ORF */
  CHECK(!longest_orf("CCCTAA", 6, &start, &end));    /* no ATG -> no ORF */
  CHECK(is_stop("TAA") && is_stop("TAG") && is_stop("TGA") && !is_stop("TGG"));
  /* exon run without ORF is not written */
  AnnoSequence s = {"x", "CCCCCCCCCCCCCCCCCCCC", 20, NULL};
  uint8_t path[20] = {0};
  for (int i = 5; i < 15; i++)
    path[i] = ANNO_PLUS_START_SINGLE_EXON;
  FILE *file = tmpfile();
  unsigned long genes = 0;
  CHECK(anno_write_gff(file, &s, path, &genes) == 0 && genes == 0);
  CHECK(ftell(file) == 0);
  fclose(file);
}

/** @brief Production softmax stays finite and normalised for extreme logits. */
static void test_softmax(void) {
  UnetConfig config = unet_default_config();
  Unet *u = unet_create(&config);
  Layer *output = &u->layer[u->nlayer - 1];
  memset(output->w, 0,
         (size_t)output->in * output->out * output->k * sizeof(float));
  float *input = xcalloc((size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW,
                         sizeof(float));
  float p[ANNO_CLASSES];
  for (int j = 0; j < ANNO_CLASSES; j++)
    output->b[j] = 1000.0f - (float)j;
  CHECK(unet_predict_channels(u, input, 1, p) == 0);
  double sum = 0.0, z = 0.0;
  for (int j = 0; j < ANNO_CLASSES; j++) {
    CHECK(isfinite(p[j]));
    sum += p[j];
    z += exp(-(double)j);
  }
  CHECK(fabs(sum - 1.0) < 1e-6);
  CHECK(fabs(p[0] - 1.0 / z) < 1e-6 && fabs(p[1] - exp(-1.0) / z) < 1e-6);
  for (int j = 0; j < ANNO_CLASSES; j++)
    output->b[j] = -1e30f;
  CHECK(unet_predict_channels(u, input, 1, p) == 0);
  CHECK(fabs(p[3] - 1.0 / ANNO_CLASSES) < 1e-6);
  free(input);
  unet_destroy(u);
}

/** @brief Morlet CWT: strand conjugate identity and block independence. */
static void test_cwt(void) {
  enum { N = 3000 };
  char *bases = xcalloc(N + 1, 1), *reverse = xcalloc(N + 1, 1);
  uint64_t seed = 1;
  for (int i = 0; i < N; i++)
    bases[i] = "ACGT"[random_below(&seed, 4)];
  for (int i = 0; i < N; i++)
    reverse[i] = complement(bases[N - 1 - i]);
  float *a = xcalloc((size_t)N * ANNO_CWT_CHANNELS * ANNO_COMPLEX_COMPONENTS,
                     sizeof(float));
  float *b = xcalloc((size_t)N * ANNO_CWT_CHANNELS * ANNO_COMPLEX_COMPONENTS,
                     sizeof(float));
  anno_cwt(bases, N, 0, N, a);
  anno_cwt(reverse, N, 0, N, b);
  double worst = 0.0;
  for (int i = 0; i < N; i++)
    for (int s = 0; s < ANNO_CWT_CHANNELS; s++) {
      const float *x =
          a + ((size_t)i * ANNO_CWT_CHANNELS + s) * ANNO_COMPLEX_COMPONENTS;
      const float *y = b + ((size_t)(N - 1 - i) * ANNO_CWT_CHANNELS + s) *
                               ANNO_COMPLEX_COMPONENTS;
      worst = fmax(worst, fabs(x[0] - y[0]) + fabs(x[1] + y[1]));
    }
  CHECK(worst < 1e-4);
  /* power via parallel 1024-bp blocks equals the one-shot transform */
  AnnoSequence sequence = {"r", bases, N, NULL};
  float *power = test_sequence_power(&sequence);
  worst = 0.0;
  for (int i = 0; i < N * ANNO_CWT_CHANNELS; i++)
    worst =
        fmax(worst, fabs(power[i] - (a[ANNO_COMPLEX_COMPONENTS * i] *
                                         a[ANNO_COMPLEX_COMPONENTS * i] +
                                     a[ANNO_COMPLEX_COMPONENTS * i + 1] *
                                         a[ANNO_COMPLEX_COMPONENTS * i + 1])));
  CHECK(worst < 1e-5);
  /* a constant (poly-A) signal far from the ends has ~zero response */
  memset(bases, 'A', N);
  anno_cwt(bases, N, N / 2, N / 2 + 1, a);
  for (int s = 0; s < ANNO_CWT_CHANNELS; s++)
    CHECK(fabs(a[ANNO_COMPLEX_COMPONENTS * s]) < 1e-3 &&
          fabs(a[ANNO_COMPLEX_COMPONENTS * s + 1]) < 1e-3);
  free(power);
  free(a);
  free(b);
  free(bases);
  free(reverse);
}

/** @brief Short-window predictions ignore values in the padded tail. */
static void test_padding_mask(void) {
  enum { N = 137 };
  size_t input_size = (size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW;
  float *zero_padded = xcalloc(input_size, sizeof(float));
  float *noisy_padded = xcalloc(input_size, sizeof(float));
  float a[N * ANNO_CLASSES], b[N * ANNO_CLASSES];
  uint64_t seed = 19;
  for (int ch = 0; ch < ANNO_INPUT_CHANNELS; ++ch)
    for (int p = 0; p < ANNO_UNET_WINDOW; ++p) {
      size_t index = (size_t)ch * ANNO_UNET_WINDOW + p;
      if (p < N) {
        float value = anno_random_symmetric(&seed);
        zero_padded[index] = noisy_padded[index] = value;
      } else {
        noisy_padded[index] = anno_random_symmetric(&seed);
      }
    }
  Unet *unet = unet_create(&(UnetConfig){
      .depth = 4, .kernel = 9, .pool = 2, .channels = {16, 32, 64, 128}});
  CHECK(unet_predict_channels(unet, zero_padded, N, a) == 0);
  CHECK(unet_predict_channels(unet, noisy_padded, N, b) == 0);
  CHECK(!memcmp(a, b, sizeof(a)));
  unet_destroy(unet);
  free(zero_padded);
  free(noisy_padded);
}

/** @brief model round trip and identical predictions. */
static void test_model_file(void) {
  UnetConfig config = {
      .depth = 4, .kernel = 9, .pool = 2, .channels = {16, 32, 64, 128}};
  Unet *unet = unet_create(&config), *loaded = NULL;
  const char *path = "test-model.ann";
  CHECK(unet_save(path, unet) == 0);
  CHECK(unet_load(path, &loaded) == 0 && loaded);
  float *input = xcalloc((size_t)ANNO_INPUT_CHANNELS * ANNO_UNET_WINDOW,
                         sizeof(float));
  float p[100 * ANNO_CLASSES], q[100 * ANNO_CLASSES];
  uint64_t seed = 3;
  for (int c = 0; c < ANNO_INPUT_CHANNELS; c++)
    for (int i = 0; i < 100; i++)
      input[(size_t)c * ANNO_UNET_WINDOW + i] = anno_random_unit(&seed);
  CHECK(unet_predict_channels(unet, input, 100, p) == 0 &&
        unet_predict_channels(loaded, input, 100, q) == 0);
  CHECK(!memcmp(p, q, sizeof(p)));
  for (int i = 0; i < 100; i++) {
    double sum = 0.0;
    for (int j = 0; j < ANNO_CLASSES; j++)
      sum += p[i * ANNO_CLASSES + j];
    CHECK(fabs(sum - 1.0) < 1e-5);
  }
  free(input);
  remove(path);
  unet_destroy(unet);
  unet_destroy(loaded);
}

/** @brief Whether anno_read_labels() rejects an annotation (runs in a child).
 */
static int labels_rejected(const char *fasta, const char *line) {
  write_text("test-bad.gff3", line);
  fflush(NULL);
  pid_t child = fork();
  if (!child) {
    if (!freopen("/dev/null", "w", stderr))
      _exit(2);
    int count = 0;
    AnnoSequence *sequences = anno_read_fasta(fasta, &count);
    anno_read_labels("test-bad.gff3", sequences, count);
    _exit(0);
  }
  int status = 0;
  waitpid(child, &status, 0);
  remove("test-bad.gff3");
  return WIFEXITED(status) && WEXITSTATUS(status) != 0;
}

/** @brief Invalid annotations are errors. */
static void test_label_errors(const char *fasta) {
  CHECK(labels_rejected(fasta, "t\ta\texon\t10\t2001\t.\t+\t.\tParent=x\n"));
  CHECK(labels_rejected(fasta, "t\ta\texon\t20\t10\t.\t+\t.\tParent=x\n"));
  CHECK(labels_rejected(fasta, "t\ta\texon\t10\t20\t.\t.\t.\tParent=x\n"));
  CHECK(labels_rejected(fasta, "u\ta\texon\t10\t20\t.\t+\t.\tParent=x\n"));
  CHECK(!labels_rejected(fasta, "t\ta\texon\t10\t20\t.\t+\t.\tParent=x\n"));
}

/** @brief Tiny train -> save -> load -> predict pipeline. */
static void test_pipeline(AnnoSequence *s) {
  Unet *unet = train_unet(s, 1, 3, 2);
  CHECK(unet_save("test-pipeline.ann", unet) == 0);
  Unet *loaded = NULL;
  CHECK(unet_load("test-pipeline.ann", &loaded) == 0);
  FILE *file = tmpfile();
  predict_genome(loaded, s, 1, file, NULL);
  char *text = read_text(file);
  /* every line is a 9-column record of a known type */
  int ok = 1;
  for (char *line = strtok(text, "\n"); line; line = strtok(NULL, "\n")) {
    int tabs = 0;
    for (char *c = line; *c; c++)
      tabs += *c == '\t';
    ok &= tabs == 8;
  }
  CHECK(ok);
  free(text);
  fclose(file);
  remove("test-pipeline.ann");
  unet_destroy(unet);
  unet_destroy(loaded);
}

int main(void) {
  AnnoSequence *sequences = test_sequences("test.fasta", "test.gff3");
  test_labels(&sequences[0]);
  test_strand_swap();
  test_round_trip(&sequences[0]);
  test_single_exons();
  test_nested_transcripts(&sequences[0]);
  test_grammar();
  test_grammar_paths();
  test_viterbi();
  test_invalid_signals();
  test_minimum_decode();
  test_utr_intron_edges();
  test_orf();
  test_softmax();
  test_cwt();
  test_model_file();
  test_padding_mask();
  test_gradient();
  test_label_errors("test.fasta");
  test_pipeline(sequences);
  anno_sequences_free(sequences, 1);
  remove("test.fasta");
  remove("test.gff3");
  printf("%d/%d checks passed\n", checks - failures, checks);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
