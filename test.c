/**
 * @file test.c
 * @brief Unit and pipeline tests for anno. 
 */
#define ANNO_NO_MAIN
#include "anno.c"

#include <sys/wait.h>
#include <unistd.h>

static int failures = 0, checks = 0;

/** @brief Return a deterministic random integer below the given bound. */
static unsigned random_below(uint64_t *state, unsigned bound) {
  return (unsigned)(anno_random_next(state) % bound);
}

/** @brief Return a deterministic pseudo-random float in [0, 1). */
static float test_random_unit(uint64_t *state) {
  return anno_random_unit(state);
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
  CHECK(y[39] == ANNO_INTERGENIC && y[40] == ANNO_PLUS_5UTR);
  CHECK(y[70] == ANNO_PLUS_INTRON && y[99] == ANNO_PLUS_INTRON);
  CHECK(y[120] == ANNO_PLUS_5UTR && y[121] == ANNO_PLUS_EXON);
  CHECK(y[390] == ANNO_PLUS_EXON && y[391] == ANNO_PLUS_3UTR);
  CHECK(y[419] == ANNO_PLUS_3UTR && y[420] == ANNO_INTERGENIC);
  CHECK(y[1579] == ANNO_INTERGENIC && y[1580] == ANNO_MINUS_3UTR);
  CHECK(y[1609] == ANNO_MINUS_EXON && y[1640] == ANNO_MINUS_INTRON);
  CHECK(y[1879] == ANNO_MINUS_5UTR && y[1959] == ANNO_MINUS_5UTR);
  CHECK(y[1960] == ANNO_INTERGENIC);
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

/** @brief GFF -> labels -> GFF and GFF -> labels -> probabilities -> Viterbi ->
 * GFF. */
static void test_round_trip(const AnnoSequence *s) {
  char *direct = gff_of_path(s, s->labels);
  CHECK(!strcmp(direct, test_gff));
  if (strcmp(direct, test_gff))
    fprintf(stderr, "got:\n%s", direct);
  float *probability = xcalloc((size_t)s->length * ANNO_CLASSES, sizeof(float));
  for (int64_t i = 0; i < s->length; i++)
    for (int j = 0; j < ANNO_CLASSES; j++)
      probability[i * ANNO_CLASSES + j] = j == s->labels[i] ? 0.6f : 0.05f;
  uint8_t *path = xcalloc((size_t)s->length, 1);
  CHECK(anno_decode(s->bases, s->length, probability, path) == 0);
  CHECK(!memcmp(path, s->labels, (size_t)s->length));
  char *decoded = gff_of_path(s, path);
  CHECK(!strcmp(decoded, test_gff));
  free(direct);
  free(decoded);
  free(path);
  free(probability);
}

/** @brief Grammar constraints on a probability matrix. */
static void test_viterbi(void) {
  enum { N = 60 };
  char bases[N + 1];
  memset(bases, 'C', N);
  bases[N] = 0;
  float p[N * ANNO_CLASSES];
  uint8_t path[N];
  /* an "intron" signal: it may only be reached through a >= 6 bp exon */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < ANNO_CLASSES; j++)
      p[i * ANNO_CLASSES + j] =
          j == (i >= 10 && i < 50 ? ANNO_PLUS_INTRON : 0) ? 0.9f : 0.0125f;
  CHECK(anno_decode(bases, N, p, path) == 0);
  CHECK(path[0] == ANNO_INTERGENIC && path[N - 1] == ANNO_INTERGENIC);
  for (int i = 1; i < N; i++)
    CHECK(!(path[i - 1] == ANNO_INTERGENIC && path[i] == ANNO_PLUS_INTRON) &&
          !(path[i - 1] == ANNO_PLUS_INTRON && path[i] == ANNO_INTERGENIC));
  /* a 4-bp exon signal is widened to the 6-bp minimum (gain 4*6.68 > 8 +
   * 2*6.68) */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < ANNO_CLASSES; j++)
      p[i * ANNO_CLASSES + j] =
          j == (i >= 20 && i < 24 ? ANNO_PLUS_EXON : 0) ? 0.99f : 0.00125f;
  CHECK(anno_decode(bases, N, p, path) == 0);
  int exon = 0;
  for (int i = 0; i < N; i++)
    exon += path[i] == ANNO_PLUS_EXON;
  CHECK(exon == ANNO_MIN_EXON_LENGTH);
  /* 5'UTR followed directly by 3'UTR is not a valid transcript */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < ANNO_CLASSES; j++)
      p[i * ANNO_CLASSES + j] = j == (i >= 10 && i < 30   ? ANNO_PLUS_5UTR
                                      : i >= 30 && i < 50 ? ANNO_PLUS_3UTR
                                                          : 0)
                                    ? 0.9f
                                    : 0.0125f;
  CHECK(anno_decode(bases, N, p, path) == 0);
  for (int i = 1; i < N; i++)
    CHECK(!(path[i - 1] == ANNO_PLUS_5UTR && path[i] == ANNO_PLUS_3UTR));
  /* a strong 30-bp exon is decoded exactly (gain 30*log(0.99/0.00125) >
   * penalty) */
  for (int i = 0; i < N; i++)
    for (int j = 0; j < ANNO_CLASSES; j++)
      p[i * ANNO_CLASSES + j] =
          j == (i >= 20 && i < 50 ? ANNO_MINUS_EXON : 0) ? 0.99f : 0.00125f;
  CHECK(anno_decode(bases, N, p, path) == 0);
  for (int i = 0; i < N; i++)
    CHECK(path[i] == (i >= 20 && i < 50 ? ANNO_MINUS_EXON : ANNO_INTERGENIC));
  /* NaN emissions are treated as zero probability, not propagated */
  p[5 * ANNO_CLASSES] = NAN;
  CHECK(anno_decode(bases, N, p, path) == 0);
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
    path[i] = ANNO_PLUS_EXON;
  FILE *file = tmpfile();
  unsigned long genes = 0;
  CHECK(anno_write_gff(file, &s, path, &genes) == 0 && genes == 0);
  CHECK(ftell(file) == 0);
  fclose(file);
}

/** @brief Softmax of extreme logits stays finite and normalised. */
static void softmax_row(const float *logits, float *probability) {
  float top = logits[0], sum = 0.0f;
  for (int j = 1; j < ANNO_CLASSES; j++)
    if (logits[j] > top)
      top = logits[j];
  for (int j = 0; j < ANNO_CLASSES; j++) {
    probability[j] = expf(logits[j] - top);
    sum += probability[j];
  }
  for (int j = 0; j < ANNO_CLASSES; j++)
    probability[j] /= sum;
}

static void test_softmax(void) {
  float logits[ANNO_CLASSES], p[ANNO_CLASSES];
  for (int j = 0; j < ANNO_CLASSES; j++)
    logits[j] = 1000.0f - (float)j;
  softmax_row(logits, p);
  double sum = 0.0, z = 0.0;
  for (int j = 0; j < ANNO_CLASSES; j++) {
    CHECK(isfinite(p[j]));
    sum += p[j];
    z += exp(-(double)j);
  }
  CHECK(fabs(sum - 1.0) < 1e-6);
  CHECK(fabs(p[0] - 1.0 / z) < 1e-6 && fabs(p[1] - exp(-1.0) / z) < 1e-6);
  for (int j = 0; j < ANNO_CLASSES; j++)
    logits[j] = -1e30f;
  softmax_row(logits, p);
  CHECK(fabs(p[3] - 1.0 / ANNO_CLASSES) < 1e-6);
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

/** @brief model round trip and identical predictions. */
static void test_model_file(void) {
  UnetConfig config = {
      .depth = 3, .kernel = 9, .pool = 2, .channels = {16, 32, 64}};
  Unet *unet = unet_create(&config), *loaded = NULL;
  const char *path = "test-model.ann";
  CHECK(unet_save(path, unet) == 0);
  CHECK(unet_load(path, &loaded) == 0 && loaded);
  float input[100 * ANNO_INPUT_CHANNELS], p[100 * ANNO_CLASSES],
      q[100 * ANNO_CLASSES];
  uint64_t seed = 3;
  for (int i = 0; i < 100 * ANNO_INPUT_CHANNELS; i++)
    input[i] = test_random_unit(&seed);
  CHECK(unet_predict(unet, input, 100, p) == 0 &&
        unet_predict(loaded, input, 100, q) == 0);
  CHECK(!memcmp(p, q, sizeof(p)));
  for (int i = 0; i < 100; i++) {
    double sum = 0.0;
    for (int j = 0; j < ANNO_CLASSES; j++)
      sum += p[i * ANNO_CLASSES + j];
    CHECK(fabs(sum - 1.0) < 1e-5);
  }
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
  predict_genome(loaded, s, 1, 2, file, NULL);
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
  test_round_trip(&sequences[0]);
  test_viterbi();
  test_orf();
  test_softmax();
  test_cwt();
  test_model_file();
  test_label_errors("test.fasta");
  test_pipeline(sequences);
  anno_sequences_free(sequences, 1);
  remove("test.fasta");
  remove("test.gff3");
  printf("%d/%d checks passed\n", checks - failures, checks);
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
