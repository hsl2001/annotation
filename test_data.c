#define _POSIX_C_SOURCE 200809L
#include "anno_data.h"
#include "anno_hmm.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc == 3) {
    Dataset *input = dataset_read(argv[1]);
    FILE *annotation = fopen(argv[2], "r");
    assert(annotation);
    annotation_read(input, annotation);
    fclose(annotation);
    dataset_free(input);
    puts("annotation validation passed");
    return 0;
  }
  Dataset *dataset = anno_alloc(1, sizeof(*dataset));
  dataset->count = 1;
  dataset->sequences = anno_alloc(1, sizeof(Sequence));
  dataset->sequences[0].name = strdup("chr");
  dataset->sequences[0].bases = strdup("ACGT");
  dataset->sequences[0].length = 100;
  FILE *stream = tmpfile();
  assert(stream);
  fputs("chr\tx\tCDS\t1\t9\t.\t+\t0\tParent=short\n"
      "chr\tx\tCDS\t1\t5\t.\t+\t0\tParent=long\n"
      "chr\tx\tCDS\t15\t21\t.\t+\t1\tParent=long\n"
      "chr\tx\tCDS\t15\t21\t.\t+\t1\tParent=long\n"
      "chr\tx\tmRNA\t1\t21\t.\t+\t.\tID=long;Parent=gene1\n"
        "chr\tx\tmRNA\t1\t9\t.\t+\t.\tID=short;Parent=gene1\n"
      "chr\tx\tmRNA\t35\t40\t.\t+\t.\tID=long;Parent=gene1;extra_copy_number=1\n"
      "chr\tx\tCDS\t35\t40\t.\t+\t0\tParent=long;extra_copy_number=1\n"
      "chr\tx\tCDS\t70\t76\t.\t-\t1\tgene_id \"gene2\"; transcript_id \"minus\";\n"
      "chr\tx\tCDS\t90\t94\t.\t-\t0\tgene_id \"gene2\"; transcript_id \"minus\";\n", stream);
  rewind(stream);
  annotation_read(dataset, stream);
  fclose(stream);
  assert(dataset->selected_transcripts == 3);
  Sequence *sequence = dataset->sequences;
  assert(sequence->coding[0][34] == 1 && sequence->coding[0][39] == 1);
  assert(sequence->coding[0][4] == 1 && sequence->coding[0][5] == 0);
  assert(sequence->coding[0][14] == 1 && sequence->coding[0][20] == 1);
  assert(sequence->states[0][5] == HMM_DONOR(2));
  assert(sequence->states[0][13] == HMM_DONOR(2) + 4);
  assert(sequence->states[0][14] == HMM_C0 + 2);
  assert(sequence->coding[1][6] == 1 && sequence->coding[1][10] == 1);
  assert(sequence->coding[1][11] == 0 && sequence->coding[1][24] == 1);
  assert(sequence->states[1][11] == HMM_DONOR(2));
  assert(sequence->states[1][23] == HMM_DONOR(2) + 4);
  assert(sequence->states[1][24] == HMM_C0 + 2);
  dataset_free(dataset);
  puts("longest-transcript, GFF3/GTF, splice and reverse-strand tests passed");
  return 0;
}