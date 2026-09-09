# anno: Complex Morlet CNN Gene Predictor

CPU implementation of masked-nucleotide pretraining, representative-transcript
CDS fine-tuning, and a learned splice-aware HMM. No new external libraries are
required: the build uses the existing genann, kseq, zlib, and C math facilities.
The embedded genann implementation by Lewis Van Winkle is used unchanged.

## Build and Run

```sh
make
./anno pretrain genome.fasta pretrained 3 0.01 1
./anno train genome.fasta annotation.gff3 pretrained coding_model 30 0.001 5 1
./anno predict new_genome.fasta coding_model > new_predictions.gff3
```

Arguments after the output prefix:

| Command | Optional arguments and defaults |
| --- | --- |
| `pretrain` | `epochs=3 learning_rate=0.01 seed=1` |
| `train` | `epochs=3 learning_rate=0.001 em_iterations=5 seed=1` |

`train` requires a pretrained model prefix **and a separate output prefix**.
Prediction requires a fine-tuned model containing both CNN and HMM parameters.
Epochs and EM iterations must be positive; the learning rate must be in `(0, 1]`.
FASTA can be plain text or gzip-compressed; annotations are plain text.
Progress and losses go to stderr. Only GFF3 is written to prediction stdout.

Start with a small, separately prepared genomic region. The 4096-base direct
wavelet transform and CPU CNN make chromosome-scale training expensive.
There is no implicit subsampling or sequence truncation.

## Signal and CNN

- Mapping: `A=1`, `C=i`, `G=-i`, `T=-1`; lowercase is supported.
  Ambiguous bases contribute zero to the complex signal.
- Wavelet sizes mean **finite kernel widths in bases**, exactly
  `8, 64, 512, 4096`, not Morlet periods or CNN tile sizes.
- For width `w`, samples use `u=(tap-w/2+0.5)/(w/8)` and the corrected
  complex Morlet `exp(-u*u/2) * (exp(6*i*u) - exp(-18))`.
  The discrete mean is removed and energy normalized to one; CWT correlates
  the signal with the conjugated kernel. Values outside the contig are zero.
  Even-width kernels use offsets `-w/2 ... w/2-1`.
- The CWT is a position-by-four-scales complex matrix. Both real and imaginary
  components are retained as two input channels; magnitude-only input is not used.
- Eight actual shared 3-by-3 convolution layers span position and scale, with
  eight output channels, tanh activations, and residual connections after the
  first layer. Position dilations are `1,2,4,8,16,1,2,4`.
- Each kernel is a zero-hidden-layer linear genann mapping reused across all
  positions and scales. Cross-entropy gradients are accumulated across the
  tile before any weight changes. Backpropagation is explicit because ordinary
  `genann_train` does not support convolution weight sharing or accumulated
  softmax cross-entropy gradients. Pretraining uses SGD; fine-tuning uses Adam
  with beta1=0.9, beta2=0.999, epsilon=1e-8 and bias correction. Both clip the
  global gradient norm at 5.
- Output tiles contain 128 positions with a 38-base CNN halo on each side.
  CWT reads its larger context directly from the full contig, including during
  tiled prediction. There is no pooling or loss of single-base resolution.

## Training Labels

Pretraining masks 15% of valid bases independently on each oriented contig
before computing **any** CWT coefficients. Only masked A/C/G/T positions receive
four-class reconstruction loss. The original base cannot leak into a neighboring
wavelet coefficient. Masks change each epoch; tiles are shuffled. Both forward
and reverse-complement strands are mixed in the same tile order.

Fine-tuning accumulates eight tiles per update, weighting each tile's mean
gradient by its actual target count, including a shorter final tile. The final
partial batch is flushed. Its learning rate decreases exponentially from the
specified rate to one tenth of that rate over the requested epochs. Adam moments
are reset for each `train` invocation; a stage-2 checkpoint can be supplied to
continue CNN fine-tuning, but the HMM is initialized and fitted again. The
four-class pretraining head is not updated by fine-tuning. No class weighting
or positive oversampling is applied.

For datasets of at most 2,000,000 bases, unmasked fine-tuning CWT coefficients
are cached for both strands (128 bytes per genomic base). Larger datasets use
the tiled transform. Pretraining never reuses an unmasked cache. The training
log reports the actual positive count divided by **twice** the genomic length
and its binary entropy baseline. For the current 100 kb subset this is
28,656 / 200,000 = 0.14328, with baseline CE 0.410873, not 0.599043.

Fine-tuning updates the full backbone and a separate two-class head. For each
gene, it selects the transcript with the largest **sum of CDS segment lengths**,
not the longest single segment, genomic span, or exon length. Duplicate identical
segments count once. Ties are broken by transcript ID. Only that transcript's
CDS bases are positive on its strand; every other position is negative, including
UTRs and CDS belonging exclusively to unselected isoforms.

The parser supports standard GFF3 `ID`/`Parent` hierarchies, multiple CDS parents,
and GTF-style `gene_id`/`transcript_id` attributes, including the Liftoff-style
files provided in this workspace. Transcript rows may follow their CDS rows.
Liftoff `extra_copy_number` distinguishes separate gene copies with reused IDs.
FASTA identifiers must match annotation sequence IDs. Unresolvable parents,
invalid coordinates, inconsistent specified CDS phases, and overlapping CDS
segments within one transcript are rejected rather than silently relabeled.
Missing phases are inferred in transcript order from the first CDS phase
(zero if unspecified). Conflicting overlapping genes on the same strand retain
binary CDS labels but their ambiguous structural positions are excluded from
HMM supervised initialization. Reverse-strand labels use transcript-oriented
reverse-complement coordinates throughout training.

## Learned Splice HMM

There are 19 hidden states:

- One noncoding/background state.
- Three coding states for consecutive positions in a codon.
- For each next coding phase, five intron states: first donor base, second donor
  base, intron interior, first acceptor base, second acceptor base.

Allowed transitions preserve coding phase across each intron. An intron must
contain at least five bases and connect two CDS segments. Annotation gaps shorter
than five bases are not used as supervised splice examples. Gaps are inferred
only between CDS segments of the selected transcript, never from a gene-wide
mask. Donor/acceptor nucleotide distributions are learned from sequence, rather
than hard-coded GT/AG rules. A warning is emitted if no splice examples exist.

The observation is the nucleotide and one of eight bins of the CNN CDS
probability. Each state has normalized nucleotide and score-bin categorical
emissions. This uses a proper discrete HMM observation model, rather than
treating discriminative CNN posterior probabilities as generative emissions.
Missing nucleotides are marginalized. Supervised counts initialize the state
identities, transition probabilities, and emissions; then scaled forward-backward
**partially supervised Baum-Welch EM** learns initial, transition, and emission
distributions while constraining known annotation positions to their structural
class: background, CDS, donor base 1/2, intron interior, or acceptor base 1/2.
Codon phase remains latent within compatible classes. Unknown/conflicting
annotation positions and observations without labels remain unconstrained.
This prevents unsupervised EM from reassigning coding states to arbitrary
sequence patterns and destroying annotated splice structure. The topology is
fixed, not its probabilities. Tests compare constrained and unconstrained
forward-backward probabilities against exhaustive path enumeration.
Separate contigs and strands do not acquire artificial cross-boundary transitions.
Viterbi uses log-space scores and the **saved learned HMM**, with no transitions
directly between background and intron states. Viterbi never applies annotation
constraints; prediction requires only FASTA and the saved model. No annotation,
coordinate lookup, or reference transcript is stored in the checkpoint.

## Output and Limits

GFF3 contains `gene -> mRNA -> exon/CDS`, with `gene_biotype=protein_coding`,
unique IDs, 1-based inclusive coordinates, both strands, and transcript-oriented
CDS phases. Each decoded intron connects CDS segments into one predicted gene.
The `exon` features cover **coding portions only**; UTR-inclusive exon boundaries
cannot be learned from a CDS-only binary target.

These are protein-coding **candidates**, not verified complete ORFs. Start/stop
codons and an ORF without internal stops are not hard constraints. Partial coding
starts/ends are allowed, and `prediction_status=unvalidated` is recorded on genes.
The model does not reconstruct all alternative isoforms. EM can alter the initial
state distributions; independent biological accuracy evaluation is necessary.

Training loads FASTA records and labels in memory. CWT/CNN tensors are tiled, but
HMM EM stores approximately `160 * longest_contig_length` bytes for the forward
table and scaling factors, in addition to FASTA, labels, and CNN probabilities.
That table alone is about 2.08 GB for a 13 Mb contig. Viterbi stores approximately
20 bytes per position for path/backpointers, plus CNN probabilities and FASTA.
Direct CWT work is proportional to `N * (8+64+512+4096)` per strand, with extra halo
work. No FFT, GPU, multithreading, or external training framework is introduced.

Model files are versioned text checkpoints named `<prefix>.bin`. They contain
architecture and wavelet metadata, both CNN heads, weights, training counters,
seed, and (after fine-tuning) the HMM. Writes use a temporary file and rename.
Old MLP-Mixer checkpoints are incompatible and are explicitly rejected.
The architecture and checkpoint layout are unchanged by the optimizer and EM
fixes. Existing Morlet CNN checkpoints load, but old HMM parameters only improve
after running `train` again.

## Overfit Evaluation

```sh
make test_overfit
./test_overfit subset_genome.fasta subset_annotation.gff3 subset_100k_fine
./test_overfit --check subset_genome.fasta subset_annotation.gff3 MODEL_PREFIX representative_cds.gff3
./anno predict subset_genome.fasta MODEL_PREFIX > predictions.gff3
./gffcompare -r representative_cds.gff3 -o cds_evaluation predictions.gff3
```

`test_overfit` uses the production annotation parser and evaluates the final
frozen CNN and saved HMM against their actual strand-specific training labels.
It also reports seeded and refitted HMMs for diagnosis, but only the **saved HMM**
is used for the gate and exported prediction. Every decode is passed an unlabeled
observation. The integration test checks that its exported prediction is byte
identical to the normal `predict` command.

`--check` requires CDS base F1 >= 0.99, exact CDS interval F1 >= 0.80, and exact
intron interval F1 >= 0.90; otherwise it exits nonzero. Exact interval matching
requires both boundaries and the strand to match. Optional output arguments
export a coding-only representative reference and a prediction GFF, respectively.
The representative reference is derived from structural labels, so it is intended
for unambiguous loci; overlapping same-strand annotation conflicts are not a
general-purpose transcript reference export.

Comparing coding-only predictions directly against all GFF exons mixes UTRs,
noncoding transcripts and alternative isoforms into the denominator. Keep that
full-annotation evaluation as a separate task, not the CDS overfit criterion.
Neither same-region performance nor an overfit gate measures generalization;
held-out genomic regions still need independent evaluation.

The current 100 kb run passes all three gates: saved-HMM CDS base F1 0.999110,
exact CDS interval F1 0.869198, and exact intron interval F1 0.973545. The original
checkpoint had CDS base F1 0.724415 and zero exact interval matches. This includes
210 additional fine-tuning epochs, not an equal-budget comparison. See the
[100 kb report](overfit_100k/REPORT.md) for preserved checkpoints, strict interval
counts, gffcompare results, training provenance, replay commands and limitations.

## Verification

```sh
make test
./test_data chr1_13m.fasta chr1_13m.gff3
make -B test CFLAGS='-O1 -g -std=c11 -Wall -Wextra -fsanitize=address,undefined -fno-omit-frame-pointer'
make -B all
```

Tests use C and Python's standard library only. No data or tools are downloaded.
Temporary integration fixtures do not overwrite existing FASTA, GFF, or models.
Coverage includes complex mapping and impulse response, masked-input invariance,
finite-difference CNN gradients, loss decrease, full/tiled equivalence, both-strand
transcript selection and phases, Liftoff copies, learned splice motifs, brute-force
HMM likelihood/Viterbi comparison, EM likelihood improvement, checkpoint round-trip,
the complete CLI workflow, GFF hierarchy, and invalid-input rejection.

Passing these tests establishes implementation consistency, not predictive
accuracy. Use nonoverlapping chromosomes or genomic partitions for pretraining,
fine-tuning, validation, and a held-out final test; avoid overlapping CWT contexts
between partitions. No full-genome training or held-out accuracy measurement has
been performed as part of these tests.