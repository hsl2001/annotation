
import argparse
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from itertools import product
from pathlib import Path
from threading import Timer

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.accuracy import (
    Counts,
    calculate,
    chain_keys,
    coding_sites,
    intron_keys,
    junction_sites,
    read_gff3,
)


def run(*command, data=None, timeout=60, **kwargs):
    return subprocess.check_output(
        list(map(str, command)), cwd=ROOT, input=data, text=True, timeout=timeout, **kwargs,
    )


def load(path, text=None):
    if text is not None:
        path.write_text(text, encoding="utf-8")
    return read_gff3(path)[0]


def models(path, forms, kind="CDS", utr=False):
    rows = []
    for tid, gid, strand, spans in forms:
        parent = f"Parent={tid}"
        features = [(kind, a, b, parent) for a, b in spans]
        if gid:
            features.append(("mRNA", spans[0][0], spans[-1][1], f"ID={tid};Parent={gid}"))
        if utr:
            features += [("exon", a - 3, b + 3, parent) for a, b in spans] + [
                ("five_prime_UTR", 1, 3, parent), ("three_prime_UTR", 90, 99, parent)]
        rows.extend(f"chr\t.\t{k}\t{a}\t{b}\t.\t{strand}\t0\t{attrs}\n"
                    for k, a, b, attrs in features)
    return load(path, "".join(rows))


def regression_checks(directory):
    path = directory / "models.gff3"
    forms = [("a", "g", "+", ((10, 15), (30, 35))),
             ("b", "g", "+", ((10, 18), (30, 35))),
             ("c", "h", "-", ((50, 55), (70, 75)))]
    ref = models(path, forms)
    perfect = calculate(ref, ref)
    utr = models(path, forms, utr=True)
    assert calculate(ref, utr) == calculate(utr, ref) == perfect
    assert calculate(ref, models(path, forms, kind="exon")) == perfect
    duplicates = models(path, forms + [("duplicate-a", "different-gene", "+", forms[0][3]),
                                       ("duplicate-b", "g", "+", forms[1][3])])
    assert calculate(ref, duplicates) == calculate(duplicates, ref) == perfect
    assert perfect["transcript-CDS-chain"].reference == 3
    assert perfect["alternative-isoform"].reference == 2
    chains = chain_keys(ref)
    assert coding_sites(chains, "start") == {("chr", "+", 10), ("chr", "-", 75)}
    assert coding_sites(chains, "stop") == {("chr", "+", 35), ("chr", "-", 50)}
    introns = intron_keys(chains)
    assert introns == {("chr", "+", 16, 29), ("chr", "+", 19, 29), ("chr", "-", 56, 69)}
    assert junction_sites(introns, "donor") == {("chr", "+", 16), ("chr", "+", 19), ("chr", "-", 69)}
    assert junction_sites(introns, "acceptor") == {("chr", "+", 29), ("chr", "-", 56)}
    wrong = models(path, [forms[0], ("wrong", "g", "+", ((10, 21), (30, 35))), forms[2]])
    assert calculate(ref, wrong)["alternative-isoform"] == Counts(1, 1, 1)
    assert calculate(ref[:1], wrong[:2])["transcript-CDS-chain"].fp == 1
    extra = models(path, forms + [("extra", "h", "-", ((50, 55), (70, 78)))])
    assert calculate(ref, extra)["alternative-isoform"].fp == 1
    print("Evaluator regressions passed", flush=True)


def predictions(path, text=None):
    transcripts = load(path, text)
    coding, genes = defaultdict(list), defaultdict(set)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            assert "utr" not in fields[2].lower(), (path, fields[2])
            if fields[2] == "CDS":
                attrs = dict(p.split("=", 1) for p in fields[8].split(";") if "=" in p)
                for parent in attrs["Parent"].split(","):
                    coding[parent].append(tuple(map(int, (fields[3], fields[4], fields[7]))))
    assert len(coding) == len(transcripts), path
    for t in transcripts:
        length = 0
        for a, b, phase in sorted(coding[t.tid], reverse=t.strand == "-"):
            assert phase == (-length) % 3, (path, t.tid, phase, length)
            length += b - a + 1
        assert length >= 6 and length % 3 == 0, (path, t.tid, length)
        assert t.spans["CDS"] == t.spans["exon"]
        genes[t.seqid, t.strand, t.gid].add(tuple(t.spans["CDS"]))
    return transcripts, sum(len(forms) - 1 for forms in genes.values())


def decoder_checks(directory):
    fixture = r'''
#include "src/predict.c"
#include <assert.h>

int main(int argc, char **argv) {
  char dna[20001], genomic[20001];
  assert(argc >= 2 && scanf("%20000s", dna) == 1);
  Sequence s = {.name = "chr", .bases = dna, .length = (int64_t)strlen(dna)};
  Model model = {0};
  if (argc == 2) {
    Sample sample = {.type = START};
    double expected = model.weights[START][SITE_FEATURES - 1] = 0.25;
    for (int k = 0; k < SITE_WIDTH; ++k) {
      model.weights[START][k * 5] = k + 1;
      model.weights[START][k * 5 + 4] = 0.5;
      sample.power[k] = 1;
      expected += k + 1.5;
    }
    assert(SITE_FEATURES == SITE_WIDTH * 5 + 1);
    assert(JUNCTION_FEATURES == SITE_WIDTH * 9 + 3 && MODEL_FEATURES == JUNCTION_FEATURES);
    assert(fabs(site_logit(&model, &sample) - expected) < 1e-9);
    model = (Model){0};
    for (int bias = -1; bias <= 0; ++bias) {
      for (int type = 0; type < SITES; ++type)
        model.weights[type][SITE_FEATURES - 1] = bias;
      Candidates found = detect_sites(&model, &s);
      int n = 0;
      for (int64_t p = 0; !bias && p <= s.length; ++p)
        for (int type = 0; type < SITES; ++type)
          if (site_candidate(s.bases, s.length, p, type)) {
            assert(n < found.count && found.points[n].position == p);
            assert((found.points[n].kind & 3) == (unsigned)type && found.points[n].probability == 0.5f);
            ++n;
          }
      assert(found.count == n);
      free(found.points); free(found.openings);
    }
    s.length = 180;
    read_gff3(argv[1], &s, 1);
    Interval pairs[] = {{30, 60}, {90, 120}, {30, 120}, {45, 150},
                        {120, 150}, {90, 150}, {45, 120}};
    assert(s.introns[0].count == 4 && s.introns[1].count == 1);
    for (int i = 0; i < 7; ++i) {
      IntervalList list = s.introns[i == 4];
      assert((bsearch(&pairs[i], list.items, list.count, sizeof(Interval), compare_interval) != NULL) == (i < 5));
    }
    free(s.labels); free(s.introns[0].items); free(s.introns[1].items);
    return 0;
  }
  assert(argc == 5);
  int reverse = atoi(argv[1]);
  for (int64_t p = 0; p < s.length; ++p)
    genomic[p] = reverse ? complement(dna[s.length - 1 - p]) : dna[p];
  genomic[s.length] = 0;
  Sequence genome = s;
  genome.bases = genomic;
  model.max_intron = atoll(argv[3]);
  model.weights[JUNCTION][JUNCTION_FEATURES - 1] = atof(argv[4]);

  Point points[128];
  Candidates c = {.points = points};
  for (int type = 0; type < SITES; ++type) {
    int count;
    assert(scanf("%d", &count) == 1 && count >= 0 && c.count + count <= 128);
    for (int e = 0; e < count; ++e) {
      int p;
      assert(scanf("%d", &p) == 1);
      points[c.count++] = (Point){p, atof(argv[2]), (unsigned)type};
    }
  }
  unsigned long genes = 0;
  prepare_openings(&c, &s);
  assemble(&model, &s, &genome, reverse ? '-' : '+', &c, stdout, &genes);
  free(c.openings);
  return 0;
}
'''
    normal_donors = (1, 2, 30, 31, 60, 90, 94, 120)
    normal_acceptors = (31, 32, 60, 62, 90, 120, 123, 126, 150)
    complete_orf = re.compile(r"ATG(?:(?!TAA|TAG|TGA)[ACGT]{3})*(?:TAA|TAG|TGA)")

    def possibilities(begin, maximum):
        yield frozenset(range(begin, len(dna)))
        for donor, acceptor in product(donors, acceptors):
            if begin < donor and 20 <= acceptor - donor <= maximum:
                for tail in possibilities(acceptor, maximum):
                    yield frozenset(range(begin, donor)) | tail


    annotation = directory / "paired.gff3"
    models(annotation, [("A,B", "", "+", ((1, 30), (121, 150))),
                        ("A", "", "+", ((61, 90),)), ("C", "", "+", ((41, 45), (151, 180))),
                        ("D", "", "-", ((1, 30), (61, 90)))])
    libraries = [p for p in (ROOT / "src").glob("*.c") if p.stem not in ("predict", "surf")]
    executable = directory / "decoder-check"
    path = directory / "decoder-output.gff3"
    for flags, threshold in (([], 0.5), (["-DSITE_RADIUS=7", "-DCDS_RECIPROCAL_OVERLAP=0.75"], 0.75)):
        run(os.environ.get("CC", "cc"), "-O1", "-Wall", "-Wextra", "-Werror", "-fopenmp",
            "-I.", "-Isrc", *flags, "-x", "c", "-", *libraries, "-lz", "-lm", "-o", executable, data=fixture)
        run(executable, annotation, data="ATGGTAG" * 2500,
            timeout=30, env={**os.environ, "OMP_NUM_THREADS": "1"})
        for mode, maximum in ((0, 1000), (1, 1000), (2, 1000), (3, 1000), (4, 1000), (2, 50)):
            dna = list(("ATGATG" + "A" * 171 if mode else "ATG" + "A" * 294) + "TAA")
            for p, base in {0: {140: "N"}, 2: {30: "T"}, 3: {80: "N"}, 4: {31: "T", 32: "G"}}.get(mode, {}).items():
                dna[p] = base
            dna = "".join(dna)
            starts = (0, 3) if mode else (0,)
            donors = normal_donors if mode else (*range(3, 121, 3), *range(4, 62, 3))
            acceptors = normal_acceptors if mode else (181, 240)
            possible = {(s, len(dna)): {c for c in possibilities(s, maximum)
                                       if complete_orf.fullmatch("".join(dna[p] for p in sorted(c)))} for s in starts}
            data = dna + " " + " ".join(f"{len(b)} " + " ".join(map(str, b))
                                       for b in (starts, donors, acceptors, (len(dna),)))
            for reverse, gate in product((0, 1), (0, -1)):
                allowed = {pair: keep for pair, forms in possible.items()
                           if (keep := {s for s in forms if not gate or len(s) == pair[1] - pair[0]})}
                output = run(executable, reverse, 0.01, maximum, gate, data=data, timeout=10)
                assert output == run(executable, reverse, 0.99, maximum, gate, data=data, timeout=10)
                groups = defaultdict(list)
                for t in predictions(path, output)[0]:
                    chain = frozenset(len(dna) - p if reverse else p - 1
                                      for a, b in t.spans["CDS"] for p in range(a, b + 1))
                    groups[min(chain), max(chain) + 1].append(chain)
                assert groups.keys() == allowed.keys(), (mode, reverse, gate)
                for pair, actual in groups.items():
                    reference = actual[0]
                    assert len(reference) == max(map(len, allowed[pair]))
                    expected = {s for s in allowed[pair]
                                if len(s & reference) >= threshold * max(len(s), len(reference))}
                    assert set(actual) == expected and len(actual) == len(expected), (mode, reverse, gate, pair)
                    if not mode:
                        assert len(actual) > 4
                        assert any(len(s & reference) == threshold * len(reference) for s in allowed[pair])
                        assert any(len(s) == len(reference) and len(s & reference) < 0.75 * len(reference)
                                   for s in allowed[pair])
    print("Decoder regressions passed", flush=True)


def predict_to_file(command, path):
    with path.open("wb") as output, subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE) as process:
        assert process.stdout is not None
        timer = Timer(60, process.kill)
        timer.daemon = True
        timer.start()
        size = 0
        try:
            while chunk := process.stdout.read(65536):
                size += len(chunk)
                if size > 64 * 1024 * 1024:
                    raise AssertionError(f"{path}: output limit exceeded")
                output.write(chunk)
            if process.wait():
                raise subprocess.CalledProcessError(process.returncode, command)
        finally:
            timer.cancel()
            if process.poll() is None:
                process.kill()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--iterations", type=int, default=100)
    parser.add_argument("-p", "--threads", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--regressions-only", action="store_true")
    args = parser.parse_args()
    if min(args.iterations, args.threads) < 1:
        parser.error("iterations and threads must be positive")
    with tempfile.TemporaryDirectory(prefix="surf-integration-") as directory:
        directory = Path(directory)
        regression_checks(directory)
        decoder_checks(directory)
        if args.regressions_only:
            return
        run("make", "-s")
        command = ["./surf", "-p", str(args.threads), "-m", str(directory / "model.crf")]
        truth = ROOT / "data/truth.nip.chr5"
        started = time.perf_counter()
        run(*command, "-i", args.iterations, f"{truth}.100k.fasta", f"{truth}.100k.gff3", timeout=6000)
        if time.perf_counter() - started > 6000:
            raise AssertionError("Training exceeded 6000 seconds")
        for size, minima in (("100k", (0.85, 0.70)), ("3m", (0.80, 0.65))):
            path = directory / f"{size}.gff3"
            predict_to_file(command + [f"{truth}.{size}.fasta"], path)
            transcripts, additional = predictions(path)
            assert transcripts, path
            result = calculate(load(f"{truth}.{size}.gff3"), transcripts)
            for metric, minimum in zip(("nt", "exon-segment"), minima):
                f1 = result[metric].as_dict()["f1"]
                if f1 is None or f1 < minimum:
                    raise AssertionError((size, metric, f1, minimum))
            assert size != "3m" or additional > 0, "3m lacks alternative coding forms"


if __name__ == "__main__":
    main()
