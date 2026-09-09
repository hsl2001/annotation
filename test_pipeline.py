import pathlib
import subprocess
import tempfile


EXECUTABLE = pathlib.Path(__file__).resolve().parent / "anno"


def run(*arguments, succeeds=True):
    result = subprocess.run(
        [str(EXECUTABLE), *map(str, arguments)], capture_output=True, text=True
    )
    assert (result.returncode == 0) == succeeds, result.stderr
    return result


def main():
    with tempfile.TemporaryDirectory(prefix="anno-test-") as temporary:
        root = pathlib.Path(temporary)
        fasta, gff = root / "tiny.fa", root / "tiny.gff3"
        fasta.write_text(">chr\n" + "ACGT" * 64 + "\n>edge\nANCGT\n")
        gff.write_text(
            "##gff-version 3\n"
            "chr\ttest\tmRNA\t33\t88\t.\t+\t.\tID=plus;Parent=gene1\n"
            "chr\ttest\tCDS\t33\t50\t.\t+\t0\tParent=plus\n"
            "chr\ttest\tCDS\t65\t88\t.\t+\t0\tParent=plus\n"
            "chr\ttest\tmRNA\t33\t60\t.\t+\t.\tID=short;Parent=gene1\n"
            "chr\ttest\tCDS\t33\t60\t.\t+\t0\tParent=short\n"
            "chr\ttest\tmRNA\t161\t216\t.\t-\t.\tID=minus;Parent=gene2\n"
            "chr\ttest\tCDS\t161\t178\t.\t-\t0\tParent=minus\n"
            "chr\ttest\tCDS\t193\t216\t.\t-\t0\tParent=minus\n"
        )
        pretrained, trained = root / "pre", root / "trained"
        assert not run("pretrain", fasta, pretrained, 1, 0.01, 7).stdout
        run("predict", fasta, pretrained, succeeds=False)
        training = run("train", fasta, gff, pretrained, trained, 2, 0.001, 2, 7)
        assert not training.stdout
        assert "84 positive / 522 strand-positions" in training.stderr
        assert "Selected 2 transcripts" in training.stderr
        assert "2 annotated splice junctions" in training.stderr
        prediction = run("predict", fasta, trained).stdout
        assert prediction == run("predict", fasta, trained).stdout
        assert prediction.startswith("##gff-version 3\n")
        evaluator = EXECUTABLE.with_name("test_overfit")
        if evaluator.exists():
            exported = root / "evaluated.gff3"
            evaluation = subprocess.run(
                [str(evaluator), str(fasta), str(gff), str(trained),
                 str(root / "reference.gff3"), str(exported)],
                capture_output=True, text=True,
            )
            assert evaluation.returncode == 0, evaluation.stderr
            assert "positive=84 total=522" in evaluation.stdout
            assert exported.read_text() == prediction
        features = {}
        for line in prediction.splitlines():
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            assert len(fields) == 9
            name, _, kind, start, end, _, strand, phase, attributes = fields
            assert 1 <= int(start) <= int(end) <= {"chr": 256, "edge": 5}[name]
            assert strand in ("+", "-")
            metadata = dict(item.split("=", 1) for item in attributes.split(";"))
            assert metadata["ID"] not in features
            features[metadata["ID"]] = (kind, int(start), int(end))
            if kind == "gene":
                assert metadata["gene_biotype"] == "protein_coding"
            else:
                parent = features[metadata["Parent"]]
                assert parent[1] <= int(start) <= int(end) <= parent[2]
            if kind == "CDS":
                assert phase in ("0", "1", "2")
        broken = root / "broken"
        broken.with_suffix(".bin").write_text("legacy mixer weights\n")
        run("predict", fasta, broken, succeeds=False)
        broken.with_suffix(".bin").write_bytes(trained.with_suffix(".bin").read_bytes()[:300])
        run("predict", fasta, broken, succeeds=False)
        run("pretrain", fasta, root / "unused", 0, succeeds=False)
        run("train", fasta, gff, pretrained, root / "unused", 1, "nan", succeeds=False)
        gff.write_text("chr\ttest\tCDS\t0\t999\t.\t+\t0\tParent=bad\n")
        run("train", fasta, gff, pretrained, root / "bad", 1, succeeds=False)
    print("CLI training, constrained EM, prediction and invalid-input tests passed")


if __name__ == "__main__":
    main()