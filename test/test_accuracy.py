import tempfile
import unittest
from pathlib import Path

from scripts.accuracy import metrics


class AccuracyTest(unittest.TestCase):
    def test_partial_and_unannotated_predictions(self):
        with tempfile.TemporaryDirectory() as directory:
            reference = Path(directory) / "reference.gff3"
            query = Path(directory) / "query.gff3"
            reference.write_text("chr1\t.\texon\t10\t19\t.\t+\t.\tParent=r\n")
            query.write_text(
                "chr1\t.\texon\t5\t14\t.\t+\t.\tParent=q\n"
                "chr2\t.\texon\t1\t10\t.\t+\t.\tParent=other\n"
            )
            result = metrics(reference, query)
            self.assertEqual(result["exon-segment"], (0.0, 0.0, 0.0))
            self.assertEqual(result["nt"], (0.25, 0.5, 1 / 3))

            query.write_text(
                "chr1\t.\texon\t10\t19\t.\t+\t.\tParent=q\n"
                "chr1\t.\texon\t30\t39\t.\t+\t.\tParent=other\n"
            )
            result = metrics(reference, query)
            self.assertEqual(result["exon-segment"], (0.5, 1.0, 2 / 3))
            self.assertEqual(result["nt"], (0.5, 1.0, 2 / 3))


if __name__ == "__main__":
    unittest.main()
