import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "10_statistical_diversity.py"
PLOT_SCRIPT = PROJECT_ROOT / "scripts" / "11_plot_statistical_diversity.py"


class StatisticalDiversityTests(unittest.TestCase):
    def test_analysis_writes_both_comparison_families_reproducibly(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            results_dir = Path(temporary_directory)
            alignments = {
                "env": (
                    ">HIV__a\nAAAACCCCGGGGTTTT\n"
                    ">HIV__b\nAAAACTCCGGGATTTT\n"
                    ">SIV__a\nGGGGTTTTAAAACCCC\n"
                    ">SIV__b\nGGGGTTTTAAAACCCC\n"
                ),
                "pol": (
                    ">HIV__a\nAAAACCCCGGGGTTTT\n"
                    ">HIV__b\nAAAACCCCGGGGTTTT\n"
                    ">SIV__a\nGGGGTTTTAAAACCCC\n"
                    ">SIV__b\nGGGGTTTTAAAACCCA\n"
                ),
                "gag": (
                    ">HIV__a\nAAAACCCCGGGGTTTT\n"
                    ">HIV__b\nAAAACCCCGGGGTTTA\n"
                    ">SIV__a\nGGGGTTTTAAAACCCC\n"
                    ">SIV__b\nGGGGTTTTAAAACCCC\n"
                ),
            }
            for gene, alignment in alignments.items():
                (results_dir / f"{gene}.aln.fasta").write_text(
                    alignment,
                    encoding="utf-8",
                )

            command = [
                sys.executable,
                str(SCRIPT),
                "--results-dir",
                str(results_dir),
                "--block-size",
                "4",
                "--bootstrap-replicates",
                "100",
                "--seed",
                "7",
            ]
            subprocess.run(command, check=True, capture_output=True, text=True)
            output_path = results_dir / "statistical_diversity_comparisons.csv"
            first_output = output_path.read_bytes()
            subprocess.run(command, check=True, capture_output=True, text=True)
            self.assertEqual(first_output, output_path.read_bytes())

            with output_path.open(newline="", encoding="utf-8") as output_file:
                rows = list(csv.DictReader(output_file))
            self.assertEqual(len(rows), 9)
            comparison_types = {row["comparison_type"] for row in rows}
            self.assertEqual(
                comparison_types,
                {"viral_groups_within_gene", "genes_within_viral_group"},
            )
            for row in rows:
                self.assertLessEqual(
                    float(row["p_value"]),
                    float(row["holm_adjusted_p_value"]),
                )
                self.assertEqual(row["valid_bootstrap_replicates"], "100")

            env_group_comparison = next(
                row
                for row in rows
                if row["comparison_type"] == "viral_groups_within_gene"
                and row["gene_a"] == "env"
                and row["group_a"] == "HIV"
                and row["group_b"] == "SIV"
            )
            self.assertAlmostEqual(
                float(env_group_comparison["pi_per_site_a"]),
                0.125,
            )
            self.assertEqual(float(env_group_comparison["pi_per_site_b"]), 0.0)
            self.assertAlmostEqual(
                float(env_group_comparison["difference_a_minus_b"]),
                0.125,
            )

    def test_missing_alignment_fails_explicitly(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--results-dir",
                    temporary_directory,
                    "--bootstrap-replicates",
                    "10",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Required alignment not found", result.stderr)

    def test_plot_script_renders_comparison_intervals_and_significance(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary_path = Path(temporary_directory)
            input_path = temporary_path / "comparisons.csv"
            output_path = temporary_path / "plot.svg"
            rows = [
                {
                    "comparison_type": "viral_groups_within_gene",
                    "gene_a": "env",
                    "group_a": "HIV",
                    "gene_b": "env",
                    "group_b": "SIV",
                    "difference_a_minus_b": "0.12",
                    "ci95_lower": "0.04",
                    "ci95_upper": "0.20",
                    "holm_adjusted_p_value": "0.01",
                },
                {
                    "comparison_type": "genes_within_viral_group",
                    "gene_a": "env",
                    "group_a": "HIV",
                    "gene_b": "gag",
                    "group_b": "HIV",
                    "difference_a_minus_b": "-0.03",
                    "ci95_lower": "-0.08",
                    "ci95_upper": "0.02",
                    "holm_adjusted_p_value": "0.20",
                },
            ]
            with input_path.open("w", newline="", encoding="utf-8") as input_file:
                writer = csv.DictWriter(input_file, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)

            subprocess.run(
                [
                    sys.executable,
                    str(PLOT_SCRIPT),
                    "--input",
                    str(input_path),
                    "--output",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            plot = output_path.read_text(encoding="utf-8")
            self.assertIn("Statistical diversity comparisons", plot)
            self.assertIn("HIV vs SIV", plot)
            self.assertIn("env vs gag", plot)
            self.assertIn("#b91c1c", plot)
            self.assertIn("#2563eb", plot)


if __name__ == "__main__":
    unittest.main()
