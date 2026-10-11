"""Real CLI smoke test, output integrity, and multiseed validation."""

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from xml.etree import ElementTree

import pandas as pd

from src.core import ROOT, PropensityFitError
from src.report import RunConfig, run_pipeline, simulation_summary, simulation_validation


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        (ROOT / "results").mkdir(exist_ok=True)

    def test_real_cli_and_manifest(self) -> None:
        # A temporary direct child of the validated generated-results directory.
        with tempfile.TemporaryDirectory(prefix="test-", dir=ROOT / "results") as temporary:
            output = Path(temporary).resolve()
            self.assertEqual(output.parent, (ROOT / "results").resolve())
            process = subprocess.run(
                [sys.executable, "-m", "src.run", "--bootstrap", "100", "--simulations", "2",
                 "--output-dir", str(output)], cwd=ROOT, capture_output=True, text=True,
                encoding="utf-8", timeout=90,
            )
            self.assertEqual(process.returncode, 0, process.stderr)
            manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["summary"]["n"], 2196)
            self.assertEqual(manifest["simulation_summary"]["successful"], 2)
            for name, expected in manifest["artifact_sha256"].items():
                self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), expected)
            for name, expected in manifest["source_sha256"].items():
                canonical = (ROOT / name).read_text(encoding="utf-8").encode("utf-8")
                self.assertEqual(hashlib.sha256(canonical).hexdigest(), expected)
            for filename in ("balance", "overlap_weights", "effect"):
                self.assertGreater((output / f"{filename}.png").stat().st_size, 1000)
                self.assertIn("<svg", (output / f"{filename}.svg").read_text(encoding="utf-8"))
                ElementTree.parse(output / f"{filename}.svg")
            self.assertEqual(len(pd.read_csv(output / "bootstrap.csv")), 100)
            self.assertIn("Not clinical evidence", (output / "REPORT.md").read_text(encoding="utf-8"))
            with self.assertRaisesRegex(ValueError, "nonempty"):
                run_pipeline(output, RunConfig(bootstrap=100, simulations=0))

    def test_dangerous_output_path_is_refused_before_write(self) -> None:
        for destination in (ROOT, ROOT / "src", ROOT / "results", ROOT / "results" / "nested" / "child"):
            with self.assertRaisesRegex(ValueError, "direct subdirectory"):
                run_pipeline(destination, RunConfig(bootstrap=100, simulations=0))

    def test_invalid_configuration(self) -> None:
        for settings in ({"simulations": -1}, {"bootstrap": 99}, {"seed": -1}, {"n": True}):
            with self.assertRaises(ValueError):
                RunConfig(**settings)

    def test_all_simulation_seeds_including_failures_are_retained(self) -> None:
        with patch("src.report.estimate", side_effect=PropensityFitError("forced failure")):
            frame = simulation_validation(n=3000, repeats=3, seed=100)
        self.assertEqual(frame.seed.tolist(), [100, 101, 102])
        self.assertEqual(simulation_summary(frame)["failed"], 3)

    def test_multiseed_behavior_against_known_mechanism(self) -> None:
        frame = simulation_validation(n=3000, repeats=20, seed=20260908)
        summary = simulation_summary(frame)
        self.assertEqual(summary["failed"], 0)
        self.assertLess(summary["weighted_rmse"], summary["crude_rmse"])
        self.assertLess(abs(summary["weighted_mean_error"]), 0.01)
        self.assertTrue((frame.maximum_smd_after < 0.05).all())
        # No assertion about every estimate's sign or statistical significance.

    def test_no_simulations_is_explicit(self) -> None:
        summary = simulation_summary(simulation_validation(n=3000, repeats=0, seed=0))
        self.assertEqual(summary, {"scheduled": 0, "successful": 0, "failed": 0})


if __name__ == "__main__":
    unittest.main()
