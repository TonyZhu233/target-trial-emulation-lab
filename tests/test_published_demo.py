"""The published result must match its provenance, not just look plausible."""

import hashlib
import json
import unittest

import numpy as np
import pandas as pd

from src.core import ROOT


class PublishedDemoTests(unittest.TestCase):
    def test_all_published_hashes_match(self) -> None:
        demo = ROOT / "results" / "demo"
        manifest = json.loads((demo / "manifest.json").read_text(encoding="utf-8"))
        for name, expected in manifest["artifact_sha256"].items():
            with self.subTest(artifact=name):
                self.assertEqual(hashlib.sha256((demo / name).read_bytes()).hexdigest(), expected)
        for name, expected in manifest["source_sha256"].items():
            with self.subTest(source=name):
                canonical = (ROOT / name).read_text(encoding="utf-8").encode("utf-8")
                self.assertEqual(hashlib.sha256(canonical).hexdigest(), expected)

    def test_published_intervals_and_counts_reconcile(self) -> None:
        demo = ROOT / "results" / "demo"
        manifest = json.loads((demo / "manifest.json").read_text(encoding="utf-8"))
        boot = pd.read_csv(demo / "bootstrap.csv")
        effects = pd.read_csv(demo / "effects.csv")
        self.assertEqual(len(boot), manifest["config"]["bootstrap"])
        self.assertEqual(len(pd.read_csv(demo / "cohort.csv")), manifest["summary"]["n"])
        self.assertEqual(len(pd.read_csv(demo / "source_patients.csv")), manifest["config"]["n"])
        simulations = pd.read_csv(demo / "simulation_runs.csv")
        self.assertEqual(len(simulations), manifest["config"]["simulations"])
        for label, key in (("Crude", "crude_risk_difference"), ("IPTW", "weighted_risk_difference")):
            row = effects.loc[effects.method == label].iloc[0]
            np.testing.assert_allclose([row.ci_lower, row.ci_upper],
                                       np.quantile(boot[key], [0.025, 0.975]), rtol=1e-8, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
