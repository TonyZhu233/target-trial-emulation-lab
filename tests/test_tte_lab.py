import unittest

from src.tte_lab import analyze, build_cohort, make_synthetic_database


class TestTTELab(unittest.TestCase):
    def test_weighting_improves_measured_balance(self) -> None:
        results = analyze(build_cohort(make_synthetic_database()))
        self.assertLess(results["age_smd_after"], results["age_smd_before"])
        self.assertLess(results["comorbidity_smd_after"], results["comorbidity_smd_before"])
        self.assertLess(results["weighted_risk_difference"], 0)


if __name__ == "__main__":
    unittest.main()

