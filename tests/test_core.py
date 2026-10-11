"""Data-contract, numerical, estimator, and inference regression tests."""

import sqlite3
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from src import core


class DatabaseTests(unittest.TestCase):
    def setUp(self) -> None:
        self.db = core.make_synthetic_database()

    def tearDown(self) -> None:
        self.db.close()

    def test_generator_is_reproducible(self) -> None:
        other = core.make_synthetic_database()
        try:
            pd.testing.assert_frame_equal(core.build_cohort(self.db), core.build_cohort(other))
        finally:
            other.close()

    def test_sql_matches_independent_boolean_filter(self) -> None:
        source = pd.read_sql_query("""
            SELECT p.patient_id,p.age,p.comorbidity,e.treatment,
                   e.prior_enrollment_days,e.prior_class_use,o.outcome_day
            FROM patients p JOIN index_exposures e USING(patient_id)
            JOIN outcomes o USING(patient_id) ORDER BY p.patient_id
        """, self.db)
        expected = source.loc[source.age.between(50, 90)
                              & source.prior_enrollment_days.ge(365)
                              & source.prior_class_use.eq(0)].copy()
        expected["treated"] = expected.treatment.eq("A").astype(int)
        expected["outcome"] = expected.outcome_day.between(1, 180).astype(int)
        expected = expected[["patient_id", "age", "comorbidity", "treated", "outcome"]].reset_index(drop=True)
        pd.testing.assert_frame_equal(core.build_cohort(self.db), expected)

    def test_boundary_days_and_eligibility(self) -> None:
        self.db.execute("UPDATE patients SET age=60 WHERE patient_id <= 6")
        self.db.execute("UPDATE index_exposures SET prior_enrollment_days=365,prior_class_use=0 WHERE patient_id <= 6")
        self.db.execute("UPDATE outcomes SET followup_days=200 WHERE patient_id <= 6")
        for person, day in ((1, 1), (2, 180), (3, 181), (4, None)):
            self.db.execute("UPDATE outcomes SET outcome_day=? WHERE patient_id=?", (day, person))
        self.db.execute("UPDATE patients SET age=49 WHERE patient_id=5")
        self.db.execute("UPDATE index_exposures SET prior_enrollment_days=364 WHERE patient_id=6")
        cohort = core.build_cohort(self.db).set_index("patient_id")
        self.assertEqual(cohort.loc[[1, 2, 3, 4], "outcome"].tolist(), [1, 1, 0, 0])
        self.assertNotIn(5, cohort.index)
        self.assertNotIn(6, cohort.index)
        self.db.execute("UPDATE patients SET age=50 WHERE patient_id=5")
        self.db.execute("UPDATE patients SET age=90 WHERE patient_id=6")
        self.db.execute("UPDATE index_exposures SET prior_enrollment_days=365 WHERE patient_id=6")
        self.assertTrue({5, 6}.issubset(set(core.build_cohort(self.db).patient_id)))

    def test_flow_reconciles(self) -> None:
        flow = core.cohort_flow(self.db)
        self.assertEqual(flow.iloc[-1].remaining, len(core.build_cohort(self.db)))
        self.assertEqual(flow.excluded.sum() + flow.iloc[-1].remaining, 3000)

    def test_missing_outcome_is_not_reclassified_as_no_event(self) -> None:
        self.db.execute("DELETE FROM outcomes WHERE patient_id=1")
        with self.assertRaisesRegex(core.DataQualityError, "missing or orphaned"):
            core.build_cohort(self.db)

    def test_missing_exposure_is_not_silently_dropped(self) -> None:
        self.db.execute("DELETE FROM index_exposures WHERE patient_id=1")
        with self.assertRaises(core.DataQualityError):
            core.build_cohort(self.db)

    def test_incomplete_followup_is_rejected(self) -> None:
        self.db.execute("UPDATE outcomes SET followup_days=179,outcome_day=NULL WHERE patient_id=1")
        with self.assertRaisesRegex(core.DataQualityError, "Incomplete follow-up"):
            core.build_cohort(self.db)

    def test_duplicate_source_rows_are_rejected_even_without_constraints(self) -> None:
        self.db.execute("ALTER TABLE outcomes RENAME TO original_outcomes")
        self.db.execute("CREATE TABLE outcomes AS SELECT * FROM original_outcomes")
        self.db.execute("INSERT INTO outcomes SELECT * FROM original_outcomes WHERE patient_id=1")
        with self.assertRaisesRegex(core.DataQualityError, "unique person IDs"):
            core.build_cohort(self.db)

    def test_database_constraints(self) -> None:
        with self.assertRaises(sqlite3.IntegrityError):
            self.db.execute("UPDATE index_exposures SET treatment='C' WHERE patient_id=1")
        with self.assertRaises(sqlite3.IntegrityError):
            self.db.execute("INSERT INTO outcomes VALUES (999999,NULL,180)")

    def test_invalid_source_number_and_day_are_rejected(self) -> None:
        self.db.execute("UPDATE outcomes SET followup_days=? WHERE patient_id=1", (float("inf"),))
        with self.assertRaises(core.DataQualityError):
            core.build_cohort(self.db)
        self.db.execute("UPDATE outcomes SET followup_days=180,outcome_day=1.5 WHERE patient_id=1")
        with self.assertRaises(core.DataQualityError):
            core.build_cohort(self.db)


class EstimatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        db = core.make_synthetic_database()
        try:
            cls.cohort = core.build_cohort(db)
        finally:
            db.close()

    def test_baseline_result_and_bounds(self) -> None:
        result = core.estimate(self.cohort)
        self.assertEqual(result.summary["n"], 2196)
        self.assertAlmostEqual(result.summary["weighted_risk_difference"], -0.0306, places=4)
        for group in ("A", "B"):
            self.assertTrue(0 <= result.summary[f"risk_{group}_weighted"] <= 1)
        self.assertTrue((result.patients.weight > 0).all())
        self.assertLessEqual(result.summary["effective_sample_size"], len(self.cohort))
        self.assertEqual(result.summary["propensity_clipped_count"], 0)

    def test_logistic_score_equations(self) -> None:
        p = core.fit_propensity_score(self.cohort)
        age = self.cohort.age.to_numpy(dtype=float)
        x = np.column_stack([np.ones(len(age)), (age-age.mean())/age.std(), self.cohort.comorbidity])
        self.assertLess(np.max(np.abs(x.T @ (self.cohort.treated.to_numpy()-p))), 1e-7)

    def test_outcome_is_not_a_propensity_feature(self) -> None:
        changed = self.cohort.assign(outcome=1-self.cohort.outcome)
        np.testing.assert_array_equal(core.fit_propensity_score(self.cohort), core.fit_propensity_score(changed))

    def test_invalid_cohorts(self) -> None:
        cases = [self.cohort.iloc[:0], self.cohort.drop(columns="age"),
                 pd.concat([self.cohort, self.cohort.iloc[:1]]), self.cohort.assign(outcome=2),
                 self.cohort.assign(treated=0), self.cohort.assign(age=np.nan),
                 self.cohort.assign(age=self.cohort.age.astype(str)),
                 self.cohort.assign(patient_id=0), self.cohort.assign(patient_id=self.cohort.patient_id+0.5)]
        for cohort in cases:
            with self.subTest(case=str(cohort.shape)), self.assertRaises(core.DataQualityError):
                core.estimate(cohort)

    def test_rank_deficiency_separation_and_nonconvergence(self) -> None:
        for data in (self.cohort.assign(age=60), self.cohort.assign(comorbidity=0),
                     self.cohort.assign(treated=self.cohort.comorbidity)):
            with self.assertRaises(core.PropensityFitError):
                core.fit_propensity_score(data)
        with self.assertRaisesRegex(core.PropensityFitError, "did not converge"):
            core.fit_propensity_score(self.cohort, max_iter=0)

    def test_no_events_is_a_valid_null_observation(self) -> None:
        result = core.estimate(self.cohort.assign(outcome=0))
        self.assertEqual(result.summary["weighted_risk_difference"], 0)

    def test_solver_failure_has_an_explicit_error(self) -> None:
        with patch("src.core.np.linalg.solve", side_effect=np.linalg.LinAlgError("forced failure")):
            with self.assertRaisesRegex(core.PropensityFitError, "information matrix"):
                core.fit_propensity_score(self.cohort)

    def test_inverse_logit_is_stable(self) -> None:
        with np.errstate(over="raise"):
            np.testing.assert_allclose(core.expit(np.array([-1000, 0, 1000])), [0, 0.5, 1])

    def test_weighted_mean_by_hand_and_invalid_weights(self) -> None:
        self.assertAlmostEqual(core.weighted_mean(np.array([0, 1]), np.array([1, 3])), 0.75)
        for weights in (np.array([0, 1]), np.array([-1, 1]), np.array([np.nan, 1]), np.array([1])):
            with self.assertRaises(ValueError):
                core.weighted_mean(np.array([0, 1]), weights)

    def test_normalized_iptw_by_hand(self) -> None:
        cohort = pd.DataFrame({"patient_id": [1,2,3,4], "age": [50,60,70,80],
                               "comorbidity": [0,1,0,1], "treated": [1,1,0,0], "outcome": [1,0,1,0]})
        with patch("src.core.fit_propensity_score", return_value=np.array([0.25,0.5,0.5,0.75])):
            result = core.estimate(cohort)
        self.assertAlmostEqual(result.summary["risk_A_weighted"], 2/3)
        self.assertAlmostEqual(result.summary["risk_B_weighted"], 1/3)
        self.assertAlmostEqual(result.summary["weighted_risk_difference"], 1/3)

    def test_smd_fixed_denominator_by_hand(self) -> None:
        values = np.array([1,3,2,4])
        treated = np.array([1,1,0,0])
        self.assertAlmostEqual(core.smd(values, treated, np.ones(4)), 1)
        self.assertAlmostEqual(core.smd(values, treated, np.array([1,3,1,1])), 0.5)
        self.assertEqual(core.smd(np.ones(4), treated, np.ones(4)), 0)
        with self.assertRaises(ValueError):
            core.smd(values, treated, np.ones(3))

    def test_clipping_and_winsorization_are_explicit(self) -> None:
        result = core.estimate(self.cohort, quantile=0.01)
        self.assertGreater(result.summary["weight_winsorized_count"], 0)
        self.assertLess(result.summary["maximum_weight"], core.estimate(self.cohort).summary["maximum_weight"])
        for clip in (0, 0.5, np.nan):
            with self.assertRaises(ValueError):
                core.estimate(self.cohort, clip=clip)
        for quantile in (0, 0.5, np.inf):
            with self.assertRaises(ValueError):
                core.estimate(self.cohort, quantile=quantile)

    def test_bootstrap_reproducible_and_refits_every_time(self) -> None:
        with patch("src.core.fit_propensity_score", wraps=core.fit_propensity_score) as fit:
            first = core.bootstrap_estimates(self.cohort, resamples=100, seed=5)
            self.assertEqual(fit.call_count, 100)
        second = core.bootstrap_estimates(self.cohort, resamples=100, seed=5)
        pd.testing.assert_frame_equal(first, second)
        self.assertGreater(first.weighted_risk_difference.std(), 0)

    def test_bootstrap_failure_is_not_discarded(self) -> None:
        with patch("src.core.fit_propensity_score", side_effect=core.PropensityFitError("forced failure")):
            with self.assertRaisesRegex(core.PropensityFitError, "replicate 0 failed; CI not produced"):
                core.bootstrap_estimates(self.cohort, resamples=100)

    def test_invalid_generator_and_bootstrap_arguments(self) -> None:
        for n in (0, 9, True, 30.5):
            with self.assertRaises(ValueError):
                core.make_synthetic_database(n=n)
        with self.assertRaises(ValueError):
            core.bootstrap_estimates(self.cohort, resamples=99)


if __name__ == "__main__":
    unittest.main()
