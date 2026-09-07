"""Minimal, reproducible target-trial-emulation demonstration."""

from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
COHORT_SQL = ROOT / "queries" / "build_cohort.sql"


def expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def make_synthetic_database(n: int = 3_000, seed: int = 20260907) -> sqlite3.Connection:
    """Create claims-like tables with measured treatment confounding."""
    rng = np.random.default_rng(seed)
    patient_id = np.arange(1, n + 1)
    age = np.clip(np.rint(rng.normal(66, 10, n)), 40, 90).astype(int)
    comorbidity = rng.binomial(1, expit(-0.8 + 0.035 * (age - 60)))
    treatment = rng.binomial(1, expit(-0.25 + 0.035 * (age - 65) + 0.9 * comorbidity))
    outcome = rng.binomial(
        1,
        expit(-2.7 + 0.035 * (age - 65) + 0.9 * comorbidity - 0.45 * treatment),
    )
    outcome_day = np.where(outcome == 1, rng.integers(1, 181, n), None)
    enrollment = rng.integers(300, 731, n)
    prior_use = rng.binomial(1, 0.08, n)

    db = sqlite3.connect(":memory:")
    db.executescript(
        """
        CREATE TABLE patients (
            patient_id INTEGER PRIMARY KEY,
            age INTEGER NOT NULL,
            comorbidity INTEGER NOT NULL
        );
        CREATE TABLE index_exposures (
            patient_id INTEGER PRIMARY KEY,
            treatment TEXT NOT NULL,
            prior_enrollment_days INTEGER NOT NULL,
            prior_class_use INTEGER NOT NULL
        );
        CREATE TABLE outcomes (
            patient_id INTEGER PRIMARY KEY,
            outcome_day INTEGER
        );
        """
    )
    db.executemany("INSERT INTO patients VALUES (?, ?, ?)", zip(patient_id.tolist(), age.tolist(), comorbidity.tolist()))
    db.executemany(
        "INSERT INTO index_exposures VALUES (?, ?, ?, ?)",
        zip(patient_id.tolist(), np.where(treatment == 1, "A", "B").tolist(), enrollment.tolist(), prior_use.tolist()),
    )
    db.executemany("INSERT INTO outcomes VALUES (?, ?)", zip(patient_id.tolist(), outcome_day.tolist()))
    return db


def build_cohort(db: sqlite3.Connection) -> pd.DataFrame:
    return pd.read_sql_query(COHORT_SQL.read_text(encoding="utf-8"), db)


def fit_propensity_score(cohort: pd.DataFrame) -> np.ndarray:
    """Fit treatment ~ age + comorbidity with Newton's method."""
    age_z = (cohort["age"].to_numpy() - cohort["age"].mean()) / cohort["age"].std()
    x = np.column_stack([np.ones(len(cohort)), age_z, cohort["comorbidity"]])
    y = cohort["treated"].to_numpy(dtype=float)
    beta = np.zeros(x.shape[1])
    for _ in range(25):
        p = np.clip(expit(x @ beta), 1e-6, 1 - 1e-6)
        information = x.T @ (x * (p * (1 - p))[:, None])
        step = np.linalg.solve(information, x.T @ (y - p))
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return np.clip(expit(x @ beta), 0.01, 0.99)


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(values, weights=weights))


def smd(values: np.ndarray, treated: np.ndarray, weights: np.ndarray) -> float:
    """Absolute standardized mean difference using weighted group moments."""
    means, variances = [], []
    for group in (1, 0):
        mask = treated == group
        mean = weighted_mean(values[mask], weights[mask])
        means.append(mean)
        variances.append(weighted_mean((values[mask] - mean) ** 2, weights[mask]))
    return abs(means[0] - means[1]) / np.sqrt((variances[0] + variances[1]) / 2)


def analyze(cohort: pd.DataFrame) -> dict[str, float]:
    treated = cohort["treated"].to_numpy()
    outcome = cohort["outcome"].to_numpy(dtype=float)
    propensity = fit_propensity_score(cohort)
    treatment_prevalence = treated.mean()
    weights = np.where(
        treated == 1,
        treatment_prevalence / propensity,
        (1 - treatment_prevalence) / (1 - propensity),
    )

    def risk(group: int, use_weights: np.ndarray) -> float:
        mask = treated == group
        return weighted_mean(outcome[mask], use_weights[mask])

    unweighted = np.ones(len(cohort))
    return {
        "n": float(len(cohort)),
        "crude_risk_difference": risk(1, unweighted) - risk(0, unweighted),
        "weighted_risk_difference": risk(1, weights) - risk(0, weights),
        "age_smd_before": smd(cohort["age"].to_numpy(), treated, unweighted),
        "age_smd_after": smd(cohort["age"].to_numpy(), treated, weights),
        "comorbidity_smd_before": smd(cohort["comorbidity"].to_numpy(), treated, unweighted),
        "comorbidity_smd_after": smd(cohort["comorbidity"].to_numpy(), treated, weights),
        "effective_sample_size": float(weights.sum() ** 2 / np.square(weights).sum()),
        "maximum_weight": float(weights.max()),
    }


def main() -> None:
    cohort = build_cohort(make_synthetic_database())
    results = analyze(cohort)
    for name, value in results.items():
        print(f"{name}: {value:.4f}")

    assert results["n"] > 1_000
    assert results["age_smd_after"] < 0.05
    assert results["comorbidity_smd_after"] < 0.05
    assert results["weighted_risk_difference"] < 0
    print("Self-check passed.")


if __name__ == "__main__":
    main()

