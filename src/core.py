"""Synthetic point-treatment analysis; no clinical evidence or real data."""

from dataclasses import dataclass
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COHORT_SQL = ROOT / "queries" / "build_cohort.sql"
DEFAULT_SEED = 20260907
DEFAULT_CLIP = 0.01


class DataQualityError(ValueError):
    """The supported one-row-per-person data contract was violated."""


class PropensityFitError(ValueError):
    """A finite, converged unpenalized propensity model could not be fitted."""


@dataclass(frozen=True)
class Estimate:
    summary: dict[str, float]
    patients: pd.DataFrame
    balance: pd.DataFrame
    descriptives: pd.DataFrame


def expit(x: np.ndarray) -> np.ndarray:
    """Numerically stable inverse logit."""
    return np.exp(-np.logaddexp(0.0, -np.asarray(x, dtype=float)))


def make_synthetic_database(n: int = 3_000, seed: int = DEFAULT_SEED) -> sqlite3.Connection:
    """Independent people, measured confounding, and complete follow-up.

    NULL outcome_day means observed event-free follow-up, NOT missing outcome
    ascertainment. No real patient records are accepted by the CLI.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 10:
        raise ValueError("n must be an integer of at least 10")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    rng = np.random.default_rng(seed)
    patient_id = np.arange(1, n + 1)
    age = np.clip(np.rint(rng.normal(66, 10, n)), 40, 90).astype(int)
    comorbidity = rng.binomial(1, expit(-0.8 + 0.035 * (age - 60)))
    treatment = rng.binomial(1, expit(-0.25 + 0.035 * (age - 65) + 0.9 * comorbidity))
    outcome = rng.binomial(
        1, expit(-2.7 + 0.035 * (age - 65) + 0.9 * comorbidity - 0.45 * treatment)
    )
    outcome_day = np.where(outcome == 1, rng.integers(1, 181, n), None)
    enrollment = rng.integers(300, 731, n)
    prior_use = rng.binomial(1, 0.08, n)
    db = sqlite3.connect(":memory:")
    db.execute("PRAGMA foreign_keys = ON")
    db.executescript("""
        CREATE TABLE patients (
            patient_id INTEGER PRIMARY KEY,
            age INTEGER NOT NULL CHECK (age BETWEEN 0 AND 120),
            comorbidity INTEGER NOT NULL CHECK (comorbidity IN (0, 1))
        );
        CREATE TABLE index_exposures (
            patient_id INTEGER PRIMARY KEY REFERENCES patients(patient_id),
            treatment TEXT NOT NULL CHECK (treatment IN ('A', 'B')),
            prior_enrollment_days INTEGER NOT NULL CHECK (prior_enrollment_days >= 0),
            prior_class_use INTEGER NOT NULL CHECK (prior_class_use IN (0, 1))
        );
        CREATE TABLE outcomes (
            patient_id INTEGER PRIMARY KEY REFERENCES patients(patient_id),
            outcome_day INTEGER CHECK (outcome_day >= 1),
            followup_days INTEGER NOT NULL CHECK (followup_days >= 0)
        );
    """)
    db.executemany("INSERT INTO patients VALUES (?, ?, ?)",
                   zip(patient_id.tolist(), age.tolist(), comorbidity.tolist()))
    db.executemany("INSERT INTO index_exposures VALUES (?, ?, ?, ?)",
                   zip(patient_id.tolist(), np.where(treatment == 1, "A", "B").tolist(),
                       enrollment.tolist(), prior_use.tolist()))
    db.executemany("INSERT INTO outcomes VALUES (?, ?, ?)",
                   zip(patient_id.tolist(), outcome_day.tolist(), [180] * n))
    return db


def validate_database(db: sqlite3.Connection) -> None:
    """Reject broken joins, ambiguous NULLs, and unsupported incomplete follow-up."""
    try:
        tables = {name: pd.read_sql_query(f"SELECT * FROM {name}", db)
                  for name in ("patients", "index_exposures", "outcomes")}
        expected = {
            "patients": {"patient_id", "age", "comorbidity"},
            "index_exposures": {"patient_id", "treatment", "prior_enrollment_days", "prior_class_use"},
            "outcomes": {"patient_id", "outcome_day", "followup_days"},
        }
        for name, frame in tables.items():
            if not expected[name].issubset(frame.columns):
                raise DataQualityError(f"{name}: missing required columns")
            if frame[list(expected[name] - {"outcome_day"})].isna().any().any():
                raise DataQualityError(f"{name}: missing required values")
            if frame.empty or frame.patient_id.duplicated().any():
                raise DataQualityError(f"{name}: expected nonempty, unique person IDs")
            numeric = sorted(expected[name] - {"outcome_day", "treatment"})
            values = frame[numeric].to_numpy(dtype=float)
            if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
                raise DataQualityError(f"{name}: expected finite integer codes/days")
            if not frame.patient_id.gt(0).all():
                raise DataQualityError(f"{name}: person IDs must be positive")
        ids = set(tables["patients"].patient_id)
        for name in ("index_exposures", "outcomes"):
            if set(tables[name].patient_id) != ids:
                raise DataQualityError(f"{name}: missing or orphaned person IDs")
        patients, exposures, outcomes = (tables[name] for name in expected)
        if not patients.age.between(0, 120).all() or not patients.comorbidity.isin([0, 1]).all():
            raise DataQualityError("Invalid baseline covariates")
        if (not exposures.treatment.isin(["A", "B"]).all()
                or not exposures.prior_class_use.isin([0, 1]).all()
                or not exposures.prior_enrollment_days.ge(0).all()):
            raise DataQualityError("Invalid exposure or enrollment values")
        if not outcomes.followup_days.ge(180).all():
            raise DataQualityError("Incomplete follow-up is unsupported; do not recode it as no event")
        days = outcomes.outcome_day.dropna()
        if not np.isfinite(days).all() or not (days.ge(1) & days.eq(np.floor(days))).all():
            raise DataQualityError("Event days must be positive integers or explicit no-event NULL")
        if (outcomes.outcome_day > outcomes.followup_days).any():
            raise DataQualityError("Event occurs after observed follow-up")
    except (pd.errors.DatabaseError, sqlite3.Error, TypeError) as exc:
        raise DataQualityError(f"Invalid source schema or values: {exc}") from exc


def cohort_flow(db: sqlite3.Connection) -> pd.DataFrame:
    """Sequential exclusions; counts add up without double-counting."""
    validate_database(db)
    source = pd.read_sql_query("""
        SELECT p.age, e.prior_enrollment_days, e.prior_class_use
        FROM patients p JOIN index_exposures e USING (patient_id)
    """, db)
    rows = [{"step": "Generated people", "remaining": len(source), "excluded": 0}]
    conditions = [
        ("Age 50-90", source.age.between(50, 90)),
        ("Prior enrollment >=365 days", source.prior_enrollment_days.ge(365)),
        ("No prior class use", source.prior_class_use.eq(0)),
    ]
    keep = pd.Series(True, index=source.index)
    for name, condition in conditions:
        previous = int(keep.sum())
        keep &= condition
        rows.append({"step": name, "remaining": int(keep.sum()),
                     "excluded": previous - int(keep.sum())})
    return pd.DataFrame(rows)


def build_cohort(db: sqlite3.Connection) -> pd.DataFrame:
    validate_database(db)
    cohort = pd.read_sql_query(COHORT_SQL.read_text(encoding="utf-8"), db)
    validate_cohort(cohort)
    return cohort


def validate_cohort(cohort: pd.DataFrame) -> None:
    required = ["patient_id", "age", "comorbidity", "treated", "outcome"]
    if not set(required).issubset(cohort.columns):
        raise DataQualityError("Cohort is missing required columns")
    if cohort.empty or cohort[required].isna().any().any():
        raise DataQualityError("Cohort is empty or has missing required values")
    if cohort.patient_id.duplicated().any():
        raise DataQualityError("Join multiplied person rows")
    if any(not pd.api.types.is_numeric_dtype(cohort[column]) for column in required):
        raise DataQualityError("Cohort columns must have numeric dtypes, not numeric strings")
    try:
        if not np.isfinite(cohort[required].to_numpy(dtype=float)).all():
            raise DataQualityError("Cohort contains nonfinite values")
    except (TypeError, ValueError) as exc:
        raise DataQualityError("Cohort values must be finite numeric values") from exc
    if not cohort.age.between(50, 90).all():
        raise DataQualityError("Cohort violates age eligibility")
    if not (cohort.patient_id.gt(0) & cohort.patient_id.eq(np.floor(cohort.patient_id))).all():
        raise DataQualityError("Person IDs must be positive integers")
    for column in ("comorbidity", "treated", "outcome"):
        if not cohort[column].isin([0, 1]).all():
            raise DataQualityError(f"{column} must be binary")
    if cohort.treated.nunique() != 2:
        raise DataQualityError("Both treatment groups are required")


def fit_propensity_score(cohort: pd.DataFrame, max_iter: int = 100) -> np.ndarray:
    """Unpenalized logistic MLE with rank, separation and convergence guards.

    Return RAW probabilities; clipping is a separate reported policy. A damped
    Newton step cannot lower log likelihood. Outcome is never a feature.
    """
    try:
        age = cohort.age.to_numpy(dtype=float)
        comorbidity = cohort.comorbidity.to_numpy(dtype=float)
        y = cohort.treated.to_numpy(dtype=float)
    except (AttributeError, TypeError, ValueError) as exc:
        raise PropensityFitError("Missing or invalid propensity inputs") from exc
    if (len(y) < 4 or not np.isfinite(np.column_stack([age, comorbidity, y])).all()
            or not np.isin(comorbidity, [0, 1]).all() or set(y) != {0, 1}):
        raise PropensityFitError("Finite inputs and both treatment groups are required")
    scale = age.std(ddof=0)
    if scale <= 0:
        raise PropensityFitError("Age has no variation")
    x = np.column_stack([np.ones(len(y)), (age - age.mean()) / scale, comorbidity])
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise PropensityFitError("Rank-deficient design matrix")
    beta = np.zeros(x.shape[1])

    def likelihood(coefficients: np.ndarray) -> float:
        eta = x @ coefficients
        return float(np.sum(y * eta - np.logaddexp(0.0, eta)))

    for _ in range(max_iter):
        p = expit(x @ beta)
        information = x.T @ (x * (p * (1 - p))[:, None])
        if np.linalg.cond(information) > 1e12 or np.max(np.abs(beta)) > 30:
            raise PropensityFitError("Separation or ill-conditioned information matrix")
        try:
            step = np.linalg.solve(information, x.T @ (y - p))
        except np.linalg.LinAlgError as exc:
            raise PropensityFitError("Could not solve the propensity information matrix") from exc
        fraction = 1.0
        current = likelihood(beta)
        while likelihood(beta + fraction * step) < current - 1e-10:
            fraction *= 0.5
            if fraction < 1e-8:
                raise PropensityFitError("Newton line search failed")
        beta += fraction * step
        if np.max(np.abs(fraction * step)) < 1e-9:
            result = expit(x @ beta)
            if np.min(result) < 1e-8 or np.max(result) > 1 - 1e-8:
                raise PropensityFitError("Near-separated propensity fit; analysis refused")
            return result
    raise PropensityFitError("Propensity model did not converge")


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    if (values.ndim != 1 or values.shape != weights.shape or values.size == 0
            or not np.isfinite(values).all() or not np.isfinite(weights).all()
            or not (weights > 0).all()):
        raise ValueError("Expected aligned finite vectors and strictly positive weights")
    return float(np.average(values, weights=weights))


def smd(values: np.ndarray, treated: np.ndarray, weights: np.ndarray) -> float:
    """Absolute SMD with FIXED unweighted pooled population-SD denominator.

    This isolates changes in means, rather than in the denominator. For binary
    covariates the group variance is the Bernoulli variance p(1-p).
    """
    values, treated, weights = np.asarray(values, float), np.asarray(treated), np.asarray(weights, float)
    if values.shape != treated.shape or values.shape != weights.shape or set(treated) != {0, 1}:
        raise ValueError("Aligned vectors and two treatment groups are required")
    means, variances = [], []
    for group in (1, 0):
        mask = treated == group
        means.append(weighted_mean(values[mask], weights[mask]))
        variances.append(float(values[mask].var(ddof=0)))
    denominator = float(np.sqrt(np.mean(variances)))
    if denominator == 0:
        if means[0] == means[1]:
            return 0.0
        raise ValueError("SMD is undefined for distinct group constants")
    return float(abs(means[0] - means[1]) / denominator)


def oracle_risk_difference(cohort: pd.DataFrame) -> float:
    """Finite-cohort expected A-B contrast from the KNOWN synthetic mechanism.

    These coefficients are simulation parameters, not a fitted outcome model
    or observed counterfactuals. Never apply this benchmark to real records.
    """
    baseline = -2.7 + 0.035 * (cohort.age.to_numpy() - 65) + 0.9 * cohort.comorbidity.to_numpy()
    return float(np.mean(expit(baseline - 0.45) - expit(baseline)))


def _estimate(cohort: pd.DataFrame, clip: float, quantile: float | None) -> Estimate:
    """Internal estimator also accepts repeated IDs in bootstrap samples."""
    if not np.isfinite(clip) or not 0 < clip < 0.5:
        raise ValueError("Propensity clip must lie strictly between 0 and 0.5")
    if quantile is not None and (not np.isfinite(quantile) or not 0 < quantile < 0.5):
        raise ValueError("Winsorization tail quantile must lie between 0 and 0.5")
    treated = cohort.treated.to_numpy()
    outcome = cohort.outcome.to_numpy(dtype=float)
    raw = fit_propensity_score(cohort)
    propensity = np.clip(raw, clip, 1 - clip)
    prevalence = float(treated.mean())
    weights = np.where(treated == 1, prevalence / propensity, (1 - prevalence) / (1 - propensity))
    changed = 0
    if quantile is not None:
        lower, upper = np.quantile(weights, [quantile, 1 - quantile])
        capped = np.clip(weights, lower, upper)
        changed = int(np.count_nonzero(capped != weights))
        weights = capped
    ones = np.ones(len(cohort))
    rows = []
    for name in ("age", "comorbidity"):
        values = cohort[name].to_numpy(dtype=float)
        rows.append({"covariate": name, "smd_before": smd(values, treated, ones),
                     "smd_after": smd(values, treated, weights)})
    balance = pd.DataFrame(rows)
    summary = {
        "n": float(len(cohort)), "propensity_min_raw": float(raw.min()),
        "propensity_max_raw": float(raw.max()),
        "propensity_clipped_count": float(np.count_nonzero(raw != propensity)),
        "weight_winsorized_count": float(changed),
        "effective_sample_size": float(weights.sum() ** 2 / np.square(weights).sum()),
        "maximum_weight": float(weights.max()), "mean_weight": float(weights.mean()),
        "weight_p99": float(np.quantile(weights, 0.99)),
    }
    descriptives = []
    for group, label in ((1, "A"), (0, "B")):
        mask = treated == group
        summary[f"risk_{label}_crude"] = weighted_mean(outcome[mask], ones[mask])
        summary[f"risk_{label}_weighted"] = weighted_mean(outcome[mask], weights[mask])
        summary[f"ess_{label}"] = float(weights[mask].sum() ** 2 / np.square(weights[mask]).sum())
        descriptives.append({
            "treatment": label, "n": int(mask.sum()), "events": int(outcome[mask].sum()),
            "age_mean": float(cohort.age.to_numpy()[mask].mean()),
            "age_sd": float(cohort.age.to_numpy()[mask].std(ddof=1)),
            "comorbidity_prevalence": float(cohort.comorbidity.to_numpy()[mask].mean()),
            "crude_risk": summary[f"risk_{label}_crude"],
            "weighted_risk": summary[f"risk_{label}_weighted"],
            "effective_sample_size": summary[f"ess_{label}"],
        })
    summary["crude_risk_difference"] = summary["risk_A_crude"] - summary["risk_B_crude"]
    summary["weighted_risk_difference"] = summary["risk_A_weighted"] - summary["risk_B_weighted"]
    for row in rows:
        for when in ("before", "after"):
            summary[f"{row['covariate']}_smd_{when}"] = row[f"smd_{when}"]
    patients = cohort.assign(propensity_raw=raw, propensity_used=propensity, weight=weights)
    return Estimate(summary, patients, balance, pd.DataFrame(descriptives))


def estimate(cohort: pd.DataFrame, clip: float = DEFAULT_CLIP,
             quantile: float | None = None) -> Estimate:
    validate_cohort(cohort)
    return _estimate(cohort, clip, quantile)


def analyze(cohort: pd.DataFrame) -> dict[str, float]:
    """Compatibility entry point for the original demonstration."""
    return estimate(cohort).summary


def bootstrap_estimates(cohort: pd.DataFrame, resamples: int = 1_000,
                        seed: int = 8675309, clip: float = DEFAULT_CLIP) -> pd.DataFrame:
    """Patient-level pairs bootstrap; refit propensity AND weights every time.

    If any scheduled fit fails, raise rather than silently redrawing or issuing
    a confidence interval from a selected subset of replicates.
    """
    validate_cohort(cohort)
    if isinstance(resamples, bool) or not isinstance(resamples, int) or resamples < 100:
        raise ValueError("Use at least 100 bootstrap resamples (1,000 for the showcase)")
    rng = np.random.default_rng(seed)
    rows = []
    for replicate in range(resamples):
        sample = cohort.iloc[rng.integers(0, len(cohort), len(cohort))]
        try:
            result = _estimate(sample, clip, None).summary
        except (ValueError, np.linalg.LinAlgError) as exc:
            raise PropensityFitError(
                f"Bootstrap replicate {replicate} failed; CI not produced: {exc}"
            ) from exc
        rows.append({"replicate": replicate,
                     "crude_risk_difference": result["crude_risk_difference"],
                     "weighted_risk_difference": result["weighted_risk_difference"]})
    return pd.DataFrame(rows)
