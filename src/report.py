"""Offline pipeline, transparent report, and accessible diagnostic figures."""

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use("Agg")  # Headless on Windows, Linux, and CI; no browser required.
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .core import (
    DEFAULT_CLIP, DEFAULT_SEED, ROOT, Estimate, bootstrap_estimates,
    build_cohort, cohort_flow, estimate, make_synthetic_database, oracle_risk_difference,
)

BLUE, ORANGE = "#0072B2", "#D55E00"


@dataclass(frozen=True)
class RunConfig:
    n: int = 3_000
    seed: int = DEFAULT_SEED
    bootstrap: int = 1_000
    bootstrap_seed: int = 8675309
    simulations: int = 100
    simulation_seed: int = 20260908

    def __post_init__(self) -> None:
        for name in ("n", "seed", "bootstrap", "bootstrap_seed", "simulations", "simulation_seed"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if self.n < 10 or self.bootstrap < 100:
            raise ValueError("Require n >=10 and bootstrap >=100")


def simulation_validation(n: int, repeats: int, seed: int) -> pd.DataFrame:
    """All consecutive seeds, including failures; no selection by effect sign."""
    if isinstance(repeats, bool) or not isinstance(repeats, int) or repeats < 0:
        raise ValueError("Simulation repetitions must be a non-negative integer")
    rows = []
    for index in range(repeats):
        run_seed = seed + index
        db = make_synthetic_database(n=n, seed=run_seed)
        try:
            cohort = build_cohort(db)
            result = estimate(cohort).summary
            oracle = oracle_risk_difference(cohort)
            rows.append({
                "seed": run_seed, "status": "ok", "error": "", "n": len(cohort),
                "oracle_risk_difference": oracle,
                "crude_risk_difference": result["crude_risk_difference"],
                "weighted_risk_difference": result["weighted_risk_difference"],
                "crude_error": result["crude_risk_difference"] - oracle,
                "weighted_error": result["weighted_risk_difference"] - oracle,
                "maximum_smd_after": max(result["age_smd_after"], result["comorbidity_smd_after"]),
                "effective_sample_size": result["effective_sample_size"],
            })
        except (ValueError, np.linalg.LinAlgError) as exc:
            rows.append({"seed": run_seed, "status": "failed", "error": str(exc)})
        finally:
            db.close()
    return pd.DataFrame(rows, columns=[
        "seed", "status", "error", "n", "oracle_risk_difference", "crude_risk_difference",
        "weighted_risk_difference", "crude_error", "weighted_error", "maximum_smd_after",
        "effective_sample_size",
    ])


def simulation_summary(frame: pd.DataFrame) -> dict[str, float | int]:
    successful = frame.loc[frame.status == "ok"]
    result: dict[str, float | int] = {"scheduled": len(frame), "successful": len(successful),
                                    "failed": int((frame.status == "failed").sum())}
    for method in ("crude", "weighted"):
        errors = successful[f"{method}_error"].to_numpy(dtype=float)
        if len(errors):
            result[f"{method}_mean_error"] = float(errors.mean())
            result[f"{method}_rmse"] = float(np.sqrt(np.mean(errors ** 2)))
        if len(errors) > 1:
            result[f"{method}_mean_error_mcse"] = float(errors.std(ddof=1) / np.sqrt(len(errors)))
    return result


def markdown_table(frame: pd.DataFrame) -> str:
    """A small table formatter without a hidden tabulate dependency."""
    def cell(value: object) -> str:
        if isinstance(value, (float, np.floating)):
            return f"{value:.4f}"
        return str(value).replace("|", "\\|").replace("\n", " ")
    lines = ["| " + " | ".join(map(str, frame.columns)) + " |",
             "| " + " | ".join(["---"] * len(frame.columns)) + " |"]
    lines.extend("| " + " | ".join(cell(value) for value in row) + " |"
                 for row in frame.itertuples(index=False, name=None))
    return "\n".join(lines)


def save_figures(destination: Path, result: Estimate, effects: pd.DataFrame,
                 oracle: float, resamples: int) -> list[str]:
    """OO figures, redundant color/shape/style, common scales, defined intervals."""
    filenames = []
    style = {
        "font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 13,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "svg.fonttype": "none",
        "svg.hashsalt": "tte-lab-v1",
    }

    def export(fig: plt.Figure, name: str) -> None:
        for extension in ("svg", "png"):
            filename = f"{name}.{extension}"
            path = destination / filename
            metadata = {"Date": None} if extension == "svg" else None
            fig.savefig(path, dpi=180, facecolor="white", metadata=metadata)
            if extension == "svg":
                # Mechanical XML whitespace normalization, not a data/figure edit.
                text = "\n".join(line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()) + "\n"
                path.write_text(text, encoding="utf-8", newline="\n")
            filenames.append(filename)
        plt.close(fig)

    with plt.rc_context(style):
        fig, ax = plt.subplots(figsize=(8, 3.7), layout="constrained")
        positions = np.arange(len(result.balance))
        for index, row in enumerate(result.balance.itertuples()):
            ax.plot([row.smd_before, row.smd_after], [index, index], color="0.7", lw=1)
        ax.scatter(result.balance.smd_before, positions, color=ORANGE, marker="s", s=65, label="Unweighted")
        ax.scatter(result.balance.smd_after, positions, color=BLUE, marker="o", s=65, label="IPTW")
        ax.axvline(0.1, color="0.35", ls="--", label="0.10 reference (not proof)")
        ax.set(yticks=positions, yticklabels=["Age", "Comorbidity"],
               xlabel="Absolute standardized mean difference (fixed denominator)",
               xlim=(-0.01, max(0.15, result.balance.smd_before.max() * 1.15)),
               ylim=(-0.5, 1.7), title="Measured baseline balance")
        ax.legend(loc="upper right", fontsize=9)
        export(fig, "balance")

        fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
        people = result.patients
        for group, label, color, linestyle in ((1, "A", BLUE, "-"), (0, "B", ORANGE, "--")):
            subset = people.loc[people.treated == group]
            axes[0].hist(subset.propensity_raw, bins=np.linspace(0, 1, 21), density=True,
                         histtype="step", color=color, linestyle=linestyle, linewidth=2,
                         label=f"{label} (n={len(subset)})")
            axes[1].hist(subset.weight, bins=np.linspace(people.weight.min(), people.weight.max(), 21),
                         density=True, histtype="step", color=color, linestyle=linestyle,
                         linewidth=2, label=label)
        axes[0].set(xlabel="Raw estimated P(A | baseline)", ylabel="Density", xlim=(0, 1),
                    title="Propensity-score overlap")
        axes[1].set(xlabel="Stabilized weight", ylabel="Density", title="IPTW weight distribution")
        for ax in axes:
            ax.legend(fontsize=9)
        export(fig, "overlap_weights")

        fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
        for position, row in enumerate(effects.itertuples()):
            color, marker = (ORANGE, "s") if position == 0 else (BLUE, "o")
            ax.hlines(position, row.ci_lower * 100, row.ci_upper * 100, color=color, linewidth=2)
            ax.plot(row.estimate * 100, position, color=color, marker=marker, markersize=8)
        ax.axvline(0, color="0.5", linewidth=1)
        ax.axvline(oracle * 100, color="0.2", linestyle="--", label="Known DGP benchmark (no CI)")
        ax.set(yticks=[0, 1], yticklabels=["Crude", "IPTW"], ylim=(-0.65, 1.8),
               xlabel="180-day risk difference, A minus B (percentage points)",
               title=f"Synthetic estimate and 95% percentile bootstrap CI (B={resamples:,})")
        ax.legend(loc="upper left", fontsize=9)
        export(fig, "effect")
    return filenames


def run_pipeline(destination: Path, config: RunConfig, overwrite: bool = False) -> dict[str, object]:
    """Run fully offline and export ONLY generated synthetic artifacts."""
    destination = destination.resolve()
    results_root = (ROOT / "results").resolve()
    if destination.parent != results_root:
        raise ValueError("Output must be a direct subdirectory of this repository's results/ folder")
    if destination.exists() and any(destination.iterdir()) and not overwrite:
        raise ValueError("Output is nonempty; choose another folder or explicitly use --overwrite")
    db = make_synthetic_database(config.n, config.seed)
    try:
        flow = cohort_flow(db)
        cohort = build_cohort(db)
        sources = {name: pd.read_sql_query(f"SELECT * FROM {name} ORDER BY patient_id", db)
                   for name in ("patients", "index_exposures", "outcomes")}
    finally:
        db.close()
    result = estimate(cohort)
    boot = bootstrap_estimates(cohort, config.bootstrap, config.bootstrap_seed)
    oracle = oracle_risk_difference(cohort)
    effects = []
    for label, key in (("Crude", "crude_risk_difference"), ("IPTW", "weighted_risk_difference")):
        lower, upper = np.quantile(boot[key], [0.025, 0.975])
        effects.append({"method": label, "estimate": result.summary[key],
                        "ci_lower": float(lower), "ci_upper": float(upper)})
    effects_frame = pd.DataFrame(effects)
    sensitivity = []
    for label, quantile in (("No winsorization", None), ("1st/99th percentiles", 0.01),
                            ("5th/95th percentiles", 0.05)):
        stats = estimate(cohort, quantile=quantile).summary
        sensitivity.append({"policy": label, "risk_difference": stats["weighted_risk_difference"],
                            "maximum_weight": stats["maximum_weight"],
                            "effective_sample_size": stats["effective_sample_size"],
                            "maximum_smd_after": max(stats["age_smd_after"], stats["comorbidity_smd_after"]),
                            "changed_weights": int(stats["weight_winsorized_count"])})
    sensitivity_frame = pd.DataFrame(sensitivity)
    simulations = simulation_validation(config.n, config.simulations, config.simulation_seed)
    simulation_stats = simulation_summary(simulations)
    metadata: dict[str, object] = {
        "data": "entirely synthetic; independent people; complete 180-day follow-up",
        "config": asdict(config), "propensity_clip": [DEFAULT_CLIP, 1 - DEFAULT_CLIP],
        "bootstrap_method": "patient-level pairs; refitted PS; 95% percentile; all scheduled fits required",
        "smd_denominator": "fixed unweighted pooled population SD",
        "python": platform.python_version(), "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__},
        "summary": result.summary, "oracle_risk_difference": oracle,
        "simulation_summary": simulation_stats,
        "source_hash_normalization": "UTF-8 text with LF line endings (as committed by .gitattributes)",
        "source_sha256": {path: hashlib.sha256((ROOT / path).read_text(encoding="utf-8").encode("utf-8")).hexdigest()
                          for path in ("src/core.py", "src/report.py", "src/run.py", "queries/build_cohort.sql", "PROTOCOL.md")},
    }
    destination.mkdir(parents=True, exist_ok=True)
    tables = {
        "cohort.csv": result.patients, "cohort_flow.csv": flow,
        "balance.csv": result.balance, "descriptives.csv": result.descriptives,
        "effects.csv": effects_frame, "bootstrap.csv": boot,
        "sensitivity.csv": sensitivity_frame, "simulation_runs.csv": simulations,
        **{f"source_{name}.csv": frame for name, frame in sources.items()},
    }
    for name, frame in tables.items():
        frame.to_csv(destination / name, index=False, float_format="%.10g", lineterminator="\n")
    figures = save_figures(destination, result, effects_frame, oracle, config.bootstrap)
    report = render_report(config, result, flow, effects_frame, sensitivity_frame, oracle, simulation_stats)
    (destination / "REPORT.md").write_text(report, encoding="utf-8", newline="\n")
    files = [*tables, *figures, "REPORT.md"]
    metadata["artifact_sha256"] = {
        name: hashlib.sha256((destination / name).read_bytes()).hexdigest() for name in files
    }
    (destination / "manifest.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n"
    )
    return metadata


def render_report(config: RunConfig, result: Estimate, flow: pd.DataFrame,
                  effects: pd.DataFrame, sensitivity: pd.DataFrame, oracle: float,
                  simulations: dict[str, float | int]) -> str:
    stats = result.summary
    weighted = effects.loc[effects.method == "IPTW"].iloc[0]
    simulation_text = "Simulation validation was not requested."
    if simulations["scheduled"]:
        simulation_text = f"Scheduled {simulations['scheduled']} consecutive seeds; {simulations['successful']} successful, {simulations['failed']} failed."
        if simulations["successful"]:
            simulation_text += (f" Crude mean error {simulations['crude_mean_error']:.4f}, RMSE {simulations['crude_rmse']:.4f};"
                                f" IPTW mean error {simulations['weighted_mean_error']:.4f}, RMSE {simulations['weighted_rmse']:.4f}.")
    return f"""# Reproducible synthetic target-trial report

**Not clinical evidence. All source rows and outcomes were generated locally.**

## Question and estimand

What is the 180-day risk difference for initiating A versus B among eligible new
users? The primary estimator is the difference of normalized, stabilized-IPTW
group risks (Hajek-style). See [the protocol](../../PROTOCOL.md) for the design.
No real treatment, institution, employer, or research participant is represented.

## Cohort and data contract

Source seed: {config.seed}; generated people: {config.n:,}; eligible people: {int(stats['n']):,}.
One row per person in all three source tables; no orphan/missing IDs, required
values, or invalid binary codes; complete follow-up required. NULL event day
means confirmed event-free follow-up, not missing ascertainment. Validation
fails rather than silently dropping broken joins or incomplete follow-up.

Sequential exclusions (a person is counted once):

{markdown_table(flow)}

Group descriptives; age SD is the sample SD, events are observed synthetic events:

{markdown_table(result.descriptives)}

## Main effect and uncertainty

Risks and differences in the table are proportions. Multiply by 100 for percentage points.

{markdown_table(effects)}

IPTW estimate: **{weighted['estimate'] * 100:.2f} percentage points**, 95% percentile
bootstrap CI **[{weighted['ci_lower'] * 100:.2f}, {weighted['ci_upper'] * 100:.2f}]**.
The known DGP benchmark for these eligible baseline covariates is
{oracle * 100:.2f} percentage points. It is computed from generating probabilities,
not observed counterfactual outcomes; it is not an exact population effect.

Bootstrap: {config.bootstrap:,} scheduled independent-person pairs resamples,
seed {config.bootstrap_seed}; propensity model, treatment prevalence, and weights
refitted for every resample. No resamples were silently discarded or redrawn.
This interval does not incorporate unmeasured confounding, measurement error,
or uncertainty about the chosen causal model; coverage is not validated here.

![Crude and IPTW risk differences with 95% bootstrap intervals; dashed line shows the known synthetic benchmark.](effect.png)

## Balance, overlap, and weights

{markdown_table(result.balance)}

SMD uses a fixed unweighted pooled population-SD denominator for both comparisons.
The 0.10 guide is a diagnostic convention, not evidence that exchangeability holds.
Raw estimated propensity range: [{stats['propensity_min_raw']:.4f}, {stats['propensity_max_raw']:.4f}].
Numerical clipping policy: [0.01, 0.99]; {int(stats['propensity_clipped_count'])} clipped scores.
Mean weight: {stats['mean_weight']:.4f}; maximum: {stats['maximum_weight']:.4f};
99th percentile: {stats['weight_p99']:.4f}; pooled ESS: {stats['effective_sample_size']:.1f}
of {int(stats['n'])} people. ESS is a weight diagnostic, not an inferential sample size.
Rank deficiency, separation, or nonconvergence stops estimation.

![Age and comorbidity SMDs before and after IPTW; shape and color distinguish estimates.](balance.png)

![Raw propensity and stabilized-weight distributions by observed treatment, shown with solid and dashed outlines.](overlap_weights.png)

Both histograms use 20 shared bins per panel and density normalization within
each treatment group; input observations are in cohort.csv. Unit: independent
synthetic person. These graphs do not establish positivity in unseen populations.

## Prespecified weight sensitivity (descriptive)

{markdown_table(sensitivity)}

Quantile winsorization caps pooled weights without removing people. It can trade
variance for bias and alter the weighted target; it is not a cure for failed
positivity. No sensitivity-analysis CIs or hypothesis tests were computed.
Do not select a policy by whether its estimate is significant.

## Repeated simulation check

{simulation_text}

Every scheduled seed and any failure is in simulation_runs.csv. Errors are
relative to the corresponding finite-cohort DGP benchmark. Monte Carlo SEs of
mean errors are recorded in manifest.json when at least two runs succeed.
This evaluates point-estimator behavior in ONE favorable, correctly specified,
measured-confounding mechanism. It is not external validation or a CI-coverage study.

## Reproduce and limitations

Configuration, library versions, code hashes, and artifact hashes are in
manifest.json. CSVs include all synthetic source rows, cohort/weights, balance,
effects, sensitivity results, bootstrap estimates, and repeated-seed outcomes.

No censoring, competing risks, adherence changes, clustering, longitudinal
treatment, unmeasured confounding, or real-world missing-data mechanism is
implemented. Day 0 is a relative initiation index, not a validated longitudinal
claims algorithm. The script does not import real records. Conclusions apply
only to this simulation. See [methods and references](../../docs/METHODS.md).
"""
