# Target Trial Emulation Lab

[![test](https://github.com/TonyZhu233/target-trial-emulation-lab/actions/workflows/test.yml/badge.svg)](https://github.com/TonyZhu233/target-trial-emulation-lab/actions/workflows/test.yml)

**Protocol → validated SQL cohort → Python IPTW → uncertainty and diagnostics → reproducible report.**

An independent portfolio demonstration of causal-analysis reasoning and
reproducible analytical engineering. **Entirely synthetic data; no patient
records, restricted schemas, unpublished research, or clinical claims.**
Built with AI-assisted development; see [methods/disclosure](docs/METHODS.md).

## Inspect the result first

[Complete checked-in report](results/demo/REPORT.md) ·
[Target-trial protocol](PROTOCOL.md) · [Methods](docs/METHODS.md) ·
[Interview guide](docs/INTERVIEW_GUIDE.md)

![Age and comorbidity imbalance before weighting versus after IPTW; squares and circles distinguish estimates.](results/demo/balance.png)

![Crude and IPTW 180-day risk differences with 95% patient-level bootstrap intervals and a known synthetic benchmark.](results/demo/effect.png)

Negative differences mean lower **simulated** risk under A than B. These are
not estimates for actual treatments. Tables, overlap/weight plots, all
bootstrap draws, and source/code hashes accompany the report.

## What the project demonstrates

- **Study design:** eligibility, baseline time zero, strategies, outcome window,
  estimand and assumptions documented before the version-1 runs.
- **SQL and data quality:** one-person grain, aligned ID sets, sequential
  exclusions, inclusive day/age boundaries, explicit no-event versus missing
  ascertainment, and rejection of incomplete follow-up.
- **Estimation:** logistic propensity model with convergence/separation guards,
  stabilized weights, normalized weighted risks and a marginal risk difference.
- **Statistical judgment:** fixed-denominator balance diagnostics, empirical
  overlap, group/pooled ESS, extreme weights, refitted-propensity bootstrap CIs,
  and prespecified descriptive weight-winsorization sensitivity.
- **Verification:** hand-calculated estimators, numerical score equations,
  failure-path tests, real CLI integration, and repeated-seed point-estimator
  checks against the known data-generating mechanism.
- **Reproducibility:** offline run after installation, pinned direct dependencies,
  seeds/configuration, raw generated CSVs, transformed cohort, SVG/PNG figures,
  machine-readable summaries and SHA-256 manifests.

## Run from a clean environment

Supported/tested by CI: Python 3.12 and 3.13 on Windows and Linux. Clone or
download this repository and open a terminal **in its root directory**.
No API keys, accounts, databases, GPUs, or patient-data access are needed.

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe -m src.run
```

macOS/Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python -m src.run
```

The default analysis generates 3,000 people, uses 1,000 patient-level bootstrap
draws, and checks 100 consecutive simulation seeds. Open `results/local/REPORT.md`
after it completes. Duration depends on the machine. The analysis itself does
not make network calls; dependency installation needs a package source.

Using your selected Python environment:

```bash
# Faster smoke run; not the showcase uncertainty analysis
python -m src.run --bootstrap 100 --simulations 5 --output-dir results/smoke

# Recreate the checked-in full demonstration
python -m src.run --output-dir results/demo --overwrite

# Change the source seed; use a fresh destination
python -m src.run --seed 42 --output-dir results/seed42

python -m src.run --help
```

Output must be a direct child of `results/`. A nonempty destination is refused
unless `--overwrite` is explicit. That flag replaces generated filenames only;
it does not recursively delete the folder. Locally generated results are
ignored by Git except the deliberately published synthetic `results/demo/`.

`python src/tte_lab.py` remains a console-only compatibility example.
Both explicit test discovery above and `python -m unittest -v` discover tests.

## Output contract

| Artifact | What it makes inspectable |
|---|---|
| `source_*.csv` | All generated baseline/exposure/outcome source rows |
| `cohort_flow.csv`, `cohort.csv` | Sequential exclusions, eligible people, raw/used propensities and weights |
| `descriptives.csv`, `balance.csv` | Group sizes, events, covariate moments, risks and measured balance |
| `effects.csv`, `bootstrap.csv` | Crude/IPTW estimates, percentile CIs and every bootstrap draw |
| `sensitivity.csv` | All three planned weight policies and diagnostics |
| `simulation_runs.csv` | Every scheduled consecutive seed, point error and failure status |
| `*.png`, `*.svg`, `REPORT.md` | Captioned result, overlap, weight and balance displays |
| `manifest.json` | Configuration, versions, point-validation summary, code/data/artifact hashes |

The CI matrix runs tests plus the full analysis and uploads inspectable reports.
Direct dependency pins are not a complete transitive lock file; the manifest
records the actual versions used. Numeric equivalence is tested within
tolerance, not by demanding identical figure bytes across renderers.

## Scope and limitations

This is a **point-treatment, complete-follow-up method demonstration**, not a
production healthcare system. Time zero is simulated baseline; there is no
longitudinal claim-code algorithm, censoring model, competing-risk analysis,
adherence change, cluster handling, real-data ingestion or external validation.
Measured balance cannot prove no unmeasured confounding. Clipping or capping
weights cannot repair structural positivity violations. Repeated seeds assess
point estimation under one correctly specified mechanism; CI coverage is not
validated. No offer, deployment, employer project or clinical result is implied.

The intentionally favorable mechanism makes confounding and its adjustment
visible. It is useful evidence of a workflow one can explain and reproduce,
not a substitute for experience with actual governed healthcare data.

## Repository map

```text
PROTOCOL.md                design, estimand, assumptions, fixed analysis plan
queries/build_cohort.sql   baseline eligibility and 180-day outcome
src/core.py                source contract, generator, propensity, IPTW, bootstrap
src/report.py              simulation checks, tables, plots and manifest
src/run.py                 offline CLI
src/tte_lab.py             original console-entry compatibility
tests/                     numerical, SQL, data-quality and CLI tests
docs/                      methods, references, interview preparation
results/demo/              deliberately published synthetic demonstration
```

## License

[MIT](LICENSE). Methods and development-tool references are in
[docs/METHODS.md](docs/METHODS.md); citations are not endorsements or certifications.
