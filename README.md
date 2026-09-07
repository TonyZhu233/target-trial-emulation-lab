# Target Trial Emulation Lab

A small, reproducible demonstration of how an observational healthcare study
moves from a target-trial protocol to a SQL cohort and an adjusted estimate.

[![test](https://github.com/TonyZhu233/target-trial-emulation-lab/actions/workflows/test.yml/badge.svg)](https://github.com/TonyZhu233/target-trial-emulation-lab/actions/workflows/test.yml)

The repository uses **only synthetic data**. It contains no patient data,
restricted schemas, unpublished study results, or clinical claims.

## What it demonstrates

- An explicit target-trial specification: eligibility, strategies, time zero,
  follow-up, outcome, causal contrast, and estimand.
- A one-row-per-new-user cohort built in SQL.
- A logistic propensity-score model and stabilized inverse-probability weights.
- Covariate-balance, effective-sample-size, and extreme-weight diagnostics.
- A crude and weighted 180-day risk difference.

## Run it

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python src/tte_lab.py
python -m unittest
```

The deterministic self-check verifies that weighting improves measured
baseline balance and that the estimated treatment contrast has the expected
direction in the synthetic data-generating process.

Example deterministic result:

| Diagnostic | Before weighting | After weighting |
|---|---:|---:|
| Age absolute SMD | 0.346 | 0.0004 |
| Comorbidity absolute SMD | 0.425 | 0.0014 |

The weighted risk difference is -0.031 and the weighted effective sample size
is 2,048 of 2,196 eligible synthetic patients.

## Repository map

```text
PROTOCOL.md                 target-trial specification and assumptions
queries/build_cohort.sql    eligibility, time zero, treatment, and outcome
src/tte_lab.py              synthetic data, propensity score, IPTW, diagnostics
tests/test_tte_lab.py       minimal reproducibility check
```

## Interpretation

Treatment A is assigned more often to older, more comorbid synthetic patients,
creating measured confounding. The outcome-generating process gives treatment A
a protective effect. IPTW should reduce age and comorbidity imbalance and move
the estimate toward that synthetic effect.

This is an educational method demonstration, not evidence about any treatment.
See [PROTOCOL.md](PROTOCOL.md) for assumptions and planned extensions.

## License

[MIT](LICENSE)
