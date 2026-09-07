# Target trial protocol

This educational example emulates a simple comparative-effectiveness trial
using entirely synthetic data. It is not a clinical or causal claim.

| Component | Specification |
|---|---|
| Eligibility | Age 50–90, at least 365 days of prior enrollment, no prior use of either study treatment |
| Treatment strategies | Initiate treatment A versus initiate treatment B |
| Assignment | Observational; modeled with baseline age and comorbidity |
| Time zero | Date of first eligible initiation |
| Follow-up | 180 days from time zero |
| Outcome | Synthetic binary adverse outcome during follow-up |
| Causal contrast | Observational analogue of intention-to-treat |
| Estimand | 180-day risk difference, A minus B |
| Analysis | Logistic propensity score, stabilized IPTW, balance and effective-sample-size diagnostics |

## Assumptions

Interpretation requires consistency, positivity, conditional exchangeability,
correct measurement, and correct model specification. These assumptions are
plausible by construction only in this synthetic demonstration.

## Planned extensions

- Confidence intervals and bootstrap uncertainty.
- Weight truncation sensitivity analysis.
- Time-to-event outcomes, censoring weights, and survival curves.
- Negative-control and outcome-definition sensitivity analyses.

