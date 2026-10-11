# Methods, data contract, and references

## 1. Synthetic mechanism

With independent random draws (NumPy Generator; seeds are in the manifest):

```text
age = clip(round(Normal(66,10)),40,90)
C ~ Bernoulli(expit(-0.8 + 0.035*(age-60)))
A ~ Bernoulli(expit(-0.25 + 0.035*(age-65) + 0.9*C))
Y ~ Bernoulli(expit(-2.7 + 0.035*(age-65) + 0.9*C - 0.45*A))
event day ~ discrete Uniform(1,...,180) if Y=1, otherwise explicit NULL
prior enrollment ~ discrete Uniform(300,...,730)
prior use ~ Bernoulli(0.08)
```

Everyone is followed for 180 days. Enrollment/prior-use draws do not introduce
additional outcome confounding. The treatment coefficient **-0.45 is on the
log-odds scale**, not a risk difference. For eligible baseline covariates,
the DGP benchmark is the average of `expit(baseline-0.45)-expit(baseline)`.
This finite-cohort expected contrast is not an observed individual causal
effect or the exact population ATE. It is available only because the mechanism
is known; it must not be substituted into analyses of real data.

## 2. Data grain and SQL

Each source table has one row per `patient_id`. `patients` holds age and
comorbidity; `index_exposures` holds baseline A/B initiation and eligibility
flags; `outcomes` holds first event day and ascertainment duration. This is a
compact person-level schema, **not** a genuine longitudinal claims database.

The validator checks uniqueness, numeric/binary values, ID set equality,
event-day validity and complete follow-up before SQL runs. An absent outcome
row is an error. An explicit NULL event day accompanied by complete follow-up
is a known no-event observation. SQL applies baseline eligibility and event
days 1–180 inclusively; source data and sequential exclusion counts are saved.
The CLI has no real-record import option.

## 3. Propensity and effect estimator

The propensity model is unpenalized logistic regression:

```text
logit(e_i) = beta_0 + beta_age * standardized_age_i + beta_C * C_i
w_i = A_i*p/e_i + (1-A_i)*(1-p)/(1-e_i)
risk_A = sum(A_i*w_i*Y_i) / sum(A_i*w_i)
risk_B = sum((1-A_i)*w_i*Y_i) / sum((1-A_i)*w_i)
RD = risk_A - risk_B
```

Here `p` is observed treatment prevalence. Stabilization constants cancel
within each normalized group risk; they still affect pooled weight diagnostics.
The solver uses a damped Newton step, stable log-likelihood arithmetic, rank/
conditioning checks, and explicit convergence/separation failures. Clipping
to [0.01,0.99] is a reported numerical policy, not a solution to poor positivity.
No outcome is used to fit treatment probabilities. Score-equation and
hand-calculated estimator tests guard against silent arithmetic errors.

## 4. Diagnostics and uncertainty

Absolute SMD = absolute weighted group-mean difference divided by the **fixed
unweighted pooled population SD**. For a binary covariate, the group variance
is `p*(1-p)`. This fixed denominator is deliberate; software packages may use
different SD conventions, so compare definitions before comparing numbers.
ESS = `sum(w)^2/sum(w^2)`; it summarizes concentration of weights rather than
the number of independent observations for confidence-interval calculations.

The patient-level pairs bootstrap resamples all cohort rows with replacement.
Every resample refits the propensity model and stabilization prevalence,
recalculates weights, and estimates crude/adjusted risks. The 2.5th/97.5th
percentiles of 1,000 draws form the showcase 95% interval. All scheduled fits
must succeed; a failure stops interval production instead of redrawing.
The CLI's minimum 100 draws is for fast verification. No CI-coverage guarantee
is claimed; interval uncertainty does not include identification bias.

Normality/variance tests for continuous two-group outcomes are not applicable
to this binary weighted-risk analysis. Relevant checks here are ascertainment,
grain, finite/converged propensity fits, overlap, weight concentration and
measured covariate balance. Favorable diagnostics do not prove exchangeability.

Weight winsorization is descriptive sensitivity, not an alternative selected
for a smaller p-value. Pooled 1/99 and 5/95 caps retain all people but may change
bias, balance and the weighted target. No interval or significance test is
computed for these sensitivity estimates.

## 5. Reproducibility and repeated-seed validation

`manifest.json` records configuration, versions, and SHA-256 hashes. Source-code
hashes normalize UTF-8 text to LF, matching committed Git blobs; artifact hashes
use exact saved bytes. Source
CSVs preserve generated originals, while `cohort.csv` adds probabilities and
weights. SVG/PNG plots have tabular companions and report captions/alt text.
Figures use common histogram bins, accessible colors plus marker/line styles,
and explicit uncertainty labels. Vector exports retain editable text.

The repeated-seed check reports mean point-estimation error, its Monte Carlo
SE, and RMSE relative to each cohort's DGP benchmark. All scheduled seeds,
including failures, are retained. This is one correctly specified mechanism,
not a broad robustness benchmark, external evaluation, or coverage study.
Floating-point results should agree within tolerances across platforms;
figure byte hashes may differ due to renderer metadata. The exact recorded
environment and source hashes identify a particular generated artifact set.

## References and development disclosure

- Austin, P. C., & Stuart, E. A. (2015). Moving towards best practice when using
  inverse probability of treatment weighting (IPTW) using the propensity score
  to estimate causal treatment effects in observational studies. *Statistics
  in Medicine*, 34, 3661–3679. [doi:10.1002/sim.6607](https://doi.org/10.1002/sim.6607).
  Provides methodological context for measured-balance/weight diagnostics and
  weighting-aware uncertainty; this repository is an independent simplified
  implementation, not a replication of the paper's study.
- Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026).
  *Scientific Agent Skills: A Library of Procedural Knowledge for Research
  Agents.* [arXiv:2609.00065](https://arxiv.org/abs/2609.00065).
  Scientific-analysis and visualization procedures informed reporting and
  reproducibility checks; this citation is not methodological validation.

AI-assisted development was used to implement and review this portfolio
demonstration. No clinical deployment, employer deliverable, patient-data
access, or publication result is claimed. The project owner should be able to
explain and reproduce the design, calculations, diagnostics, and limitations
before using it as interview evidence.
