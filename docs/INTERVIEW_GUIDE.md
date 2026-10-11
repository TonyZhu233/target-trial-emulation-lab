# Explain the project without overstating it

## A concise introduction

“This is an independent synthetic target-trial demonstration connecting an
explicit protocol to a validated SQLite cohort and a Python analysis. It
estimates a 180-day risk difference using normalized IPTW, reports bootstrap
uncertainty and balance/overlap diagnostics, and checks estimator behavior
against a known simulation mechanism across repeated seeds. It is a method
and engineering portfolio project, not evidence from patients.”

## Questions to prepare

1. **Why might a crude comparison be biased?** Older/comorbid people receive A
   more often and have higher baseline risk; the mechanism makes these measured
   confounders. Weighting attempts to balance those covariates.
2. **Why specify time zero?** Eligibility, treatment assignment and outcome
   follow-up must align. This toy schema aligns them at simulated baseline;
   it does not solve real-world longitudinal immortal-time problems.
3. **Why a risk difference rather than the logistic coefficient?** A marginal
   probability contrast is interpretable in percentage points. -0.45 is a
   conditional generating log-odds coefficient, not the marginal RD.
4. **What exactly does the bootstrap refit?** The propensity model, treatment
   prevalence, weights and both risk estimates. Resampling with fixed weights
   would not reproduce the full estimation pipeline.
5. **Does SMD below 0.10 make the effect causal?** No. It describes measured
   balance under a stated SD convention; it cannot test no unmeasured confounding.
6. **Why save raw scores and clipped-score counts?** Clipping can hide extreme
   probabilities. Raw overlap and weight distributions make that policy visible.
7. **Why report ESS by group?** Pooled ESS can hide one treatment group's weight
   concentration; neither ESS is automatically the right inferential degrees of freedom.
8. **Why not force every seed to show a negative effect?** Sampling variation is
   expected. Evaluate aggregate point-estimation error, not favorable signs.
9. **What does weight winsorization change?** The bias/variance trade-off,
   covariate balance and potentially the weighted target. It is not a positivity fix.
10. **How would this need to change for real claims?** Longitudinal new-user
    definitions, validated code lists, ascertainment/missingness, time-varying
    censoring, clustering, and governance. None is claimed implemented here.

## A hands-on ownership check

Before describing the project as a personal skill, run the pipeline in a new
folder, reconcile the cohort counts with SQL, reproduce one normalized risk
by hand, explain a bootstrap draw, change a seed without selecting the sign,
and demonstrate the missing-outcome/incomplete-follow-up tests. Review the
AI-assistance disclosure and keep the synthetic-data qualifier in your CV.

Do not describe this project as real-world clinical effectiveness, production
MLE deployment, large-scale data engineering, or work completed for a lab or
employer. Those require different evidence.
