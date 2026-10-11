# Version-1 verification record

Local verification date: 2026-10-11. This records engineering/simulation
checks, not a statistical certification or clinical validation.

## Executed locally

- A new isolated Windows Python 3.13.9 environment was created and dependencies
  installed from `requirements.txt`. `pip check` found no broken requirements.
- Explicit and default unittest discovery both run actual tests; the original
  zero-discovery README problem is repaired. The final suite includes source
  integrity, SQL boundaries/grain, numerical calculations, model failures,
  bootstrap refitting, CLI execution, and published-artifact provenance.
- The full showcase ran with 3,000 generated / 2,196 eligible people, 1,000
  bootstrap draws and 100 consecutive simulation seeds. All 100 simulation
  fits succeeded; no draws/seeds were selected by treatment-effect sign.
- An additional source seed 42 ran successfully with 100 bootstrap draws and
  simulations explicitly disabled. The console compatibility entry also ran.
- Generated PNGs were visually inspected for legibility, legend/label clipping,
  common scales and faithful uncertainty display. SVGs are parsed in tests.
- Blue/orange graphical contrast on white was screened at 5.185:1 / 3.867:1
  against a 3:1 reference. Their grayscale L* separation (8.214) is below the
  helper's 10-unit heuristic, so color is deliberately redundant: square/circle
  markers and solid/dashed histogram outlines. This is not a WCAG certification.
- Staged source and generated text pass whitespace checks. The demo's code/
  artifact hashes and its interval/count reconciliation are tested directly.

## Checked-in synthetic results

IPTW risk difference: -0.030565; 95% percentile bootstrap interval
[-0.054697, -0.007585]. Known finite-cohort DGP benchmark: -0.033630.
Across the scheduled 100 seeds, crude point-error RMSE was 0.024888 and IPTW
RMSE 0.009996. These values are in risk-proportion units, not percentage points.
They are not externally validated performance metrics or actual treatment effects.

## Remote verification

The [GitHub Actions workflow](https://github.com/TonyZhu233/target-trial-emulation-lab/actions/workflows/test.yml)
is configured for Windows/Linux and Python 3.12/3.13. The workflow runs all tests
and a full 1,000-bootstrap/100-seed analysis, then uploads its output. Consult
the linked live run status rather than treating configuration as proof of a
successful run. Local checks alone do not establish macOS compatibility.

## Evidence boundaries

No clinical study, real patient data, production deployment, multi-site
validation, broad model-misspecification benchmark or CI-coverage evaluation
was performed. See the protocol and methods for unsupported settings.
