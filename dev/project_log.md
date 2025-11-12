# cpppPrototype — Development Log  


Encountered Nimble compilation issue
  • Error:
“Invalid project argument; models and nimbleFunctions must be compiled before being used to specify a project.”
  • Logged for future debugging; not critical for today’s structural progress.
  • We decided to postpone this debugging until the simulator architecture is locked in.

## 2025-11-12 — Calibration engine finalized + discrepancy unification

✔ Consolidated discrepancy interface (disc_fun)
  • Unified offline and online discrepancies under a common interface:
  • Every disc_fun() now returns:
  • obs: discrepancies computed on the current dataset
  • sim: discrepancies computed on posterior predictive replications
  • Removed ambiguity between “discrepancy” and “disc_fun”.
  • Adopted helpers:
  • make_col_disc_fun() (online, uses columns already in MCMC samples)
  • make_offline_disc_fun() (offline, computes discrepancies from user functions)

Result:
A single, consistent contract for discrepancies across all engines (Nimble, Stan, JAGS, etc.).

⸻

✔ Updated runCalibrationNIMBLE() interface
  • Replaced the discrepancy argument with a user-supplied disc_fun.
  • Removed all internal discrepancy construction logic.
  • Ensured that Nimble backend simply passes:
  • MCMC_samples
  • new_data_fun()
  • disc_fun()
  • to the generic runCalibration() engine.

This enforces the design philosophy:
The backend handles sampling. The user controls the discrepancy.

⸻

✔ runCalibration() updated to use the unified discrepancy
  • Now checks explicitly for obs and sim fields.
  • Properly computes:
  • PPP_obs = P(sim ≥ obs) for the observed world
  • PPP_rep = PPP in each calibration world
  • CPPP = P(PPP_rep ≤ PPP_obs)
  • Added validation:
  • vector lengths match
  • structure of discrepancy return object
  • row-selection logic is stable

runCalibration() is now the correct, backend-neutral engine.

⸻

✔ Created Option B plan for Nimble posterior predictive simulator
  • Designed a fully dependency-driven Nimble simulator (make_nimble_pp_simulator()).
  • Identified its usage:
  • as inner simulator for offline discrepancy (new_data_fun inside control)
  • as outer simulator for calibration worlds
  • Confirmed this structure matches the statistical definition of CPPP.

Though not yet fully implemented inside the wrapper, we now have a clear, correct blueprint.

⸻

✔ Newcomb example: skeleton assembled
  • Built a working pipeline for:
  • Nimble model (Newcomb)
  • Offline discrepancy (min observation)
  • Inner posterior predictive simulation
  • disc_fun built via make_offline_disc_fun()
  • Ready to plug into runCalibrationNIMBLE()

This example will serve as the first testbed once the Nimble simulator is finalized.

⸻

---

## **2025-11-12 — cPPP Computation Utilities**

### ✔ Implemented `computeCppp()`
Computes:
\[
\hat p_{cPPP} = \frac{1}{r} \sum I\{\text{pppCal}_j \le \text{pppObs}\}.
\]

Details:
- Validates scalar and vector inputs
- Handles empty calibration vectors gracefully

### ✔ Added tests for `computeCppp()`
Tests confirm:
- Correct numeric results for known inputs
- Monotonicity in `pppObs`
- Bounds in `[0,1]`
- Handling of empty calibration vectors

---

## **2025-11-12 — Roxygen, Math, and Documentation Fixes**

### ✔ Fixed Rd cross-reference warnings
- `[0,1]` incorrectly parsed as a link
- Replaced with `\eqn{[0,1]}` or escaped brackets `$begin:math:display$0,1$end:math:display$`

### ✔ Standardized mathematical expressions
Replaced `$...$` with:
- `\eqn{}` for inline math  
- `\deqn{}` for display math

---

## **Earlier — S3 Class and Base Infrastructure**

### ✔ Created S3 constructor: `new_cppresults()`
- Preserved original comments  
- Added `validate_cpppResult()`
- Checks:
  - `cppp ∈ [0,1]` when finite
  - `ppp` and `obs_ppp` numeric + within bounds  

### ✔ Constructor tests added
Ensures:
- Correct class structure  
- Proper lengths and numeric type  
- Bounds enforcement  

---

## **Earlier — Test Framework and Interface Scaffolding**

### ✔ testthat infrastructure set up
- Added `tests/testthat.R`
- Global helper for reproducible seeds

### ✔ Interface presence tests
Checked existence (not behavior) of:
- `computeCppp()`
- `transferAutocorrelation()`
- `runCalibration()`

These remain green even before full implementation.

---

## **Earlier — Package Setup & Cleanup**

### ✔ Cleaned and corrected `DESCRIPTION`
- Proper `Authors@R` formatting  
- Fixed license pointer  
- Removed duplicate `Suggests` entries  
- Added `URL`, `BugReports`, and `Config/testthat/edition: 3`  
- Added required newline at end of file  

### ✔ Added minimal `LICENSE`
Two-line file as required by R: