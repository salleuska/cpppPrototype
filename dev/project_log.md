# cpppPrototype — Development Log  

---

## **2025-11-12 — draft runCalibration + started tests**

### ✔ Implemented `runCalibration()`
Built from your original version and corrected only structural issues while preserving your logic and comments.

Key fixes:
- Moved return *outside* the loop.
- `ppp <- list()` → `ppp <- numeric(num_reps)`.
- Forwarded `...` to `new_data_fun()`, `MCMC_fun()`, and `disc_fun()`.
- Added checks that `disc_fun()` returns a list with `$rep` and `$obs`.
- Safe looping using `seq_len(num_reps)`.

Result:  
A **engine-agnostic calibration loop**, ready to plug into Nimble or other engines.

---

### ✔ Added test suite for `runCalibration()`
Using toy functions:
- `new_data_fun()`
- `MCMC_fun()`
- `disc_fun()`

Tests assert:
- Output vector is numeric
- Length matches `num_reps`
- All values in \([0,1]\)
- All values finite
- Proper errors for invalid inputs or malformed discrepancy objects

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