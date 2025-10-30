# CPPP Package — Design

_Last updated: 2025‑10‑02_

## Goal
Provide a backend‑agnostic engine to compute **calibrated posterior predictive p‑values (CPPP)** using three user‑supplied functions: a data simulator, an MCMC runner, and a discrepancy extractor. 
The version 0.1 delivers a reproducible workflowsupporting plain R or NIMBLE via adapters.

---

## Scope

**To do**

  - Backend‑agnostic engine `runCalibration()`.
  - **Online discrepancies only** (the indicator is already available per MCMC draw).
  - Fixed replicate chain length `m̃` (no ESS tuning yet).
  - Serial execution (no parallel).
  - Deterministic seeding per replicate.
  - Minimal S3 result object.

**Out (later versions)**
  
  - Offline discrepancy builder (simulate predictive nodes inside `disc_fun`).
  - ESS‑based choice of `m̃` and transfer‑ESS.
  - Plugin / bootstrap SEs that account for within‑replicate MC noise.
  - Parallelization; additional adapters (Stan/JAGS/torch).

---

## Core idea

### Data shapes

- **`real_data_samples`**: `data.frame` with `n_draws × p_params`. Column names are free‑form; functions may depend on them.
- **`new_data`**: arbitrary object consumed by `MCMC_fun` (e.g., list of data nodes for NIMBLE, or a data.frame/vector in plain R).
- **`samples_df`** (from `MCMC_fun`): `data.frame` with `m̃ × q` columns. Must contain what `disc_fun` needs (e.g., a `disc_indicator` column).

### Function interfaces

- **Data simulator**  
  `new_data_fun(real_data_samples)` → `new_data`

- **MCMC runner**  
  `MCMC_fun(new_data, control = list())` → `samples_df`

- **Discrepancy extractor**  
  `disc_fun(samples_df)` → **logical/integer vector** of length `nrow(samples_df)` with values in `{0,1}` (one indicator per draw).

> **Note**: Keep these as **closures** to bake in backend specifics (e.g., compiled objects, node names). Avoid global state.

### Determinism & seeds
- `run_cppp(..., seed_base = NULL)` controls reproducibility.
- If `seed_base` is not `NULL`, replicate `j` sets `set.seed(seed_base + j)` **once** before invoking `new_data_fun` and `MCMC_fun`. Both operate under that RNG state. This makes results identical in serial vs. later parallel modes.

### Tie rule
- Define success for calibration replicate `j` as `I(ppp_rep_j ≤ ppp_obs)`.

### PPP definitions (MVP)
- `ppp_obs = mean( disc_fun(real_data_samples) )`.
- For replicate `j`, `ppp_rep_j = mean( disc_fun(samples_df_j) )`.

---

## Engine lifecycle
1. **Compute observed PPP**: `ppp_obs` from `real_data_samples` via `disc_fun`.

2. **Replicate loop** for `j = 1..r`:
   - Set RNG (if `seed_base` provided).
   - `y_tilde_j ← new_data_fun(real_data_samples)`.
   - `samples_df_j ← MCMC_fun(y_tilde_j, control)` (with `m̃ = control$niter` draws recommended).
   - `ind_j ← disc_fun(samples_df_j)` (vector of 0/1).
   - `ppp_rep_j ← mean(ind_j)`.
3. **Calibrated p‑value**: `cppp ← mean( ppp_rep_j ≤ ppp_obs )`.
4. **Return** an S3 `cppp_result` with core fields (below).

---

## Result object
Class: `cppp_result`
- `cppp` (numeric scalar)
- `r` (integer; number of replicates)
- `ppp_obs` (numeric scalar)
- `ppp_reps` (numeric vector length `r`)  
- `meta` (list): `seed_base`, `control`, optional `backend` tag (from `attr(MCMC_fun, "backend")`).

---

## Validation & diagnostics (MVP)
- **Shape checks**: `disc_fun` returns 0/1 vector; lengths match `nrow(samples_df)`; `nrow(samples_df) ≥ 1`.
- **Degeneracy warnings**: If `mean(ind_j) ∈ {0,1}` (or very near), issue a warning suggesting a larger `m̃` or a different discrepancy.
- **Minimal summaries**: Report `ppp_obs`, median/quantiles of `ppp_reps`, and `cppp`.

---

## Milestones (incremental plan)
1. **M0 – Contracts**: finalize this doc.
2. **M1 – Thin engine**: implement `runCalibration` with checks and summaries (no SE).  
3. **M2 – Toy backend**: generic R  functions 
4. **M3 – NIMBLE adapter (minimal)*
5. **M4 – SE based on transfer method: and a `print()` method.  
6. **M5 – Fit & finish**: documentation, examples, and a single vignette.

---

## Future extensions (post‑MVP)
- **Offline discrepancy builder** using predictive node simulation.
- **ESS‑aware m̃** selection; transfer‑ESS.
- **Variance estimators**: plugin; bootstrap‑normal; block bootstrap.
- **Parallelization** with per‑replicate RNG streams (L’Ecuyer).
- Additional adapters (Stan/JAGS); plotting helpers.

---

## Open questions (decide early)
- Should `new_data_fun` sample **one** θ per replicate from `real_data_samples` (default) or support multiple? _(MVP: one.)_
- Minimal `control` schema for `MCMC_fun` (recommend `list(niter = m̃, thin = 1)`; everything else baked into the closure).
- Naming convention for the online discrepancy column (e.g., `disc_indicator`).

---

## Tiny templates (illustrative; not executable)
```r
# new_data_fun(real_data_samples) -> new_data
# • May internally draw one row of real_data_samples
# • Returns whatever MCMC_fun expects (e.g., list of data nodes)

# MCMC_fun(new_data, control) -> samples_df
# • Returns data.frame with one row per draw (m̃ rows)
# • Includes a column used by disc_fun (e.g., "disc_indicator")

# disc_fun(samples_df) -> 0/1 vector of length nrow(samples_df)
# • For MVP, typically: samples_df[["disc_indicator"]]
```

