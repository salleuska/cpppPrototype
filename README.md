# cpppPrototype

The **cppp** package implements the *Calibrated Posterior p-value* (cppp) procedure described in  
*Computational Methods for Fast Bayesian Model Assessment via Calibrated Posterior p-values*.

It provides a general framework for a MCMC engine—to:

1. Compute calibrated posterior p-values (cppp),
2. Estimate their Monte Carlo variance using the idea of the **transfer effective sample size (ESS)**
3. Can handle different MCMC engines. **NIMBLE** and R for no, other MCMC engines later.

---

## Concept

Given data \(y\), a model \(p(\theta, y)\), and a discrepancy function \(D(y,\theta)\):

1. Run a long MCMC chain on \(p(\theta\mid y)\) → get \(M\) draws and compute
   \[
   \Delta_i = D(y_i^*, \theta_i) - D(y, \theta_i),
   \quad
   y_i^* \sim p(y^*\mid \theta_i).
   \]
   This $\Delta$-chain encodes the discrepancy structure of the model.

2. Generate \(r\) *calibration replicates*:
   - simulate new datasets \(\tilde y_j\) from the model,
   - run short chains of length \(\tilde m\) for \(p(\theta\mid \tilde y_j)\),
   - compute short-run posterior-predictive p-values \(\hat p_j\).

3. Combine:
   \[
   \widehat{\text{cppp}}
   = \frac{1}{r}\sum_{j=1}^r
     \mathbf{1}\{\hat p_j \le \hat p_{\text{obs}}\}.
   \]

4. Estimate the Monte Carlo variance using the *transfer ESS* idea:
   match each \(\hat p_j\) to a quantile on the observed Δ-chain,
   compute the transfer autocorrelation, and estimate the cppp variance.
  
---

## Package architecture

### Main functions 
| Function | Purpose |
|-----------|----------|
| `runCalibration()` | Generic function for calibration |
| `compute_cppp()` | Computes cppp from \(\hat p_{\text{obs}}\) and \(\hat p_j\). |
| `runCalibrationNIMBLE()` | NIMBLE-specific setup that builds the generic inputs and calls `runCalibration()`. |

### Helpers
| Function | Role |
|-----------|------|
| `transfer_ess_variance()` | Estimates cppp variance via transfer ESS. |
| `make_col_disc_fun()` | Returns a function that extracts an “online” discrepancy column. |
| `make_offline_disc_fun()` | Returns a function that computes discrepancies “offline.” |

Possible? 

| `make_data_sim_fun()` | Returns a function that simulates new datasets \(\tilde y\). |
| `make_MCMCfun()` | Returns a closure that runs one short MCMC given new data. |

### Data object

`cpppResult` (S3) 

  - cppp estimate, se, confidence interval
  - informations about the calibration procedure? number of replicates and mcmc per replicated
  - methods: `print`, `summary`, `plot`

---

## Typical workflow 

1. **Build your MCMC engine**
   - For NIMBLE: configure a model and compiled MCMC object.
   - Otherwise: supply a function `MCMC_fun(new_data, control)` that runs an MCMC and returns samples.

2. **Provide a data simulator**
   - `new_data_fun()` draws a synthetic dataset \(\tilde y\) from your model.

3. **Provide a discrepancy handler**
   - *Online*: the discrepancy or PPP is computed during MCMC; just extract the column.
   - *Offline*: compute \(D(y,\theta)\) and \(D(y^*,\theta)\) after sampling.

4. **Run calibration**
   - Call `runCalibration()` with your functions and number of replicates \(r\).
   - It loops: simulate → run short chain → compute \(\hat p_j\).

5. **Compute cppp and variance**
   - `compute_cppp(p_hat_obs, p_hat_cal)`
   - `transfer_ess_variance(delta_chain, p_hat_obs, p_hat_cal, m_tilde, c=1.3)`

---

## Algorithm 

**A. Run long chain**
1. Run long MCMC → \(\Delta\)  chain 
2. \( \hat p_{\text{obs}} = M^{-1}\sum \mathbf{1}\{\Delta_i \le 0\} \)

**B. Calibration (r short runs)**
For each replicate j:  

  - Simulate \(\tilde y_j\)
  - Run short MCMC (\(\tilde m_j\) iterations)
  - Compute \(\hat p_j = \tilde m_j^{-1}\sum \mathbf{1}\{\tilde\Delta_{j,t}\le0\}\)

**C. Estimate cppp**
\[
\widehat{\text{cppp}} = r^{-1}\sum_{j=1}^r
  \mathbf{1}\{\hat p_j \le \hat p_{\text{obs}}\}.
\]

**D. Variance via transfer ESS**
For each j:  

  1. Let \(q_j=\hat p_j\).
  2. Find threshold $q_j^{*}$ so empirical CDF of Δ at \(q_j^*\)= \(q_j\).
  3. Form indicator \(Z_i^{(j)}=\mathbf{1}\{\Delta_i\le q_j^*\}\).
  4. Estimate integrated autocorrelation \(\tau_j\), 
  5. Transfer ESS \(\widetilde{\text{ESS}}_j=\tilde m_j/\tilde\tau_j\).
  6. $\widehat{\mathrm{Var}}(\hat{cppp})=$


**E. Variance via two-piece decomposition (paper formulation)**

The cppp variance decomposes into two parts:

\[
\operatorname{Var}[\widehat{\text{cppp}}(y)]
=
\frac{1}{r}
\,\mathbb{E}_Y
\Big[
  \operatorname{Var}_{K|Y}(
    \mathbf{1}\{K \le \tilde m\,\widehat{\mathrm{ppp}}(y)\}\mid Y)
\Big]
+
\operatorname{Var}_Y
\Big[
  \mathbb{E}_{K|Y}(
    \mathbf{1}\{K \le \tilde m\,\widehat{\mathrm{ppp}}(y)\}\mid Y)
\Big].
\]

---

##   Future extensions
- Support for other MCMC frameworks (`runCalibration_other()`).

---

## Example idea

```r
# create a discrepancy extractor
disc_fun <- make_col_disc_fun("ppp_column")

# or an offline calculator
disc_fun <- make_offline_disc_fun(control = list(n_pred = 100))

# define short-chain runner
MCMC_fun <- make_MCMCfun(niter = 200, data_nodes, cmodel, cmcmc)

# define simulator for new datasets
new_data_fun <- make_data_sim_fun(...)

# orchestrate calibration
res <- runCalibration(MCMC_samples_obs, MCMC_fun, new_data_fun, disc_fun,
                      num_reps = 50, control = list(seed = 123))

# compute cppp and variance
cppp_val <- compute_cppp(res$p_hat_obs, res$p_hat_cal)
var_out  <- transfer_ess_variance(res$delta_chain, res$p_hat_obs, res$p_hat_cal,
                                  m_tilde = res$m_tilde, c = 1.3)

summary(var_out)
```
