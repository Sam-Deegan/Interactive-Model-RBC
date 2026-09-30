################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## The Real Business Cycle Model: Checks Against Whelan's RATS Program        ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   From the repo root:
##     Rscript tests/verify_against_rats.R
##   Needs shiny, bslib, ggplot2 and htmltools. About three minutes, most of
##   it the sampler, which V_08 runs and V_10 runs once more.
##
## Inputs:
##   R/model.R for V_02 to V_08; R/toolkit.R and app.R from V_09 on. The
##   arithmetic of Whelan's RATS program rbc.prg is transcribed into V_02 to
##   V_05, with his calibration on its lines 8 to 12.
##
## Outputs:
##   A pass or fail line per check and a non-zero exit status if any fails.
##
## What is checked:
##   V_02 to V_05  the solver against rbc.prg: the derived parameters, the
##                 shock loading computed his way, the structural residual,
##                 and his brute-force iteration run here against the
##                 eigenvalue solution
##   V_06, V_07    the findings Whelan states in part 7: volatile investment,
##                 output tracking technology, no positive growth
##                 autocorrelation, hours rising against Gali, no propagation
##                 from iid technology
##   V_07b, V_07c  his impulse responses read off [W7 36] and [W7 37]; the
##                 wealth-effect reading of eta; every equation the Equations
##                 panel prints, held against the solver's own path
##   V_08          the estimation stage against its own definitions, not
##                 against RATS: the scalar Kalman filter against a matrix
##                 one, the sampler's acceptance rates, recovery of a known
##                 rho from a SIMULATED sample, a tight prior returning
##                 itself, and the drawn prior densities against their exact
##                 masses
##   V_09, V_10    the app around the model: the export shapes, the equation
##                 and notation lists closing over each other, the page frame,
##                 every stage's panels through testServer, and every Save
##                 PNG handler
##
## References:
##   Whelan, K. MA Advanced Macroeconomics, part 7 (the real business cycle
##     model) and part 10 (estimating DSGE models); his RATS program rbc.prg.

#-------------------------------- Script Begin --------------------------------#

################################################################################
## V_01: Setup #################################################################
################################################################################
# Note: Options, the script's own directory, the model, and the check helper.

options(scipen = 999, digits = 7)

V_01_01_here_dir <- tryCatch({
  args_vec <- commandArgs(trailingOnly = FALSE)
  file_chr <- sub("^--file=", "", args_vec[grepl("^--file=", args_vec)])
  if (length(file_chr) == 1L && nzchar(file_chr)) {
    dirname(normalizePath(file_chr))
  } else "."
}, error = function(e) ".")

source(file.path(V_01_01_here_dir, "..", "R", "model.R"))

V_01_02_fail_int <- 0L

V_01_03_check_fn <- function(label_chr, ok_lgl, detail_chr = "") {
  message(sprintf("  [%s] %-52s %s",
                  if (isTRUE(ok_lgl)) "PASS" else "FAIL",
                  label_chr, detail_chr))
  if (!isTRUE(ok_lgl)) V_01_02_fail_int <<- V_01_02_fail_int + 1L
}

################################################################################
## V: Verification #############################################################
################################################################################
# Note: One message per section, one pass or fail line per check.

#### V_02: Whelan's Calibration ################################################
# Note: rbc.prg lines 8 to 12, and the derived parameters computed his way.

V_02_01_par_lst <- list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                        eta = 1.0, rho = 0.95)

message("\nECON42240 RBC app: verification against rbc.prg\n")
message("V_02  Whelan's calibration")

V_02_02_der_lst <- C_01_01_derived_fn(V_02_01_par_lst)

# rbc.prg: comp ishare = (alpha*delta) / (delta -1 + (beta)**(-1) )
V_02_03_ishare_num <- (V_02_01_par_lst$alpha * V_02_01_par_lst$delta) /
  (V_02_01_par_lst$delta - 1 + (V_02_01_par_lst$beta)^(-1))
# rbc.prg: comp mu = 1 - beta*(1-delta)
V_02_04_mu_num <- 1 - V_02_01_par_lst$beta * (1 - V_02_01_par_lst$delta)

V_01_03_check_fn("ishare matches rbc.prg",
                 abs(V_02_02_der_lst$s_i - V_02_03_ishare_num) < 1e-14,
                 sprintf("%.6f", V_02_03_ishare_num))
V_01_03_check_fn("mu matches rbc.prg",
                 abs(V_02_02_der_lst$mu - V_02_04_mu_num) < 1e-14,
                 sprintf("%.6f", V_02_04_mu_num))

#### V_03: The Shock Loading ###################################################
# Note: rbc.prg sets R = 0 (line 46) and D2 = 0 (line 45), so its forward
#   loop keeps only the k = 0 term, H_rats = inv(I - B C) inv(Chat) D1. The
#   app's H = inv(Chat - Bhat C) D1 is the same object, since Chat B = Bhat.

message("\nV_03  The shock loading, computed Whelan's way")

V_03_01_sol_lst <- C_01_04_solve_fn(V_02_01_par_lst)
V_03_02_mat_lst <- V_03_01_sol_lst$mat_lst
V_03_03_eye_mat <- diag(length(C_01_02_names_vec))
V_03_04_b_mat   <- solve(V_03_02_mat_lst$c_hat, V_03_02_mat_lst$b_hat)

V_03_05_h_rats_mat <- solve(V_03_03_eye_mat -
                              V_03_04_b_mat %*% V_03_01_sol_lst$c_mat) %*%
  solve(V_03_02_mat_lst$c_hat) %*% V_03_02_mat_lst$d_one

V_01_03_check_fn("H equals rbc.prg's own expression",
                 max(abs(V_03_05_h_rats_mat - V_03_01_sol_lst$h_mat)) < 1e-10,
                 sprintf("max |diff| = %.2e",
                         max(abs(V_03_05_h_rats_mat -
                                   V_03_01_sol_lst$h_mat))))

#### V_04: The Structural Equations ############################################
# Note: The solution substituted back into the system it came from.

message("\nV_04  The solution satisfies the structural system")

V_04_01_res_lst <- C_04_02_residual_fn(V_03_01_sol_lst)
V_01_03_check_fn("state residual below 1e-10",
                 V_04_01_res_lst$state_num < 1e-10,
                 sprintf("%.2e", V_04_01_res_lst$state_num))
V_01_03_check_fn("shock residual below 1e-10",
                 V_04_01_res_lst$shock_num < 1e-10,
                 sprintf("%.2e", V_04_01_res_lst$shock_num))

#### V_05: Whelan's Own Iteration ##############################################
# Note: rbc.prg iterates C <- inv(I - B C) A from C = I, which solves the
#   matrix quadratic Bhat C^2 - Chat C + Ahat = 0 and lands on whichever root
#   its arithmetic path leads to. This runs his iteration here, requires its
#   answer to be a solution, records which root it reached, and requires the
#   app's answer to be the stable one.

message("\nV_05  Whelan's brute-force iteration, run here")

V_05_01_a_mat <- solve(V_03_02_mat_lst$c_hat, V_03_02_mat_lst$a_hat)
V_05_02_c_mat <- V_03_03_eye_mat
V_05_03_iter_int <- 0L
V_05_04_crit_num <- Inf

while (V_05_04_crit_num > 1e-10 && V_05_03_iter_int < 20000L) {
  V_05_05_step_mat <- solve(V_03_03_eye_mat -
                              V_03_04_b_mat %*% V_05_02_c_mat)
  V_05_06_new_mat  <- V_05_05_step_mat %*% V_05_01_a_mat
  V_05_04_crit_num <- max(abs(V_05_06_new_mat - V_05_02_c_mat))
  V_05_02_c_mat    <- V_05_06_new_mat
  V_05_03_iter_int <- V_05_03_iter_int + 1L
}

V_05_07_iter_eigen_vec <- sort(Mod(eigen(V_05_02_c_mat,
                                         only.values = TRUE)$values),
                               decreasing = TRUE)
V_05_08_app_eigen_vec  <- sort(Mod(eigen(V_03_01_sol_lst$c_mat,
                                         only.values = TRUE)$values),
                               decreasing = TRUE)

message(sprintf("    converged after %d iterations, criterion %.2e",
                V_05_03_iter_int, V_05_04_crit_num))
message(sprintf("    iteration reaches roots : %s",
                paste(sprintf("%.4f", V_05_07_iter_eigen_vec[1:3]),
                      collapse = ", ")))
message(sprintf("    the app selects roots   : %s",
                paste(sprintf("%.4f", V_05_08_app_eigen_vec[1:3]),
                      collapse = ", ")))

V_05_09_iter_res_num <- max(abs(
  V_03_02_mat_lst$c_hat %*% V_05_02_c_mat -
    (V_03_02_mat_lst$a_hat +
       V_03_02_mat_lst$b_hat %*% V_05_02_c_mat %*% V_05_02_c_mat)))

V_01_03_check_fn("the iteration's answer solves the quadratic",
                 V_05_09_iter_res_num < 1e-08,
                 sprintf("residual %.2e", V_05_09_iter_res_num))
V_01_03_check_fn("the app's answer is the STABLE root",
                 max(V_05_08_app_eigen_vec) < 1,
                 sprintf("largest root %.4f", max(V_05_08_app_eigen_vec)))

if (max(V_05_07_iter_eigen_vec) > 1) {
  message("    NOTE: here the iteration reached the EXPLOSIVE fixed point.")
  message("          His published figures are the stable one, which is what")
  message("          the app draws. This is the reason for the departure.")
}

#### V_06: The Three Findings He States in Class ###############################
# Note: [W7 31], [W7 33] and [W7 35], with the data number of [W7 34].

message("\nV_06  The three findings Whelan states in class")

V_06_01_sim_df <- C_02_01_simulate_fn(V_03_01_sol_lst, V_02_01_par_lst,
                                      n_periods_int = 1000L)
V_06_02_mom_lst <- C_02_03_moments_fn(V_06_01_sim_df)
V_06_03_irf_df  <- C_02_02_irf_fn(V_03_01_sol_lst, n_horizon_int = 25L)

# [W7 31] RBC models can generate cycles with volatile investment
V_01_03_check_fn("[W7 31] investment far more volatile than consumption",
                 V_06_02_mom_lst$sd_rel_vec[["investment"]] >
                   3 * V_06_02_mom_lst$sd_rel_vec[["consumption"]],
                 sprintf("i/y = %.2f against c/y = %.2f",
                         V_06_02_mom_lst$sd_rel_vec[["investment"]],
                         V_06_02_mom_lst$sd_rel_vec[["consumption"]]))

# [W7 33] the cycles rely heavily on technology fluctuations
V_01_03_check_fn("[W7 33] output follows technology closely",
                 V_06_02_mom_lst$tech_corr > 0.98,
                 sprintf("correlation %.4f", V_06_02_mom_lst$tech_corr))

# [W7 35] RBCs do not generate positively autocorrelated growth.
#   Cogley and Nason's figure for the data, quoted on [W7 34], is +0.34.
V_01_03_check_fn("[W7 35] growth autocorrelation is not positive",
                 V_06_02_mom_lst$growth_acf < 0.15,
                 sprintf("model %+.3f against data +0.34",
                         V_06_02_mom_lst$growth_acf))

# [W7 36] the output response decays from impact at his rho
V_06_04_peak_int <- which.max(V_06_03_irf_df$mod_output_val) - 1L
V_01_03_check_fn("[W7 36] output response peaks on impact",
                 V_06_04_peak_int == 0L,
                 sprintf("peak at quarter %d", V_06_04_peak_int))

# [W7 38] Gali: hours fall after a positive technology shock and the model
#   says they rise; the check passes when the model gets it wrong
V_01_03_check_fn("[W7 38] model puts hours UP, against Gali",
                 V_06_03_irf_df$mod_hours_val[1L] > 0,
                 sprintf("%+.2f on impact", V_06_03_irf_df$mod_hours_val[1L]))

#### V_07: The Propagation Test ################################################
# Note: [W7 32]. The early claim was that the model would generate cycles even
#   from iid technology. Set rho to zero and read the correlation.

message("\nV_07  The propagation test at iid technology")

V_07_01_par_lst <- modifyList(V_02_01_par_lst, list(rho = 0))
V_07_02_sol_lst <- C_01_04_solve_fn(V_07_01_par_lst)
V_07_03_mom_lst <- C_02_03_moments_fn(
  C_02_01_simulate_fn(V_07_02_sol_lst, V_07_01_par_lst, n_periods_int = 1000L))

V_01_03_check_fn("[W7 32] iid technology still gives no propagation",
                 V_07_03_mom_lst$tech_corr > 0.98,
                 sprintf("correlation %.4f at rho = 0",
                         V_07_03_mom_lst$tech_corr))

#### V_07b: Whelan's Impulse Responses, and Two Prose Claims ###################
# Note: [W7 36] and [W7 37] are pictures, so the bands are read off his
#   figures: output about 2.1 on impact, investment about 9, hours about 1.7,
#   consumption about 0.45; at quarter 50 output about 0.17 against
#   technology's 0.08. The last three checks guard two claims in the app's
#   prose: eta is a wealth effect (hours on impact fall as it rises), and no
#   corner of the sliders passes all three of the Tests tab's failures.

message("\nV_07b  Whelan's impulse responses and the corrected prose")

V_07b_01_irf_df <- C_02_02_irf_fn(V_03_01_sol_lst, n_horizon_int = 51L)
V_07b_02_end_df <- V_07b_01_irf_df[51L, ]

V_01_03_check_fn("[W7 36] impact responses match his figure",
                 abs(V_07b_01_irf_df$mod_output_val[1L] - 2.13) < 0.05 &&
                   abs(V_07b_01_irf_df$mod_investment_val[1L] - 8.96) < 0.10 &&
                   abs(V_07b_01_irf_df$mod_hours_val[1L] - 1.70) < 0.05 &&
                   abs(V_07b_01_irf_df$mod_consumption_val[1L] - 0.43) < 0.05,
                 sprintf("y %.2f, i %.2f, n %.2f, c %.2f",
                         V_07b_01_irf_df$mod_output_val[1L],
                         V_07b_01_irf_df$mod_investment_val[1L],
                         V_07b_01_irf_df$mod_hours_val[1L],
                         V_07b_01_irf_df$mod_consumption_val[1L]))

V_01_03_check_fn("[W7 36] investment and hours end below zero by q50",
                 V_07b_02_end_df$mod_investment_val < 0 &&
                   V_07b_02_end_df$mod_hours_val < 0,
                 sprintf("i %.2f, n %.2f at quarter 50",
                         V_07b_02_end_df$mod_investment_val,
                         V_07b_02_end_df$mod_hours_val))

V_01_03_check_fn("[W7 37] output ends above technology at q50",
                 V_07b_02_end_df$mod_output_val >
                   V_07b_02_end_df$mod_technology_val &&
                   abs(V_07b_02_end_df$mod_technology_val - 0.95^50) < 1e-12 &&
                   abs(V_07b_02_end_df$mod_output_val - 0.16) < 0.03,
                 sprintf("y %.3f against a %.3f",
                         V_07b_02_end_df$mod_output_val,
                         V_07b_02_end_df$mod_technology_val))

V_07b_03_eta_vec <- c(0.5, 1, 2, 4)
V_07b_04_hours_vec <- vapply(V_07b_03_eta_vec, function(e_num) {
  s_lst <- C_01_04_solve_fn(modifyList(V_02_01_par_lst, list(eta = e_num)))
  C_02_02_irf_fn(s_lst, 2L)$mod_hours_val[1L]
}, numeric(1))

V_01_03_check_fn("eta is a wealth effect: impact hours fall as it rises",
                 all(diff(V_07b_04_hours_vec) < 0),
                 paste(sprintf("%.2f", V_07b_04_hours_vec), collapse = " > "))

V_07b_05_sol_lst <- C_01_04_solve_fn(
  modifyList(V_02_01_par_lst, list(eta = 4, rho = 0.99)))
V_07b_06_hours_num <- C_02_02_irf_fn(V_07b_05_sol_lst, 2L)$mod_hours_val[1L]

V_01_03_check_fn("at eta = 4 and rho = 0.99 hours fall on impact",
                 V_07b_06_hours_num < 0,
                 sprintf("%+.2f on impact", V_07b_06_hours_num))

# the Tests tab's three failing tests at its own thresholds, over the
# corners and midpoints of the sliders' ranges
local({
  grid_df <- expand.grid(alpha = c(0.20, 1 / 3, 0.45),
                         beta  = c(0.95, 0.99, 0.999),
                         delta = c(0.005, 0.015, 0.040),
                         eta   = c(0.25, 1, 4),
                         rho   = c(0, 0.5, 0.9, 0.95, 0.99))
  all_lgl <- vapply(seq_len(nrow(grid_df)), function(i) {
    p_lst <- as.list(grid_df[i, ])
    s_lst <- C_01_04_solve_fn(p_lst)
    if (!isTRUE(s_lst$converged_lgl)) return(FALSE)
    m_lst <- C_02_03_moments_fn(C_02_01_simulate_fn(s_lst, p_lst, 400L))
    h_num <- C_02_02_irf_fn(s_lst, 2L)$mod_hours_val[1L]
    m_lst$growth_acf > 0.15 && m_lst$tech_corr < 0.95 && h_num < 0
  }, logical(1))
  V_01_03_check_fn("no corner of the sliders passes all three failures",
                   !any(all_lgl),
                   sprintf("%d settings swept", nrow(grid_df)))
})

#### V_07c: The Equations Panel, Held Against the Solver #######################
# Note: Every model equation the Equations panel prints, as its 2.1b and 2.4a
#   versions write it, evaluated on the solver's impulse responses at four
#   calibrations. After impact the path is deterministic, so the Euler
#   equation holds with the expectation dropped. The 2.4a decision rules and
#   the 2.1c innovation standard deviation are checked the same way.

message("\nV_07c  The Equations panel, equation by equation, against the ",
        "solver")

local({
  cal_lst <- list(
    list(alpha = 1 / 3, beta = 0.99,  delta = 0.015, eta = 1.0, rho = 0.95),
    list(alpha = 0.40,  beta = 0.97,  delta = 0.030, eta = 2.5, rho = 0.80),
    list(alpha = 0.25,  beta = 0.995, delta = 0.010, eta = 0.5, rho = 0.60),
    list(alpha = 0.36,  beta = 0.98,  delta = 0.025, eta = 4.0, rho = 0.99))
  worst_vec <- c()
  for (p in cal_lst) {
    s_lst <- C_01_04_solve_fn(p)
    z_df  <- C_02_02_irf_fn(s_lst, 60L)
    y <- z_df$mod_output_val; c <- z_df$mod_consumption_val
    i <- z_df$mod_investment_val; k <- z_df$mod_capital_val
    n <- z_df$mod_hours_val; r <- z_df$mod_return_val
    a <- z_df$mod_technology_val
    k_lag <- c(0, head(k, -1L)); a_lag <- c(0, head(a, -1L))
    shr <- p$alpha * p$delta / (1 / p$beta + p$delta - 1)     # panel 2.1b
    res_lst <- list(
      resource   = y - ((1 - shr) * c + shr * i),
      capital    = k - ((1 - p$delta) * k_lag + p$delta * i),
      production = y - (a + p$alpha * k_lag + (1 - p$alpha) * n),
      labour     = n - (y - p$eta * c),
      euler      = head(c, -1L) - (c[-1L] - (1 / p$eta) * r[-1L]),
      return     = r - (1 - p$beta * (1 - p$delta)) * (y - k_lag),
      ar1        = (a - p$rho * a_lag)[-1L])
    # the decision rules of [W10 11], in the solver's own coefficients
    sp <- C_06_03_space_fn(p)
    res_lst$rule_k <- k - (sp[["a_kk"]] * k_lag + sp[["a_kz"]] * a)
    res_lst$rule_c <- c - (sp[["a_ck"]] * k_lag + sp[["a_cz"]] * a)
    worst_vec <- pmax(if (length(worst_vec)) worst_vec else 0,
                      vapply(res_lst, function(v) max(abs(v)), 0))
  }
  names(worst_vec) <- c("resource", "capital", "production", "labour",
                        "euler", "return", "ar1", "rule_k", "rule_c")
  V_01_03_check_fn(
    "every panel equation holds on the solver's path",
    all(worst_vec < 1e-10),
    if (all(worst_vec < 1e-10)) {
      sprintf("9 equations x 4 calibrations, max %.1e", max(worst_vec))
    } else
      paste(names(worst_vec)[worst_vec >= 1e-10], collapse = ", "))

  # the panel's steady-state descriptors are the solver's
  der <- C_01_01_derived_fn(cal_lst[[2]])
  p2  <- cal_lst[[2]]
  V_01_03_check_fn(
    "s_i and mu in the panel equal the solver's",
    abs(der$s_i - p2$alpha * p2$delta / (p2$delta - 1 + 1 / p2$beta)) < 1e-15 &&
      abs(der$mu - (1 - p2$beta * (1 - p2$delta))) < 1e-15,
    sprintf("s_i %.4f, mu %.4f", der$s_i, der$mu))

  # sigma_eps = sqrt(1 - rho^2), from the simulator's own innovations
  p3  <- modifyList(cal_lst[[1]], list(rho = 0.6))
  sim <- C_02_01_simulate_fn(C_01_04_solve_fn(p3), p3, 20000L)
  sd_num <- stats::sd(sim$mod_innovation_val[-1L])
  V_01_03_check_fn(
    "the innovation sd is sqrt(1 - rho^2), as printed",
    abs(sd_num / sqrt(1 - 0.6^2) - 1) < 0.02,
    sprintf("%.4f against %.4f", sd_num, sqrt(1 - 0.6^2)))
})

#### V_08: The Estimation Stage ################################################
# Note: Not against RATS: rbc.prg does not estimate. The sample every check
#   uses is SIMULATED from this model at a known persistence, so the sampler
#   can be held against a truth; nothing here is an estimate of a real
#   economy. The likelihood first, then the sampler, then the drawn
#   densities. The settings mirror B_01_04 in app.R and are written out so
#   this section needs only R/model.R.

message("\nV_08  Lecture 2.4: the state space, the filter and the sampler")

V_08_01_default_lst <- list(prior_mean = 0.70, prior_sd = 0.15,
                            step = 0.05, n_obs = 120L,
                            n_draw = 800L, warm_up = 200L,
                            me_frac = 0.50, sigma_z = 1)

V_08_02_obs_lst <- C_06_04_sample_fn(
  V_02_01_par_lst,
  n_obs_int   = V_08_01_default_lst$n_obs,
  me_frac_num = V_08_01_default_lst$me_frac,
  sigma_z_num = V_08_01_default_lst$sigma_z)

###### V_08a: The Filter, Against a Plain Matrix Implementation ################
# Note: C_06_05 writes the Kalman recursion in scalars; this rebuilds it in
#   matrices, textbook form, and requires the two to agree.

V_08_03_matrix_loglik_fn <- function(space_vec, obs_lst, sigma2_num) {

  f_mat <- matrix(c(space_vec[["a_kk"]], 0,
                    space_vec[["a_kz"]], space_vec[["rho"]]), 2L, 2L)
  h_mat <- matrix(c(1, space_vec[["a_ck"]],
                    0, space_vec[["a_cz"]]), 2L, 2L)
  q_mat <- matrix(c(0, 0, 0, sigma2_num), 2L, 2L)
  r_mat <- diag(c(obs_lst$var_k_num, obs_lst$var_c_num))

  p_mat <- matrix(solve(diag(4L) - kronecker(f_mat, f_mat),
                        as.vector(q_mat)), 2L, 2L)
  s_vec  <- matrix(0, 2L, 1L)
  ll_num <- 0
  n_obs  <- obs_lst$n_obs_int

  for (t in seq_len(n_obs)) {
    s_pred <- f_mat %*% s_vec
    p_pred <- f_mat %*% p_mat %*% t(f_mat) + q_mat
    v_vec  <- matrix(c(obs_lst$k_obs_vec[t], obs_lst$c_obs_vec[t]), 2L, 1L) -
      h_mat %*% s_pred
    o_mat  <- h_mat %*% p_pred %*% t(h_mat) + r_mat
    ll_num <- ll_num - 0.5 * (log(det(o_mat)) +
                                as.numeric(t(v_vec) %*% solve(o_mat, v_vec)))
    k_mat  <- p_pred %*% t(h_mat) %*% solve(o_mat)
    s_vec  <- s_pred + k_mat %*% v_vec
    p_mat  <- p_pred - k_mat %*% h_mat %*% p_pred
  }
  ll_num - n_obs * log(2 * pi)
}

V_08_04_diff_vec <- vapply(c(0.80, 0.95, 0.99), function(rho_num) {
  par_lst   <- modifyList(V_02_01_par_lst, list(rho = rho_num))
  space_vec <- C_06_03_space_fn(par_lst)
  sigma2    <- (V_08_01_default_lst$sigma_z^2) * (1 - rho_num^2)
  abs(C_06_05_loglik_fn(space_vec, V_08_02_obs_lst, sigma2) -
        V_08_03_matrix_loglik_fn(space_vec, V_08_02_obs_lst, sigma2))
}, numeric(1))

V_01_03_check_fn("the scalar filter equals the matrix filter",
                 max(V_08_04_diff_vec) < 1e-08,
                 sprintf("max |diff| = %.2e", max(V_08_04_diff_vec)))

# the Riccati recursion carries no data, so freezing it is exact;
# tol_num = -1 never freezes
V_08_05_space_vec <- C_06_03_space_fn(V_02_01_par_lst)
V_08_06_sigma2    <- 1 - V_02_01_par_lst$rho^2
V_08_07_gap_num   <- abs(
  C_06_05_loglik_fn(V_08_05_space_vec, V_08_02_obs_lst, V_08_06_sigma2) -
    C_06_05_loglik_fn(V_08_05_space_vec, V_08_02_obs_lst, V_08_06_sigma2,
                      tol_num = -1))

V_01_03_check_fn("freezing the Riccati changes nothing",
                 V_08_07_gap_num < 1e-09,
                 sprintf("|diff| = %.2e", V_08_07_gap_num))

# a_ck and a_cz are the consumption policy the diagnostics tab prints
V_08_08_sol_lst <- C_01_04_solve_fn(V_02_01_par_lst)
V_01_03_check_fn("the state space is read off the solver",
                 abs(V_08_05_space_vec[["a_ck"]] -
                       V_08_08_sol_lst$policy_k) < 1e-14 &&
                   abs(V_08_05_space_vec[["a_cz"]] -
                         V_08_08_sol_lst$policy_a) < 1e-14,
                 sprintf("a_ck = %.6f, a_cz = %.6f",
                         V_08_05_space_vec[["a_ck"]],
                         V_08_05_space_vec[["a_cz"]]))

###### V_08b: The Sampler's Acceptance Rate ####################################
# Note: Both failure modes, almost everything rejected and almost everything
#   accepted, are reachable on the slider, and the default sits between them.

V_08_09_prior_lst <- C_06_01_beta_fn(V_08_01_default_lst$prior_mean,
                                     V_08_01_default_lst$prior_sd)

V_08_10_run_fn <- function(step_num, obs_lst = V_08_02_obs_lst,
                           prior_lst = V_08_09_prior_lst,
                           par_lst = V_02_01_par_lst,
                           n_draw_int = V_08_01_default_lst$n_draw,
                           warm_up_int = V_08_01_default_lst$warm_up) {
  C_07_03_chains_fn(par_lst, obs_lst, prior_lst$a_num, prior_lst$b_num,
                    step_num = step_num, n_draw_int = n_draw_int,
                    warm_up_int = warm_up_int)
}

V_08_11_base_lst  <- V_08_10_run_fn(V_08_01_default_lst$step)
V_08_12_coarse_lst <- V_08_10_run_fn(0.40)
V_08_13_fine_lst   <- V_08_10_run_fn(0.005)

V_01_03_check_fn("acceptance rate at the defaults is in band",
                 all(V_08_11_base_lst$accept_vec > 0.15) &&
                   all(V_08_11_base_lst$accept_vec < 0.60),
                 sprintf("%.2f and %.2f",
                         V_08_11_base_lst$accept_vec[1L],
                         V_08_11_base_lst$accept_vec[2L]))

V_01_03_check_fn("a step ten times too large stalls the chain",
                 all(V_08_12_coarse_lst$accept_vec < 0.12),
                 sprintf("%.2f at step 0.40",
                         V_08_12_coarse_lst$accept_vec[1L]))

V_01_03_check_fn("a step ten times too small accepts almost all",
                 all(V_08_13_fine_lst$accept_vec > 0.80),
                 sprintf("%.2f at step 0.005",
                         V_08_13_fine_lst$accept_vec[1L]))

# the crawling chain: a high acceptance rate and two chains that never meet
V_08_14_split_num <- abs(
  mean(V_08_13_fine_lst$draw_df$mod_rho_val[
    V_08_13_fine_lst$draw_df$mod_chain_cat ==
      unique(V_08_13_fine_lst$draw_df$mod_chain_cat)[1L]]) -
    mean(V_08_13_fine_lst$draw_df$mod_rho_val[
      V_08_13_fine_lst$draw_df$mod_chain_cat ==
        unique(V_08_13_fine_lst$draw_df$mod_chain_cat)[2L]]))

V_01_03_check_fn("the crawling chains do NOT meet",
                 V_08_14_split_num > 0.20,
                 sprintf("chains %.2f apart", V_08_14_split_num))

###### V_08c: Recovering a Known Parameter #####################################
# Note: With a flat prior the posterior is the likelihood, and with a long
#   sample it should sit on the persistence the sample was drawn at.

V_08_15_long_lst  <- C_06_04_sample_fn(V_02_01_par_lst, n_obs_int = 400L,
                                       me_frac_num = 0.50)
V_08_16_flat_lst  <- list(a_num = 1, b_num = 1)
V_08_17_recov_lst <- C_07_03_chains_fn(
  V_02_01_par_lst, V_08_15_long_lst, 1, 1,
  step_num = 0.03, n_draw_int = 1200L, warm_up_int = 400L)

V_01_03_check_fn("a flat prior and a long sample recover rho",
                 abs(V_08_17_recov_lst$post_mean_num -
                       V_02_01_par_lst$rho) < 0.03,
                 sprintf("%.4f against a truth of %.2f",
                         V_08_17_recov_lst$post_mean_num,
                         V_02_01_par_lst$rho))

V_01_03_check_fn("and the two chains agree with each other",
                 diff(range(V_08_17_recov_lst$accept_vec)) < 0.15 &&
                   V_08_17_recov_lst$post_sd_num < 0.05,
                 sprintf("posterior sd %.4f",
                         V_08_17_recov_lst$post_sd_num))

###### V_08d: A Posterior That Equals Its Prior ################################
# Note: A prior tight enough that the data stop mattering returns its own
#   mean and width as the posterior.

V_08_18_tight_lst <- C_06_01_beta_fn(0.60, 0.005)
V_08_19_stuck_lst <- C_07_03_chains_fn(
  V_02_01_par_lst, V_08_02_obs_lst,
  V_08_18_tight_lst$a_num, V_08_18_tight_lst$b_num,
  step_num = 0.008, n_draw_int = 1200L, warm_up_int = 400L)

V_01_03_check_fn("a very tight prior returns its own mean",
                 abs(V_08_19_stuck_lst$post_mean_num -
                       V_08_18_tight_lst$mean_num) < 0.005,
                 sprintf("posterior %.4f against prior %.4f",
                         V_08_19_stuck_lst$post_mean_num,
                         V_08_18_tight_lst$mean_num))

V_01_03_check_fn("and its own width",
                 abs(V_08_19_stuck_lst$post_sd_num /
                       V_08_18_tight_lst$sd_num - 1) < 0.20,
                 sprintf("posterior sd %.5f against prior sd %.5f",
                         V_08_19_stuck_lst$post_sd_num,
                         V_08_18_tight_lst$sd_num))

# the prior asked for by mean and width is the prior returned
V_08_20_round_lst <- C_06_01_beta_fn(0.70, 0.15)
V_01_03_check_fn("prior mean and width survive the round trip",
                 abs(V_08_20_round_lst$mean_num - 0.70) < 1e-12 &&
                   abs(V_08_20_round_lst$sd_num - 0.15) < 1e-12,
                 sprintf("Beta(%.3f, %.3f)", V_08_20_round_lst$a_num,
                         V_08_20_round_lst$b_num))

###### V_08e: The Drawn Prior Densities ########################################
# Note: The trapezoid rule over the points C_06_02 hands the drawing layer,
#   against the exact probability on the range drawn. The gamma panel stops
#   at 2.0, so it holds 1 - e^-8 (1 + 8) of its mass; the full-support
#   integral is checked separately.

V_08_21_prior_df <- C_06_02_priors_fn()

V_08_22_area_fn <- function(curve_chr) {
  one_df <- V_08_21_prior_df[
    V_08_21_prior_df$mod_curve_cat == curve_chr &
      V_08_21_prior_df$mod_inside_lgl, , drop = FALSE]
  n_row  <- nrow(one_df)
  sum((one_df$mod_density_val[-1L] + one_df$mod_density_val[-n_row]) / 2 *
        diff(one_df$mod_value_val))
}

V_08_23_exact_vec <- c(
  "Normal(0, 1)" = stats::pnorm(4) - stats::pnorm(-4),
  "Gamma(2, 4)"  = stats::pgamma(2, shape = 2, rate = 4),
  "Beta(2, 2)"   = 1,
  "Beta(5, 2)"   = 1)

for (V_08_24_name_chr in names(V_08_23_exact_vec)) {
  V_08_25_gap_num <- abs(V_08_22_area_fn(V_08_24_name_chr) -
                           V_08_23_exact_vec[[V_08_24_name_chr]])
  V_01_03_check_fn(
    sprintf("%s integrates to its drawn mass", V_08_24_name_chr),
    # 5e-04 is the trapezoid rule's error on the 601-point grid
    V_08_25_gap_num < 5e-04,
    sprintf("%.6f, exact %.6f", V_08_22_area_fn(V_08_24_name_chr),
            V_08_23_exact_vec[[V_08_24_name_chr]]))
}

# over the full support; the gamma runs to sixty
V_08_26_full_fn <- function(dens_fn, lo_num, hi_num, n_int = 200001L) {
  x_vec <- seq(lo_num, hi_num, length.out = n_int)
  y_vec <- dens_fn(x_vec)
  sum((y_vec[-1L] + y_vec[-n_int]) / 2) * (x_vec[2L] - x_vec[1L])
}

V_08_27_beta_num  <- V_08_26_full_fn(
  function(x) stats::dbeta(x, 2, 2), 0, 1)
V_08_28_beta5_num <- V_08_26_full_fn(
  function(x) stats::dbeta(x, 5, 2), 0, 1)
V_08_29_gamma_num <- V_08_26_full_fn(
  function(x) stats::dgamma(x, shape = 2, rate = 4), 0, 60)

V_01_03_check_fn("beta and gamma integrate to one over their support",
                 max(abs(c(V_08_27_beta_num, V_08_28_beta5_num,
                           V_08_29_gamma_num) - 1)) < 1e-06,
                 sprintf("%.8f, %.8f, %.8f", V_08_27_beta_num,
                         V_08_28_beta5_num, V_08_29_gamma_num))

################################################################################
## V_09: The Panels, the Shapes and the Downloads ##############################
################################################################################
# Note: The app around the model: the equation and notation lists, the
#   aspect-ratio floor, the figures and the page frame. Loads app.R, with
#   the toolkit sourced first so app.R's guards skip. The locale is set
#   because real PNGs are written and a Greek letter needs UTF-8.

###### V_09_01: Load the App ###################################################
# Note: shiny and ggplot2, the locale, the toolkit and app.R.

suppressMessages({
  library(shiny)
  library(ggplot2)
})

invisible(Sys.setlocale("LC_CTYPE", "C.utf8"))

source(file.path(V_01_01_here_dir, "..", "R", "toolkit.R"))
source(file.path(V_01_01_here_dir, "..", "app.R"))

message("")
message("V_09  The equations, the notation and the figures")

###### V_09_02: Read a PNG's Size ##############################################
# Note: Width and height from the IHDR chunk, bytes 17 to 24, as two
#   big-endian four-byte integers.

V_09_02_pngdim_fn <- function(path_chr) {
  con <- file(path_chr, "rb")
  on.exit(close(con))
  byte_vec <- readBin(con, "raw", n = 24L)
  if (length(byte_vec) < 24L) return(c(width = NA_real_, height = NA_real_))
  int_fn <- function(i) sum(as.integer(byte_vec[i:(i + 3L)]) * 256^(3:0))
  c(width = int_fn(17L), height = int_fn(21L))
}

###### V_09_03: Every Declared Shape Is at Least 3:2 ###########################
# Note: The constants themselves, before anything is drawn from them.

local({
  ratio_vec <- vapply(B_03_02_shape_lst,
                      function(s) s$width_px / s$height_px, 0)
  V_01_03_check_fn(
    "every declared export shape is at least 3:2",
    all(ratio_vec >= B_03_03_min_ratio_num),
    paste(sprintf("%s %.2f:1", names(ratio_vec), ratio_vec),
          collapse = ", "))
})

###### V_09_04: Every Equation Item Is Complete ################################
# Note: A label, at least one version, and a note for every version.

local({
  bad_lst <- Filter(function(one_lst) {
    is.null(one_lst$label) || !nzchar(one_lst$label) ||
      length(one_lst$versions) == 0L ||
      !all(names(one_lst$versions) %in% names(one_lst$notes)) ||
      !all(vapply(one_lst$notes, nzchar, TRUE))
  }, B_04_01_equation_lst)
  V_01_03_check_fn(
    "every equation has a label, a version and a note per version",
    length(bad_lst) == 0L,
    sprintf("%d equation item(s) checked", length(B_04_01_equation_lst)))
})

###### V_09_05: Every Version Is Keyed by a Real Stage #########################
# Note: A version keyed to a stage the radio does not carry would hide the
#   item, and a group the titles do not name is never rendered.

local({
  key_vec <- unlist(lapply(B_04_01_equation_lst,
                           function(one_lst) names(one_lst$versions)))
  grp_vec <- vapply(B_04_01_equation_lst, function(one_lst) one_lst$group, "")
  V_01_03_check_fn(
    "every version is keyed by a stage and every group is titled",
    all(!is.na(B_02_01a_rank_fn(key_vec))) &&
      all(grp_vec %in% names(B_04_02_groups_vec)),
    sprintf("stages %s, groups %s",
            paste(sort(unique(key_vec)), collapse = " "),
            paste(sort(unique(grp_vec)), collapse = " ")))
})

###### V_09_06: Every Notation Entry Is Complete ###############################
# Note: grp is one of the four groups, from is a real stage, and the gloss is
#   lower case with no full stop.

local({
  ok_grp <- vapply(B_04_03_notation_lst,
                   function(x) x$grp %in% c("var", "par", "tgt", "shk"), TRUE)
  ok_frm <- vapply(B_04_03_notation_lst,
                   function(x) !is.na(B_02_01a_rank_fn(x$from)), TRUE)
  ok_txt <- vapply(B_04_03_notation_lst, function(x) {
    nzchar(x$txt) && !grepl("\\.$", x$txt) &&
      substr(x$txt, 1L, 1L) == tolower(substr(x$txt, 1L, 1L))
  }, TRUE)
  V_01_03_check_fn(
    "every notation entry has a real group, stage and gloss",
    all(ok_grp) && all(ok_frm) && all(ok_txt),
    sprintf("%d symbol(s)", length(B_04_03_notation_lst)))
})

###### V_09_07: The Head of a Symbol ###########################################
# Note: A symbol without its argument list, subscript and superscript:
#   p(\theta) is p, a_{kk} is a, \theta^{(s)} is \theta. The back-check
#   works in head space, since k_t and k_{t-1} are the same letter.

V_09_07_head_fn <- function(sym_chr) {
  out_chr <- sub("\\(.*$", "", sym_chr)
  out_chr <- gsub("\\^\\{[^}]*\\}|_\\{[^}]*\\}", "", out_chr)
  out_chr <- gsub("\\^[A-Za-z0-9*]|_[A-Za-z0-9*]", "", out_chr)
  out_chr
}

###### V_09_08: The LaTeX That Is Not a Symbol #################################
# Note: Structure rather than notation, removed longest first together with
#   the symbol heads, so \mathcal goes before c and \lambda before a.

V_09_08_ignore_vec <- c(
  "\\begin{pmatrix}", "\\end{pmatrix}", "\\mathrm{constant}",
  "\\mathrm{Beta}", "\\mathcal{N}", "\\displaystyle", "\\underbrace",
  "\\implies", "\\approx", "\\propto", "\\infty", "\\quad", "\\tfrac",
  "\\frac", "\\sqrt", "\\left", "\\right", "\\lim", "\\log", "\\min",
  "\\max", "\\sum", "\\sim", "\\big", "\\Big", "\\mid", "\\to",
  "\\,", "\\;", "\\ ", "\\\\",
  # Euler's number in the log-linear approximation, removed last
  "e"
)

###### V_09_09: The Two Lists Close Over Each Other ############################
# Note: Forward, on the raw LaTeX: every symbol in the key appears in an
#   equation exactly as the key writes it. Back, in head space: strip every
#   subscript and superscript, take out the heads and the structural tokens,
#   and nothing may be left.

local({
  tex_vec <- unlist(lapply(B_04_01_equation_lst,
                           function(one_lst) unlist(one_lst$versions)))
  sym_vec <- vapply(B_04_03_notation_lst, function(x) x$sym, "")

  used_lgl <- vapply(sym_vec,
                     function(s) any(grepl(s, tex_vec, fixed = TRUE)), TRUE)
  V_01_03_check_fn(
    "every symbol in the notation key is used in an equation",
    all(used_lgl),
    if (all(used_lgl)) sprintf("%d symbol(s)", length(sym_vec)) else
      paste("unused:", paste(sym_vec[!used_lgl], collapse = " ")))

  flat_vec <- gsub("\\^\\{[^}]*\\}|_\\{[^}]*\\}", "", tex_vec)
  flat_vec <- gsub("\\^[A-Za-z0-9*]|_[A-Za-z0-9*]", "", flat_vec)
  head_vec <- unique(vapply(sym_vec, V_09_07_head_fn, ""))

  token_vec <- c(head_vec, V_09_08_ignore_vec)
  token_vec <- token_vec[order(nchar(token_vec), decreasing = TRUE)]
  left_vec  <- flat_vec
  for (tk_chr in token_vec) {
    left_vec <- gsub(tk_chr, "", left_vec, fixed = TRUE)
  }
  left_vec <- gsub("[0-9[:space:][:punct:]\\\\]", "", left_vec)
  V_01_03_check_fn(
    "every symbol used in an equation is in the notation key",
    all(!nzchar(left_vec)),
    if (all(!nzchar(left_vec))) "nothing left over" else
      paste("left over:",
            paste(unique(left_vec[nzchar(left_vec)]), collapse = " | ")))
})

###### V_09_10: The Decks' Own Notation ########################################
# Note: The letters the four lecture decks use for the objects the app draws,
#   from their Notation frames. The app's key may carry more, not less, and
#   may not spell one of these differently.

V_09_10_deck_vec <- c(
  "y_t", "\\rho", "\\epsilon_t", "h",                       # 1.1
  "x_t", "a", "E_t", "b_t", "N", "k", "\\lambda",           # 1.5
  "c_t", "i_t", "k_t", "n_t", "a_t", "r_t",                 # 2.1
  "\\alpha", "\\beta", "\\delta", "\\eta",
  "z_t", "a_{kk}", "a_{kz}", "a_{ck}", "a_{cz}",            # 2.4
  "k^{*}_t", "c^{*}_t", "u^{k}_t", "u^{c}_t",
  "\\theta", "\\theta^{*}", "\\theta^{(s)}", "\\nu")

local({
  have_vec <- vapply(B_04_03_notation_lst, function(x) x$sym, "")
  miss_vec <- setdiff(V_09_10_deck_vec, have_vec)
  V_01_03_check_fn(
    "every symbol on the four decks' notation frames is in the app",
    length(miss_vec) == 0L,
    if (length(miss_vec) == 0L) {
      sprintf("%d deck symbols, %d app symbols", length(V_09_10_deck_vec),
              length(have_vec))
    } else {
      paste("missing:", paste(miss_vec, collapse = " "))
    })
})

###### V_09_11: Every Figure Has a File Name and a Shape #######################
# Note: The download machinery is keyed by plotOutput id, so the two lists
#   have to be the same set.

local({
  file_vec  <- sort(names(B_03_05_figfile_lst))
  shape_vec <- vapply(B_03_05_figfile_lst, function(x) x$shape_chr, "")
  name_vec  <- vapply(B_03_05_figfile_lst, function(x) x$name_chr, "")
  V_01_03_check_fn(
    "every figure has a file name, a shape and a lower-case slug",
    all(shape_vec %in% names(B_03_02_shape_lst)) &&
      all(grepl("^[a-z]+(-[a-z]+)*$", name_vec)) &&
      !anyDuplicated(name_vec),
    sprintf("%d figure(s)", length(file_vec)))
})

###### V_09_12: Every Figure's Colours Follow the Decks' Order #################
# Note: Every builder is rendered at the defaults and its line layers read
#   back from ggplot_build: every line colour is a series token, no two lines
#   in a panel share one, the first line is blue, and no figure keeps a
#   legend. Ghost layers are skipped. See CONVENTIONS.md 6 and 7.

local({
  grDevices::pdf(NULL)   # ggplotGrob needs a device
  on.exit(grDevices::dev.off(), add = TRUE)
  par_lst <- B_01_01_default_lst
  sol_lst <- C_01_04_solve_fn(par_lst)
  sim_df  <- C_02_01_simulate_fn(sol_lst, par_lst, 400L)
  irf_df  <- C_02_02_irf_fn(sol_lst, B_01_05_irf_horizon_int)
  obs_lst <- C_06_04_sample_fn(par_lst, n_obs_int = 60L)
  pri_lst <- C_06_01_beta_fn(0.70, 0.15)
  ch_lst  <- C_07_03_chains_fn(par_lst, obs_lst, pri_lst$a_num,
                               pri_lst$b_num, n_draw_int = 120L,
                               warm_up_int = 30L)
  dens_df <- C_07_04_density_fn(ch_lst$kept_vec, pri_lst$a_num,
                                pri_lst$b_num)

  plot_lst <- list(
    persistence = D_02_01_roots_fn(lambda_vec = B_01_03_persistence_vec,
                                   n_periods_int = 21L),
    roots     = D_02_01_roots_fn(),
    branch    = D_02_02_branch_fn(par_lst),
    bubble    = D_04_01_bubble_fn(),
    approx    = D_04_02_approx_fn(),
    ratios    = D_04_03_ratios_fn(par_lst),
    cycles    = D_03_01_cycles_fn(sim_df),
    tech      = D_03_02_technology_fn(sim_df),
    growth    = D_03_03_growth_fn(sim_df),
    irf       = D_03_04_irf_fn(irf_df),
    irf_tech  = D_03_05_irf_tech_fn(irf_df),
    priors    = D_05_01_priors_fn(),
    trace     = D_05_02_trace_fn(ch_lst),
    posterior = D_05_03_posterior_fn(dens_df, pri_lst, ch_lst))

  series_vec <- toupper(unname(B_03_01_palette_vec[
    c("main", "second", "third", "compare")]))
  bad_chr <- character(0)

  for (nm_chr in names(plot_lst)) {
    b_lst <- ggplot_build(plot_lst[[nm_chr]])
    for (i in seq_along(b_lst$data)) {
      lyr <- plot_lst[[nm_chr]]$layers[[i]]
      if (!inherits(lyr$geom, "GeomLine")) next
      if (!is.null(lyr$aes_params$alpha)) next   # the ghost
      d_df <- b_lst$data[[i]]
      col_vec <- toupper(d_df$colour)
      if (!all(col_vec %in% series_vec)) {
        bad_chr <- c(bad_chr, paste0(nm_chr, ": off-palette ",
                                     paste(unique(col_vec[!col_vec %in%
                                                            series_vec]),
                                           collapse = " ")))
      }
      for (p_chr in unique(d_df$PANEL)) {
        one_df <- d_df[d_df$PANEL == p_chr, , drop = FALSE]
        if (length(unique(one_df$group)) !=
            length(unique(toupper(one_df$colour)))) {
          bad_chr <- c(bad_chr, paste0(nm_chr, ": two lines share a colour"))
        }
        if (!(series_vec[1L] %in% toupper(one_df$colour))) {
          bad_chr <- c(bad_chr, paste0(nm_chr, ": no blue first series"))
        }
      }
    }
    gt_obj <- ggplotGrob(plot_lst[[nm_chr]])
    if (any(grepl("guide-box", gt_obj$layout$name) &
            !vapply(gt_obj$grobs, inherits, TRUE, "zeroGrob"))) {
      bad_chr <- c(bad_chr, paste0(nm_chr, ": keeps a legend"))
    }
  }

  V_01_03_check_fn(
    "every figure's lines follow the decks' series order",
    length(bad_chr) == 0L,
    if (length(bad_chr) == 0L) {
      sprintf("%d figures, blue first, no shared colour, no legend",
              length(plot_lst))
    } else paste(unique(bad_chr), collapse = "; "))

  # the two comparator lines: the other root and the prior
  col_fn <- function(p_obj, lvl_chr) {
    sc <- p_obj$scales$get_scales("colour")
    toupper(sc$palette(0)[[lvl_chr]])
  }
  V_01_03_check_fn(
    "the other root and the prior are the light-blue comparator",
    identical(col_fn(plot_lst$branch, "Other root"),
              toupper(B_03_01_palette_vec[["compare"]])) &&
      identical(col_fn(plot_lst$posterior, "Prior"),
                toupper(B_03_01_palette_vec[["compare"]])),
    toupper(B_03_01_palette_vec[["compare"]]))

  # a PNG through the same export function, at the pair shape: 3:2, at the
  # pair panel's declared pixel size
  png_chr <- file.path(tempdir(), "rbc-verify.png")
  D_01_01_export_fn(plot_lst$irf, png_chr, "pair", "png")
  dim_vec <- V_09_02_pngdim_fn(png_chr)
  want_vec <- c(B_03_02_shape_lst$pair$width_px,
                B_03_02_shape_lst$pair$height_px)
  V_01_03_check_fn("the pair export through T_02_03c_export_fn is 3:2",
                   all(abs(dim_vec - want_vec) <= 1),
                   sprintf("%d x %d px", dim_vec[["width"]],
                           dim_vec[["height"]]))
})

###### V_09_13: The eta Preset's Wording #######################################
# Note: The story calls eta a wealth effect, which V_07b checks on the model.

local({
  story_chr <- B_02_02_example_lst$elastic$story
  V_01_03_check_fn(
    "the eta preset calls eta a wealth effect, not a supply slope",
    !grepl("flatter", story_chr, ignore.case = TRUE) &&
      grepl("wealth effect", story_chr, fixed = TRUE),
    B_02_02_example_lst$elastic$label)
})

###### V_09_14: The Page Frame #################################################
# Note: The UI rendered to HTML once: every plotOutput sits in a toolkit
#   figure card (T_07_07f_figcard_fn), every card carries a Save PNG button
#   in its header, no local CSS class remains, and no control label or card
#   header spells a Greek letter out.

local({
  html_chr <- tryCatch(htmltools::renderTags(E_02_04_ui)$html,
                       error = function(e) NA_character_)
  V_01_03_check_fn("the page renders to HTML without error",
                   !is.na(html_chr) && nchar(html_chr) > 1000L,
                   sprintf("%d characters", nchar(html_chr)))
  if (is.na(html_chr)) return(invisible(NULL))

  # every plot output's id, and the id inside each card
  plot_re  <- 'class="shiny-plot-output[^"]*" id="([^"]+)"'
  plot_vec <- regmatches(html_chr, gregexpr(plot_re, html_chr))[[1L]]
  plot_vec <- sub(paste0("^", plot_re, "$"), "\\1", plot_vec)
  card_vec <- strsplit(html_chr,
                       'class="card bslib-card[^"]*fig-card')[[1L]][-1L]
  in_card_vec <- vapply(card_vec, function(c_chr) {
    m_chr <- regmatches(c_chr, regexpr(plot_re, c_chr))
    if (length(m_chr) == 0L) NA_character_ else
      sub(paste0("^", plot_re, "$"), "\\1", m_chr)
  }, "", USE.NAMES = FALSE)
  V_01_03_check_fn(
    "every plotOutput sits in a T_07_07f figure card",
    setequal(plot_vec, names(B_03_05_figfile_lst)) &&
      length(plot_vec) == length(B_03_05_figfile_lst) &&
      setequal(in_card_vec, plot_vec) &&
      length(in_card_vec) == length(plot_vec),
    sprintf("%d plot(s), %d card(s)", length(plot_vec), length(in_card_vec)))

  btn_ok_lgl <- vapply(names(B_03_05_figfile_lst), function(id_chr) {
    grepl(sprintf('id="%s__png"', id_chr), html_chr, fixed = TRUE)
  }, TRUE)
  V_01_03_check_fn("every figure card carries a Save PNG button",
                   all(btn_ok_lgl),
                   if (all(btn_ok_lgl)) "all fourteen" else
                     paste(names(btn_ok_lgl)[!btn_ok_lgl], collapse = ", "))

  # the page's CSS is the toolkit's
  V_01_03_check_fn("the app carries no local CSS classes of its own",
                   !grepl("dg-(card|figbox|stat|story|prompt|dl)", html_chr),
                   "no dg-* classes on the page")

  greek_re <- paste0("\\b(alpha|beta|gamma|delta|epsilon|eta|theta|lambda|",
                     "mu|nu|xi|pi|rho|sigma|tau|phi|varphi|psi|omega)\\b")
  label_vec  <- vapply(B_03_07_controls_lst, function(x) x$label, "")
  header_vec <- unlist(B_03_09_header_lst)
  head_html_vec <- regmatches(html_chr, gregexpr(
    '<div class="card-header[^"]*">[^<]*', html_chr))[[1L]]
  head_html_vec <- sub("^<div[^>]*>", "", head_html_vec)
  stage_vec  <- names(B_02_01_stage_vec)
  text_vec   <- c(label_vec, header_vec, head_html_vec, stage_vec)
  # entities stripped first, so "&rho;" is not read as the word "rho"
  bare_vec   <- gsub("&[a-zA-Z]+;", " ", text_vec)
  hit_vec    <- text_vec[grepl(greek_re, bare_vec, ignore.case = TRUE)]
  V_01_03_check_fn(
    "no slider label or card header spells out a Greek letter",
    length(hit_vec) == 0L,
    if (length(hit_vec) == 0L) {
      sprintf("%d label(s) and header(s)", length(text_vec))
    } else {
      paste(hit_vec, collapse = " | ")
    })
})

################################################################################
## V_10: The App, Driven #######################################################
################################################################################
# Note: The server run over every stage through testServer, and every
#   download handler fired and its file read back.

message("")
message("V_10  The app, driven through every stage")

###### V_10_01: The Sliders at Their Defaults ##################################
# Note: The model and estimation defaults as one input list.

V_10_01_input_lst <- c(B_01_01_default_lst, B_01_04_estimation_lst)

###### V_10_02: One Pass Over Every Stage ######################################
# Note: Two checks off one pass, since stage 2.4b runs both chains on every
#   visit: do the panels build and does the notation key grow, and does any
#   equation reach the page as raw markup.

local({
  built_lgl <- logical(0)
  row_int   <- integer(0)
  bad_chr   <- character(0)

  shiny::testServer(F_01_01_server, {
    for (s_chr in unname(B_02_01_stage_vec)) {
      do.call(session$setInputs, V_10_01_input_lst)
      session$setInputs(stage = s_chr)

      eq_chr <- paste(as.character(output$equation_ui), collapse = "")
      no_chr <- paste(as.character(output$notation_ui), collapse = "")
      ex_chr <- paste(as.character(output$explain_ui), collapse = "")

      n_fn <- function(pat_chr, in_chr) {
        lengths(regmatches(in_chr, gregexpr(pat_chr, in_chr)))
      }

      # the prompt and preset title must build; the tiles must not throw
      pr_chr <- paste(as.character(output$prompt), collapse = "")
      pt_chr <- paste(as.character(output$preset_title), collapse = "")
      invisible(output$tiles)

      built_lgl <<- c(built_lgl,
                      nchar(eq_chr) > 0L && nchar(no_chr) > 0L &&
                        nchar(ex_chr) > 0L && grepl("prompt", pr_chr) &&
                        grepl("Worked Example", pt_chr))
      # rows, not characters: a "new" badge lengthens the HTML
      row_int <<- c(row_int, n_fn("<tr>", no_chr))

      for (one_chr in c(eq_chr, no_chr, ex_chr)) {
        math_vec <- regmatches(
          one_chr, gregexpr("\\\\\\((?:(?!\\\\\\)).)*", one_chr,
                            perl = TRUE))[[1L]]
        if (any(grepl("[<>&](?!(amp|lt|gt);)", math_vec, perl = TRUE))) {
          bad_chr <<- c(bad_chr, s_chr)
        }
      }
    }
  })

  V_01_03_check_fn(
    "every stage builds its equations, notation, in-words and prompt panels",
    # non-decreasing, since 2.1b introduces no symbol of its own; it must
    # still grow overall
    all(built_lgl) && all(diff(row_int) >= 0L) &&
      row_int[length(row_int)] > row_int[1L],
    sprintf("notation panel grows %s rows",
            paste(row_int, collapse = " -> ")))

  # a "<" or "&" inside an equation marked as HTML opens a tag
  V_01_03_check_fn(
    "no equation reaches the page as raw markup",
    length(bad_chr) == 0L,
    if (length(bad_chr) == 0L) "every stage, all three panels" else
      paste("raw <, > or & at:", paste(unique(bad_chr), collapse = ", ")))
})

###### V_10_03: Every Download Handler Writes Its PNG ##########################
# Note: Each file is copied out of testServer's temporary directory before it
#   is read. The stage is set to the figure's own, since the name carries it.
#   Every file exists, is the declared pixel size and has the name D_01_03
#   builds. The app exports PNG only.

V_10_03_stage_lst <- list(
  plot_persistence = "1.1a", plot_bubble = "1.5a",
  plot_roots       = "1.5b", plot_branch = "1.5b",
  plot_approx      = "2.1a", plot_ratios = "2.1b",
  plot_cycles      = "2.1c", plot_tech   = "2.1c",
  plot_growth      = "2.1c", plot_irf    = "2.1c",
  plot_irf_tech    = "2.1c",
  plot_priors      = "2.4a", plot_trace  = "2.4b",
  plot_posterior   = "2.4b")

V_10_03_keep_dir <- file.path(tempdir(), "rbc-verify-downloads")
dir.create(V_10_03_keep_dir, showWarnings = FALSE, recursive = TRUE)

local({
  got_lst <- list()
  shiny::testServer(F_01_01_server, {
    for (id_chr in names(V_10_03_stage_lst)) {
      do.call(session$setInputs, V_10_01_input_lst)
      session$setInputs(stage = V_10_03_stage_lst[[id_chr]])
      src_chr <- output[[paste0(id_chr, "__png")]]
      dst_chr <- file.path(V_10_03_keep_dir, basename(src_chr))
      file.copy(src_chr, dst_chr, overwrite = TRUE)
      got_lst[[id_chr]] <<- dst_chr
    }
  })

  id_vec   <- names(got_lst)
  path_vec <- unlist(got_lst)
  want_lst <- lapply(id_vec, function(k_chr) {
    s_lst <- B_03_02_shape_lst[[B_03_05_figfile_lst[[k_chr]]$shape_chr]]
    c(s_lst$width_px, s_lst$height_px)
  })
  name_vec <- vapply(id_vec, function(k_chr) {
    D_01_03_figfile_fn(k_chr, V_10_03_stage_lst[[k_chr]])
  }, "")

  exist_lgl <- file.exists(path_vec)
  size_num  <- ifelse(exist_lgl, file.info(path_vec)$size, 0)
  dim_lst   <- lapply(path_vec, V_09_02_pngdim_fn)
  ok_dim    <- mapply(function(d, w) all(abs(d - w) <= 1), dim_lst, want_lst)
  ratio_num <- vapply(dim_lst, function(d) d[["width"]] / d[["height"]], 0)
  ok_name   <- basename(path_vec) == name_vec

  V_01_03_check_fn(
    "every download handler writes a named PNG at the declared size",
    all(exist_lgl) && all(size_num > 0) && all(ok_dim) && all(ok_name),
    sprintf("%d file(s), e.g. %s", length(path_vec), basename(path_vec[[1]])))

  V_01_03_check_fn(
    "every exported PNG is at least 3:2",
    all(ratio_num >= B_03_03_min_ratio_num),
    sprintf("narrowest %.2f:1", min(ratio_num)))
})

################################################################################
## V_11: Verdict ###############################################################
################################################################################
# Note: A non-zero exit status if any check failed.

message("")
if (V_01_02_fail_int == 0L) {
  message("ALL CHECKS PASSED.")
} else {
  message(V_01_02_fail_int, " CHECK(S) FAILED.")
  quit(status = 1L)
}

#--------------------------------- Script End ---------------------------------#
