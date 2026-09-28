################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## The Real Business Cycle Model: Solver, Simulator, Filter and Sampler       ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Sourced automatically by app.R. Can be sourced alone from a lecture
##   .qmd so slide figures come from the same model:
##     source("R/model.R")
##
## Inputs:
##   None. Every function is a pure function of a parameter list.
##
## Outputs:
##   C_01_* the linear rational-expectations solver
##   C_02_* simulation, impulse responses and moments
##   C_03_* the scalar root paths of stages 1.1a and 1.5b
##   C_04_* diagnostics: root counting, residuals, warnings
##   C_05_* the layers before the simulation: the log-linear approximation,
##          the steady-state ratios, the forward solution and a bubble
##   C_06_* the estimation objects: beta prior shapes, the prior densities,
##          the state-space form, a SIMULATED observed sample and the
##          Kalman log-likelihood
##   C_07_* the random-walk Metropolis-Hastings sampler
##
## The model (Whelan, MA Advanced Macroeconomics, part 7, slides 22 to 28).
##   Seven log-deviation variables, in this order, which is the order of
##   every matrix row and of every column the simulator returns:
##
##     1 output y      2 consumption c   3 investment i   4 capital k
##     5 hours n       6 return r        7 technology a
##
##   Resource      y_t = (1 - s_i) c_t + s_i i_t
##   Capital       k_t = (1 - delta) k_{t-1} + delta i_t
##   Production    y_t = alpha k_{t-1} + (1 - alpha) n_t + a_t
##   Labour        n_t = y_t - eta c_t
##   Euler         c_t = E_t c_{t+1} - (1 / eta) E_t r_{t+1}
##   Return        r_t = mu (y_t - k_{t-1})
##   Technology    a_t = rho a_{t-1} + w_t
##
##   with s_i = alpha delta / (delta - 1 + 1/beta) the investment share and
##   mu = 1 - beta (1 - delta) the return coefficient, both from the steady
##   state (slides 25 to 27). Whelan writes the resource constraint on C*/Y*
##   and I*/Y*; s_i is I*/Y*. Written as a matrix system in structural form:
##
##     Chat Z_t = Ahat Z_{t-1} + Bhat E_t Z_{t+1} + D1 w_t
##
##   Whelan's RATS program rbc.prg solves this system by iterating on
##   C <- inv(I - B C) A from C = I (the Binder and Pesaran method of part
##   6). That map solves the matrix quadratic Bhat C^2 - Chat C + Ahat = 0,
##   which has more than one root, and the iteration converges to whichever
##   root its arithmetic path leads to; run in R it reaches the explosive
##   one. C_01_04 therefore reduces the model to its three-equation core,
##   takes the eigenvalues and selects the stable manifold, which is the
##   fixed point Whelan's figures show. The shock loading is the same object
##   as rbc.prg's, since (I - B C)^-1 Chat^-1 = [Chat - Bhat C]^-1. C_04_02
##   substitutes the solution back into the structural system and reports
##   the residual; tests/verify_against_rats.R runs both methods.
##
##   The estimation objects (C_06, C_07) follow part 10, slides 11 to 12, for
##   the state-space form and measurement error, and the ECON42240 lecture
##   deck 2.4 for the random-walk Metropolis sampler. The sample the sampler
##   is run on is SIMULATED from this model; nothing here estimates a real
##   economy.
##
## Parameter list (par) elements:
##   alpha, beta, delta, eta, rho, sigma (optional; default sqrt(1 - rho^2))

#-------------------------------- Script Begin --------------------------------#

################################################################################
## A: Table of Contents ########################################################
################################################################################
# Note: C holds the model; app.R holds sections B, D, E, F and G.
#
#   C: Model
#     C_01_01  Derived parameters
#     C_01_02  Structural matrices
#     C_01_03  Reduce to the minimal system
#     C_01_04  Solve the model
#     C_02_01  Simulate a path
#     C_02_02  Impulse responses
#     C_02_03  Moments of a simulated path
#     C_03_01  Scalar paths at any vector of roots
#     C_04_01  Root counting
#     C_04_02  Solution residuals
#     C_04_03  Diagnostics and warnings
#     C_05_01  The log-linearisation trick, and its error
#     C_05_02  Steady-state ratios
#     C_05_03  Forward, backward and bubbles
#     C_06_01  Beta shape parameters from a mean and a standard deviation
#     C_06_02  The prior densities of lecture 2.4
#     C_06_03  The solved model in state-space form
#     C_06_04  A SIMULATED observed sample from this model
#     C_06_05  The Kalman log-likelihood
#     C_07_01  The log posterior for one persistence parameter
#     C_07_02  One random-walk Metropolis-Hastings chain
#     C_07_03  Two chains, the warm-up and the running means
#     C_07_04  Prior and posterior densities on one grid

################################################################################
## C: Model ####################################################################
################################################################################
# Note: Every variable is a log deviation from steady state, in per cent.
#   Nothing here touches Shiny, so the solver can be sourced by the figure
#   exporter and by a .qmd without loading the app.

#### C_01: Solver ##############################################################
# Note: The structural matrices, the reduction to a three-equation core and
#   the eigenvalue solution that selects the stable manifold.

###### C_01_01: Derived Parameters #############################################
# Note: The investment share and the return coefficient, implied by the
#   steady state [W7 25] to [W7 27] and never set directly.

C_01_01_derived_fn <- function(par) {
  list(
    s_i = (par$alpha * par$delta) / (par$delta - 1 + (par$beta)^(-1)),
    mu  = 1 - par$beta * (1 - par$delta)
  )
}

###### C_01_02: Structural Matrices ############################################
# Note: Chat, Ahat, Bhat and D1 for Chat Z_t = Ahat Z_{t-1} + Bhat E_t Z_{t+1}
#   + D1 w_t, one row per equation in the header's order. The row and column
#   order is C_01_02_names_vec; a new variable goes at the end.

C_01_02_names_vec <- c("output", "consumption", "investment", "capital",
                       "hours", "return", "technology")

C_01_02_matrices_fn <- function(par) {

  der_lst <- C_01_01_derived_fn(par)
  n_var   <- length(C_01_02_names_vec)
  zero_fn <- function() matrix(0, nrow = n_var, ncol = n_var)

  c_hat <- zero_fn()
  a_hat <- zero_fn()
  b_hat <- zero_fn()
  d_one <- matrix(0, nrow = n_var, ncol = 1)

  # 1 Resource constraint: y = (1 - s_i) c + s_i i
  c_hat[1, 1] <- 1
  c_hat[1, 2] <- -(1 - der_lst$s_i)
  c_hat[1, 3] <- -der_lst$s_i

  # 2 Capital accumulation: k_t = (1 - delta) k_{t-1} + delta i_t
  c_hat[2, 4] <- 1
  c_hat[2, 3] <- -par$delta
  a_hat[2, 4] <- 1 - par$delta

  # 3 Production: y = alpha k_{t-1} + (1 - alpha) n + a
  c_hat[3, 1] <- 1
  c_hat[3, 5] <- -(1 - par$alpha)
  c_hat[3, 7] <- -1
  a_hat[3, 4] <- par$alpha

  # 4 Labour supply against the marginal product: n = y - eta c
  c_hat[4, 5] <- 1
  c_hat[4, 1] <- -1
  c_hat[4, 2] <- par$eta

  # 5 Euler equation: c_t = E_t c_{t+1} - (1 / eta) E_t r_{t+1}
  c_hat[5, 2] <- 1
  b_hat[5, 2] <- 1
  b_hat[5, 6] <- -(par$eta)^(-1)

  # 6 Return on capital: r = mu (y - k_{t-1})
  c_hat[6, 6] <- 1
  c_hat[6, 1] <- -der_lst$mu
  a_hat[6, 4] <- -der_lst$mu

  # 7 Technology: a_t = rho a_{t-1} + w_t
  c_hat[7, 7] <- 1
  a_hat[7, 7] <- par$rho
  d_one[7, 1] <- 1

  list(c_hat = c_hat, a_hat = a_hat, b_hat = b_hat, d_one = d_one,
       s_i = der_lst$s_i, mu = der_lst$mu)
}

###### C_01_03: Reduce to the Minimal System ###################################
# Note: Four of the seven equations are static. Given k_{t-1}, a_t and c_t
#   they pin output, hours, the return and investment:
#
#     y_t = k_{t-1} + a_t / alpha - ((1 - alpha) eta / alpha) c_t
#     n_t = k_{t-1} + a_t / alpha - (eta / alpha) c_t
#     r_t = mu (a_t / alpha - ((1 - alpha) eta / alpha) c_t)
#     i_t = (y_t - (1 - s_i) c_t) / s_i
#
#   so the model is two predetermined variables and one jump around a
#   three-equation core:
#
#     k_t = A_kk k_{t-1} + A_ka a_t + A_kc c_t
#     a_t = rho a_{t-1} + w_t
#     c_t = phi_c E_t c_{t+1} - psi a_t
#
#   The Euler equation carries no capital: r_{t+1} depends on y_{t+1} - k_t
#   and the k_t terms cancel, which is why the count is one jump and not two.

C_01_03_reduce_fn <- function(par) {

  der_lst <- C_01_01_derived_fn(par)
  s_i     <- der_lst$s_i
  mu      <- der_lst$mu

  # Static block: output as a function of (k_{t-1}, a_t, c_t)
  gy_k_num <- 1
  gy_a_num <- 1 / par$alpha
  gy_c_num <- -((1 - par$alpha) * par$eta / par$alpha)

  list(
    s_i      = s_i,
    mu       = mu,
    gy_k     = gy_k_num,
    gy_a     = gy_a_num,
    gy_c     = gy_c_num,
    a_kk     = (1 - par$delta) + (par$delta / s_i) * gy_k_num,
    a_ka     = (par$delta / s_i) * gy_a_num,
    a_kc     = (par$delta / s_i) * (gy_c_num - (1 - s_i)),
    phi_c    = 1 + mu * (1 - par$alpha) / par$alpha,
    psi      = mu * par$rho / (par$eta * par$alpha)
  )
}

###### C_01_04: Solve the Model ################################################
# Note: Writes the core as M0 E_t v_{t+1} = M1 v_t in v = (k, a, c) and takes
#   the eigenvalues of W = inv(M0) M1. Blanchard and Kahn (1980): one jump
#   variable, so exactly one eigenvalue must sit outside the unit circle. The
#   row of inv(V) belonging to it annihilates any non-explosive path, and
#   setting it to zero pins consumption; that selects the stable manifold.
#   The seven-variable transition C and loading H are rebuilt from the static
#   block. Returns converged_lgl rather than stopping, so the app can say on
#   screen that a slider setting has no stable solution.
#   branch_chr = "unstable" zeroes a stable row instead, putting the model on
#   the divergent manifold for stage 1.5b. That solution still satisfies the
#   structural equations; it fails transversality, not the model.

C_01_04_solve_fn <- function(par, branch_chr = c("stable", "unstable")) {

  branch_chr <- match.arg(branch_chr)
  red_lst <- C_01_03_reduce_fn(par)
  mat_lst <- C_01_02_matrices_fn(par)
  n_var   <- length(C_01_02_names_vec)

  fail_lst <- list(c_mat = NULL, h_mat = NULL, mat_lst = mat_lst,
                   red_lst = red_lst, eigen_vec = NA_real_,
                   outside_n = NA_integer_, converged_lgl = FALSE)

  m_zero <- matrix(c(1, 0, 0,
                     -red_lst$a_ka, 1, 0,
                     -red_lst$a_kc, 0, red_lst$phi_c),
                   nrow = 3, ncol = 3)
  m_one  <- matrix(c(red_lst$a_kk, 0, 0,
                     0, par$rho, red_lst$psi,
                     0, 0, 1),
                   nrow = 3, ncol = 3)

  w_mat <- tryCatch(solve(m_zero, m_one), error = function(e) NULL)
  if (is.null(w_mat)) return(fail_lst)

  eig_lst   <- eigen(w_mat)
  mod_vec   <- Mod(eig_lst$values)
  outside_n <- sum(mod_vec > 1 + 1e-09)

  # one jump variable, so one explosive root; the unstable branch is exempt
  if (outside_n != 1L && identical(branch_chr, "stable")) {
    fail_lst$eigen_vec <- mod_vec
    fail_lst$outside_n <- outside_n
    return(fail_lst)
  }

  vec_inv <- tryCatch(solve(eig_lst$vectors), error = function(e) NULL)
  if (is.null(vec_inv)) return(fail_lst)

  # "stable" zeroes the explosive coordinate; "unstable" zeroes the stable
  # row with most weight on consumption, since a row without consumption
  # cannot be solved for it
  outside_int <- which(mod_vec > 1 + 1e-09)
  cons_vec    <- abs(Re(vec_inv[, 3L]))

  pick_int <- if (identical(branch_chr, "stable")) {
    outside_int[1L]
  } else {
    stable_int <- setdiff(seq_along(mod_vec), outside_int)
    if (length(stable_int) == 0L) return(fail_lst)
    stable_int[which.max(cons_vec[stable_int])]
  }
  if (length(pick_int) == 0L || is.na(pick_int)) return(fail_lst)

  row_vec <- Re(vec_inv[pick_int, ])
  if (abs(row_vec[3L]) < 1e-12) return(fail_lst)
  c_on_k  <- -row_vec[1L] / row_vec[3L]
  c_on_a  <- -row_vec[2L] / row_vec[3L]

  # c_t is written on k_t, which depends on c_t; resolve
  den_num <- 1 - c_on_k * red_lst$a_kc
  if (abs(den_num) < 1e-12) return(fail_lst)
  p_k_num <- c_on_k * red_lst$a_kk / den_num
  p_a_num <- (c_on_k * red_lst$a_ka + c_on_a) / den_num

  # rebuild all seven variables from (k_{t-1}, a_t)
  state_fn <- function(k_lag_num, a_num) {
    c_num <- p_k_num * k_lag_num + p_a_num * a_num
    y_num <- red_lst$gy_k * k_lag_num + red_lst$gy_a * a_num +
      red_lst$gy_c * c_num
    n_num <- k_lag_num + a_num / par$alpha - (par$eta / par$alpha) * c_num
    r_num <- red_lst$mu * (a_num / par$alpha + red_lst$gy_c * c_num)
    i_num <- (y_num - (1 - red_lst$s_i) * c_num) / red_lst$s_i
    k_num <- (1 - par$delta) * k_lag_num + par$delta * i_num
    c(y_num, c_num, i_num, k_num, n_num, r_num, a_num)
  }

  c_mat <- matrix(0, nrow = n_var, ncol = n_var)
  c_mat[, 4L] <- state_fn(1, 0)         # one unit of k_{t-1}
  c_mat[, 7L] <- state_fn(0, par$rho)   # a_{t-1} = 1 gives a_t = rho
  h_mat <- matrix(state_fn(0, 1), ncol = 1L)

  dimnames(c_mat) <- list(C_01_02_names_vec, C_01_02_names_vec)
  rownames(h_mat) <- C_01_02_names_vec

  list(
    c_mat         = c_mat,
    h_mat         = h_mat,
    mat_lst       = mat_lst,
    red_lst       = red_lst,
    eigen_vec     = sort(mod_vec, decreasing = TRUE),
    outside_n     = outside_n,
    policy_k      = p_k_num,
    policy_a      = p_a_num,
    branch_chr    = branch_chr,
    converged_lgl = TRUE
  )
}

#### C_02: Simulation and Responses ############################################
# Note: Both take a solved model and return a long-format data frame, so the
#   drawing functions never see a matrix.

###### C_02_01: Simulate a Path ################################################
# Note: Period 1 is the steady state and the first innovation lands in period
#   2, as in rbc.prg. The innovation standard deviation defaults to
#   sqrt(1 - rho^2), rbc.prg's stdeps, which holds the variance of technology
#   at one whatever rho is.

C_02_01_simulate_fn <- function(sol_lst, par, n_periods_int = 200L,
                                seed_int = 42L) {

  if (!isTRUE(sol_lst$converged_lgl)) return(NULL)

  set.seed(seed_int)
  n_var    <- length(C_01_02_names_vec)
  sigma_num <- if (is.null(par$sigma)) sqrt(1 - par$rho^2) else par$sigma
  shock_vec <- c(0, stats::rnorm(n_periods_int - 1L, 0, sigma_num))

  z_mat <- matrix(0, nrow = n_periods_int, ncol = n_var)
  z_vec <- matrix(0, nrow = n_var, ncol = 1)

  for (t in seq.int(2L, n_periods_int)) {
    z_vec      <- sol_lst$c_mat %*% z_vec + sol_lst$h_mat * shock_vec[t]
    z_mat[t, ] <- z_vec[, 1]
  }

  sim_df <- data.frame(mod_period_tm = seq_len(n_periods_int))
  for (j in seq_len(n_var)) {
    sim_df[[paste0("mod_", C_01_02_names_vec[j], "_val")]] <- z_mat[, j]
  }
  sim_df$mod_innovation_val <- shock_vec
  sim_df$mod_growth_val     <- c(NA_real_, diff(z_mat[, 1]))
  sim_df
}

###### C_02_02: Impulse Responses ##############################################
# Note: One unit innovation at horizon 0, nothing after: H, C H, C^2 H, ...
#   Whelan's [W7 36] and [W7 37] are columns of this.

C_02_02_irf_fn <- function(sol_lst, n_horizon_int = 25L) {

  if (!isTRUE(sol_lst$converged_lgl)) return(NULL)

  n_var <- length(C_01_02_names_vec)
  z_mat <- matrix(0, nrow = n_horizon_int, ncol = n_var)
  z_vec <- sol_lst$h_mat

  for (k in seq_len(n_horizon_int)) {
    z_mat[k, ] <- z_vec[, 1]
    z_vec      <- sol_lst$c_mat %*% z_vec
  }

  irf_df <- data.frame(mod_horizon_tm = seq_len(n_horizon_int) - 1L)
  for (j in seq_len(n_var)) {
    irf_df[[paste0("mod_", C_01_02_names_vec[j], "_val")]] <- z_mat[, j]
  }
  irf_df
}

###### C_02_03: Moments of a Simulated Path ####################################
# Note: Relative standard deviations, the growth autocorrelation and the
#   correlation of output with technology, after a burn-in.

C_02_03_moments_fn <- function(sim_df, burn_in_int = 20L) {

  if (is.null(sim_df)) return(NULL)

  keep_df <- sim_df[sim_df$mod_period_tm > burn_in_int, , drop = FALSE]
  sd_vec  <- vapply(C_01_02_names_vec, function(one_chr) {
    stats::sd(keep_df[[paste0("mod_", one_chr, "_val")]])
  }, numeric(1))

  growth_vec <- keep_df$mod_growth_val
  growth_vec <- growth_vec[!is.na(growth_vec)]
  n_growth   <- length(growth_vec)

  list(
    sd_vec        = sd_vec,
    sd_rel_vec    = sd_vec / sd_vec[["output"]],
    growth_acf    = if (n_growth > 2L) {
      stats::cor(growth_vec[-n_growth], growth_vec[-1L])
    } else NA_real_,
    tech_corr     = stats::cor(keep_df$mod_output_val,
                               keep_df$mod_technology_val)
  )
}

#### C_03: The Root Illustration ###############################################
# Note: A single difference equation, not the seven-variable model: the
#   smallest object that shows what persistence does to a shock (1.1a) and
#   what a root outside the unit circle does (1.5b).

###### C_03_01: Scalar Paths ###################################################
# Note: y_t = lambda y_{t-1} + w_t with one unit innovation and nothing after,
#   one block per root in long format. lambda_vec is named and its order is
#   the drawing order. n_periods_int is the number of rows, first_period_int
#   the index of the first row and shock_period_int the index the innovation
#   lands on. lambda_stable_num and lambda_explosive_num rebuild the stage
#   1.5b pair by name.

C_03_01_roots_fn <- function(lambda_vec = c("Inside the unit circle"  = 0.85,
                                            "Outside the unit circle" = 1.05),
                             n_periods_int        = 40L,
                             first_period_int     = 1L,
                             shock_period_int     = 2L,
                             lambda_stable_num    = NULL,
                             lambda_explosive_num = NULL) {

  # a stable-and-explosive pair given by name becomes the vector
  if (!is.null(lambda_stable_num) || !is.null(lambda_explosive_num)) {
    lambda_vec <- c(
      "Inside the unit circle"  =
        if (is.null(lambda_stable_num)) 0.85 else lambda_stable_num,
      "Outside the unit circle" =
        if (is.null(lambda_explosive_num)) 1.05 else lambda_explosive_num)
  }
  if (is.null(names(lambda_vec)) || any(!nzchar(names(lambda_vec)))) {
    names(lambda_vec) <- sprintf("lambda == %s", format(lambda_vec))
  }

  # a positional (stable, explosive, periods) call lands a root here
  if (length(n_periods_int) != 1L || is.na(n_periods_int) ||
      n_periods_int < 1 || n_periods_int != round(n_periods_int)) {
    stop("n_periods_int must be a single whole number of periods. The roots ",
         "are a named vector in the first argument; a positional (stable, ",
         "explosive, periods) call lands here. Name them: lambda_stable_num ",
         "= , lambda_explosive_num = .", call. = FALSE)
  }

  period_vec <- seq.int(first_period_int, length.out = n_periods_int)

  path_fn <- function(lambda_num) {
    y_vec <- numeric(n_periods_int)
    for (i in seq_len(n_periods_int)) {
      y_prev   <- if (i == 1L) 0 else y_vec[i - 1L]
      y_vec[i] <- lambda_num * y_prev +
        as.numeric(period_vec[i] == shock_period_int)
    }
    y_vec
  }

  do.call(rbind, lapply(names(lambda_vec), function(one_chr) {
    data.frame(mod_period_tm = period_vec,
               mod_path_val  = path_fn(lambda_vec[[one_chr]]),
               mod_root_val  = lambda_vec[[one_chr]],
               mod_case_cat  = one_chr)
  }))
}

#### C_04: Diagnostics #########################################################
# Note: Whether a parameter setting has a unique stable solution, and whether
#   the solution satisfies the model.

###### C_04_01: Root Counting ##################################################
# Note: Blanchard and Kahn (1980), counted on the three roots of the reduced
#   system: two predetermined variables (capital and technology) and one
#   jump (consumption). A root of modulus exactly one is reported separately.

C_04_01_roots_fn <- function(sol_lst, tol_num = 1e-08) {

  if (!isTRUE(sol_lst$converged_lgl)) return(NULL)

  # the three roots of the reduced system; the seven-variable transition
  # carries five zeros for the variables pinned by construction
  mod_vec <- sort(sol_lst$eigen_vec, decreasing = TRUE)

  list(
    modulus_vec     = mod_vec,
    inside_n        = sum(mod_vec < 1 - tol_num),
    on_circle_n     = sum(abs(mod_vec - 1) <= tol_num),
    outside_n       = sum(mod_vec > 1 + tol_num),
    predetermined_n = 2L,
    jump_n          = 1L,
    full_modulus_vec = sort(Mod(eigen(sol_lst$c_mat,
                                      only.values = TRUE)$values),
                            decreasing = TRUE)
  )
}

###### C_04_02: Solution Residuals #############################################
# Note: The solution substituted back into the structural system:
#     state  Chat C  - (Ahat + Bhat C C)
#     shock  Chat H  - (Bhat C H + D1)
#   Both are zero for an exact solution; anything above 1e-06 is a bug.

C_04_02_residual_fn <- function(sol_lst) {

  if (!isTRUE(sol_lst$converged_lgl)) return(NULL)

  mat_lst <- sol_lst$mat_lst
  c_mat   <- sol_lst$c_mat
  h_mat   <- sol_lst$h_mat

  state_mat <- mat_lst$c_hat %*% c_mat -
    (mat_lst$a_hat + mat_lst$b_hat %*% c_mat %*% c_mat)
  shock_mat <- mat_lst$c_hat %*% h_mat -
    (mat_lst$b_hat %*% c_mat %*% h_mat + mat_lst$d_one)

  list(state_num = max(abs(state_mat)), shock_num = max(abs(shock_mat)))
}

###### C_04_03: Diagnostics and Warnings #######################################
# Note: Everything the app puts on screen about the solution, and the
#   problems vector naming the slider to move.

C_04_03_diagnostics_fn <- function(sol_lst, par) {

  problems_vec <- character(0)

  if (!isTRUE(sol_lst$converged_lgl)) {
    out_n <- sol_lst$outside_n
    problems_vec <- c(problems_vec, if (!is.na(out_n) && out_n != 1L) {
      paste0(
        "No unique stable solution here. The model has one jump variable, ",
        "consumption, so it needs exactly one root outside the unit circle ",
        "and this calibration gives ", out_n, ". Bring &rho; back below one, ",
        "or move &alpha; back towards a third."
      )
    } else {
      paste0(
        "The model could not be solved at these parameters. Bring &rho; back ",
        "below one, or move alpha back towards a third."
      )
    })
    return(list(roots_lst = NULL, residual_lst = NULL,
                s_i = NA_real_, mu = NA_real_, problems = problems_vec))
  }

  roots_lst <- C_04_01_roots_fn(sol_lst)
  resid_lst <- C_04_02_residual_fn(sol_lst)

  if (roots_lst$outside_n != roots_lst$jump_n) {
    problems_vec <- c(problems_vec, paste0(
      "The reduced system carries ", roots_lst$outside_n,
      " roots outside the unit circle where the model has ",
      roots_lst$jump_n, " jump variable. The solution on screen is not ",
      "the unique stable one."
    ))
  }
  if (roots_lst$on_circle_n > 0L) {
    problems_vec <- c(problems_vec, paste0(
      "A root sits exactly on the unit circle, which the counting rule ",
      "does not cover. Move &rho; off its current value to see which side ",
      "the model falls."
    ))
  }
  if (max(resid_lst$state_num, resid_lst$shock_num) > 1e-06) {
    problems_vec <- c(problems_vec, paste0(
      "The solution does not satisfy the structural equations it came ",
      # the em dash is an escape; see CONVENTIONS.md 6
      "from. This is a bug, not a parameter choice \u2014 report it."
    ))
  }

  list(
    roots_lst    = roots_lst,
    residual_lst = resid_lst,
    s_i          = sol_lst$mat_lst$s_i,
    mu           = sol_lst$mat_lst$mu,
    problems     = problems_vec
  )
}

#### C_05: The Layers Before the Simulation ####################################
# Note: Whelan builds part 7 in five labelled parts and simulates only in the
#   last; these functions carry the earlier layers that earn a figure.

###### C_05_01: The Log-Linearisation Trick ####################################
# Note: [W7 19]: exp(x) is about 1 + x around the steady state. Returns the
#   exact curve, the approximation and the error over a range of log
#   deviations, as a share of the true level.

C_05_01_approx_fn <- function(max_dev_num = 0.6, n_point_int = 241L) {
  x_vec <- seq(-max_dev_num, max_dev_num, length.out = n_point_int)
  data.frame(
    mod_deviation_val  = x_vec,
    mod_exact_val      = exp(x_vec),
    mod_linear_val     = 1 + x_vec,
    mod_error_val      = exp(x_vec) - (1 + x_vec),
    mod_error_pct_val  = 100 * (exp(x_vec) - (1 + x_vec)) / exp(x_vec)
  )
}

###### C_05_02: Steady-State Ratios ############################################
# Note: [W7 25] to [W7 27]. From R* = 1 / beta and R* = alpha Y*/K* + 1 -
#   delta,
#     Y*/K* = (1 / beta + delta - 1) / alpha
#     I*/K* = delta
#     I*/Y* = alpha delta / (1 / beta + delta - 1)
#   The last is the solver's s_i. Returns a sweep over one parameter.

C_05_02_ratios_fn <- function(par, sweep_chr = c("delta", "beta"),
                              n_point_int = 121L) {

  sweep_chr <- match.arg(sweep_chr)
  grid_vec  <- if (identical(sweep_chr, "delta")) {
    seq(0.005, 0.040, length.out = n_point_int)
  } else {
    seq(0.950, 0.999, length.out = n_point_int)
  }

  out_df <- data.frame(mod_sweep_val = grid_vec)
  out_df$mod_y_over_k_val <- NA_real_
  out_df$mod_i_over_y_val <- NA_real_

  for (i in seq_along(grid_vec)) {
    par_i <- par
    par_i[[sweep_chr]] <- grid_vec[i]
    out_df$mod_y_over_k_val[i] <- ((par_i$beta)^(-1) + par_i$delta - 1) /
      par_i$alpha
    out_df$mod_i_over_y_val[i] <- C_01_01_derived_fn(par_i)$s_i
  }
  out_df$mod_i_over_k_val <- if (identical(sweep_chr, "delta")) {
    grid_vec
  } else {
    rep(par$delta, n_point_int)
  }
  out_df$mod_sweep_cat <- sweep_chr
  out_df
}

###### C_05_03: Forward, Backward and Bubbles ##################################
# Note: [W6 9] to [W6 14]. y_t = a E_t y_{t+1} + x_t has the forward solution
#   and, on top of it, any b_t with b_t = a E_t b_{t+1}: a rational bubble
#   growing at 1 / a. Returns both paths from a constant x.

C_05_03_bubble_fn <- function(a_num = 0.95, x_num = 1,
                              bubble_seed_num = 0.10,
                              n_periods_int = 40L) {

  fundamental_num <- x_num / (1 - a_num)
  bubble_vec      <- bubble_seed_num * (1 / a_num)^(seq_len(n_periods_int) - 1)

  rbind(
    data.frame(mod_period_tm = seq_len(n_periods_int),
               mod_value_val = rep(fundamental_num, n_periods_int),
               mod_path_cat  = "Fundamental"),
    data.frame(mod_period_tm = seq_len(n_periods_int),
               mod_value_val = fundamental_num + bubble_vec,
               mod_path_cat  = "With a bubble")
  )
}

#### C_06: The Estimation Objects ##############################################
# Note: The solved model rewritten as a state-space system, as Whelan does on
#   [W10 11] and [W10 12]: k_t = a_kk k_{t-1} + a_kz z_t, c_t = a_ck k_{t-1}
#   + a_cz z_t, z_t = rho z_{t-1} + epsilon_t, with capital and consumption
#   observed with error. Every coefficient comes out of C_01_04_solve_fn.
#   The sample the sampler is run on is SIMULATED from this model at a known
#   persistence, and C_06_04 returns simulated_lgl = TRUE and the truth so
#   no drawing layer can present it as data. There is no estimate of any
#   real economy here.

###### C_06_01: Beta Shapes From a Mean and a Standard Deviation ###############
# Note: Method of moments, with k = a + b:
#     a = mean k,   b = (1 - mean) k,   k = mean (1 - mean) / sd^2 - 1
#   A beta with this mean exists only for sd^2 < mean (1 - mean), and a shape
#   below one spikes at a boundary, so both are clamped and the prior
#   actually used is returned.

C_06_01_beta_fn <- function(mean_num = 0.70, sd_num = 0.15,
                            min_shape_num = 1.05) {

  mean_num <- min(max(mean_num, 1e-03), 1 - 1e-03)
  edge_num <- min(mean_num, 1 - mean_num)
  conc_num <- mean_num * (1 - mean_num) / (sd_num^2) - 1
  conc_num <- max(conc_num, min_shape_num / edge_num)

  a_num <- mean_num * conc_num
  b_num <- (1 - mean_num) * conc_num

  list(
    a_num    = a_num,
    b_num    = b_num,
    mean_num = a_num / (a_num + b_num),
    sd_num   = sqrt(a_num * b_num /
                      (((a_num + b_num)^2) * (a_num + b_num + 1))),
    conc_num = conc_num
  )
}

###### C_06_02: The Prior Densities ############################################
# Note: Whelan's three prior families, [W10 17] to [W10 19]: normal for an
#   unrestricted parameter, gamma for a positive one, beta for one on the
#   unit interval. Four curves, each drawn wider than its support;
#   mod_inside_lgl marks the support. The parameters are illustrative.

C_06_02_priors_fn <- function(n_point_int = 601L) {

  one_fn <- function(name_chr, family_chr, from_num, to_num,
                     dens_fn, lo_num, hi_num, mode_num) {
    x_vec <- seq(from_num, to_num, length.out = n_point_int)
    data.frame(
      mod_value_val   = x_vec,
      mod_density_val = dens_fn(x_vec),
      mod_inside_lgl  = x_vec >= lo_num & x_vec <= hi_num,
      mod_mode_val    = mode_num,
      mod_peak_val    = dens_fn(mode_num),
      mod_curve_cat   = name_chr,
      mod_family_cat  = family_chr,
      stringsAsFactors = FALSE)
  }

  out_df <- rbind(
    one_fn("Normal(0, 1)", "Normal", -4, 4,
           function(x) stats::dnorm(x, 0, 1), -Inf, Inf, 0),
    one_fn("Gamma(2, 4)", "Gamma", -0.6, 2.0,
           function(x) stats::dgamma(x, shape = 2, rate = 4), 0, Inf, 0.25),
    one_fn("Beta(2, 2)", "Beta", -0.3, 1.3,
           function(x) stats::dbeta(x, 2, 2), 0, 1, 0.5),
    one_fn("Beta(5, 2)", "Beta", -0.3, 1.3,
           function(x) stats::dbeta(x, 5, 2), 0, 1, 0.8))

  out_df$mod_family_cat <- factor(out_df$mod_family_cat,
                                  levels = c("Normal", "Gamma", "Beta"))
  out_df
}

###### C_06_03: The Solved Model in State-Space Form ###########################
# Note: [W10 12]. With the state s_t = (k_{t-1}, z_t),
#     s_t = F s_{t-1} + (0, epsilon_t)'      F = [[a_kk, a_kz], [0, rho]]
#     x_t = H s_t + (u^k_{t-1}, u^c_t)'      H = [[1, 0], [a_ck, a_cz]]
#   a_kk is c_mat["capital", "capital"], a_kz is h_mat["capital", ], and
#   a_ck and a_cz are the consumption policy. Returned as a numeric vector,
#   since C_06_05 reads it thousands of times.

C_06_03_space_fn <- function(par) {

  sol_lst <- C_01_04_solve_fn(par)
  if (!isTRUE(sol_lst$converged_lgl)) return(NULL)

  c(a_kk = unname(sol_lst$c_mat["capital", "capital"]),
    a_kz = unname(sol_lst$h_mat["capital", 1L]),
    rho  = par$rho,
    a_ck = sol_lst$policy_k,
    a_cz = sol_lst$policy_a)
}

###### C_06_04: A SIMULATED Observed Sample From This Model ####################
# Note: The sample is simulated from this model at the given parameters and
#   is not data. The innovation standard deviation is sigma_z sqrt(1 -
#   rho^2), C_02_01's convention. Two measurement errors answer the
#   stochastic singularity of [W10 9]: one shock and two observables need
#   two more random terms. Their standard deviations are a fraction of each
#   series' own and are treated as known; only rho is estimated.

C_06_04_sample_fn <- function(par, n_obs_int = 120L, me_frac_num = 0.50,
                              sigma_z_num = 1, seed_int = 20240L) {

  space_vec <- C_06_03_space_fn(par)
  if (is.null(space_vec)) return(NULL)

  set.seed(seed_int)
  sigma_num <- sigma_z_num * sqrt(max(1 - par$rho^2, 0))
  shock_vec <- stats::rnorm(n_obs_int, 0, sigma_num)

  k_vec <- numeric(n_obs_int)
  z_vec <- numeric(n_obs_int)
  c_vec <- numeric(n_obs_int)
  s_one <- 0
  s_two <- 0

  for (t in seq_len(n_obs_int)) {
    s_one <- space_vec[["a_kk"]] * s_one + space_vec[["a_kz"]] * s_two
    s_two <- space_vec[["rho"]] * s_two + shock_vec[t]
    k_vec[t] <- s_one
    z_vec[t] <- s_two
    c_vec[t] <- space_vec[["a_ck"]] * s_one + space_vec[["a_cz"]] * s_two
  }

  sd_k_num <- me_frac_num * stats::sd(k_vec)
  sd_c_num <- me_frac_num * stats::sd(c_vec)

  list(
    k_obs_vec     = k_vec + stats::rnorm(n_obs_int, 0, sd_k_num),
    c_obs_vec     = c_vec + stats::rnorm(n_obs_int, 0, sd_c_num),
    var_k_num     = sd_k_num^2,
    var_c_num     = sd_c_num^2,
    n_obs_int     = n_obs_int,
    sigma_z_num   = sigma_z_num,
    truth_num     = par$rho,
    simulated_lgl = TRUE
  )
}

###### C_06_05: The Kalman Log-Likelihood ######################################
# Note: The prediction error decomposition (Hamilton 1994, ch. 13; [W10 13]):
#     log L = -(1/2) sum [ log|Omega_t| + eps_t' Omega_t^-1 eps_t ]
#             - (mT/2) log 2pi
#   Written in scalars, since every matrix is 2 by 2. The Riccati recursion
#   carries no data, so it is frozen once it converges; tests V_08 checks
#   this against a matrix filter. The initial covariance is the
#   unconditional one, from P = F P F' + Q.

C_06_05_loglik_fn <- function(space_vec, obs_lst, sigma2_num,
                              tol_num = 1e-12) {

  if (is.null(space_vec) || is.null(obs_lst)) return(-Inf)

  f_one <- space_vec[[1L]]   # a_kk
  f_two <- space_vec[[2L]]   # a_kz
  f_rho <- space_vec[[3L]]   # rho
  h_one <- space_vec[[4L]]   # a_ck
  h_two <- space_vec[[5L]]   # a_cz

  if (abs(f_rho) >= 1 || abs(f_one) >= 1) return(-Inf)

  p_22 <- sigma2_num / (1 - f_rho * f_rho)
  p_12 <- f_two * f_rho * p_22 / (1 - f_one * f_rho)
  p_11 <- (2 * f_one * f_two * p_12 + f_two * f_two * p_22) /
    (1 - f_one * f_one)
  if (!is.finite(p_11) || !is.finite(p_12) || p_11 < 0) return(-Inf)

  x_one_vec <- obs_lst$k_obs_vec
  x_two_vec <- obs_lst$c_obs_vec
  r_one     <- obs_lst$var_k_num
  r_two     <- obs_lst$var_c_num
  n_obs     <- length(x_one_vec)

  s_one <- 0
  s_two <- 0
  ll_num <- 0
  frozen_lgl <- FALSE
  ldet_num <- 0
  k_11 <- 0; k_12 <- 0; k_21 <- 0; k_22 <- 0
  o_11 <- 0; o_12 <- 0; o_22 <- 0

  for (t in seq_len(n_obs)) {

    if (!frozen_lgl) {
      # Predicted state covariance: F P F' + Q, all three distinct entries.
      a_num <- f_one * p_11 + f_two * p_12
      b_num <- f_one * p_12 + f_two * p_22
      q_11  <- f_one * a_num + f_two * b_num
      q_12  <- f_rho * b_num
      q_22  <- f_rho * f_rho * p_22 + sigma2_num
      # Innovation covariance Omega = H Q H' + R, and its inverse.
      g_12 <- h_one * q_11 + h_two * q_12
      g_22 <- h_one * q_12 + h_two * q_22
      w_11 <- q_11 + r_one
      w_12 <- g_12
      w_22 <- h_one * g_12 + h_two * g_22 + r_two
      det_num <- w_11 * w_22 - w_12 * w_12
      if (!is.finite(det_num) || det_num <= 0) return(-Inf)
      ldet_num <- log(det_num)
      o_11 <- w_22 / det_num
      o_12 <- -w_12 / det_num
      o_22 <- w_11 / det_num
      # Kalman gain K = Q H' Omega^-1.
      k_11 <- (q_11 * w_22 - g_12 * w_12) / det_num
      k_12 <- (-q_11 * w_12 + g_12 * w_11) / det_num
      k_21 <- (q_12 * w_22 - g_22 * w_12) / det_num
      k_22 <- (-q_12 * w_12 + g_22 * w_11) / det_num
      # Updated covariance P = Q - K H Q.
      n_11 <- q_11 - (k_11 * q_11 + k_12 * g_12)
      n_12 <- q_12 - (k_11 * q_12 + k_12 * g_22)
      n_22 <- q_22 - (k_21 * q_12 + k_22 * g_22)
      if (abs(n_11 - p_11) < tol_num && abs(n_12 - p_12) < tol_num &&
          abs(n_22 - p_22) < tol_num) frozen_lgl <- TRUE
      p_11 <- n_11; p_12 <- n_12; p_22 <- n_22
    }

    t_one <- f_one * s_one + f_two * s_two
    t_two <- f_rho * s_two
    e_one <- x_one_vec[t] - t_one
    e_two <- x_two_vec[t] - (h_one * t_one + h_two * t_two)
    ll_num <- ll_num - 0.5 * (ldet_num + o_11 * e_one * e_one +
                                2 * o_12 * e_one * e_two +
                                o_22 * e_two * e_two)
    s_one <- t_one + k_11 * e_one + k_12 * e_two
    s_two <- t_two + k_21 * e_one + k_22 * e_two
  }

  ll_num - n_obs * log(2 * pi)
}

#### C_07: The Random-Walk Metropolis-Hastings Sampler #########################
# Note: Deck 2.4's algorithm, as written there:
#     1  propose theta* = theta^(s) + nu,  nu ~ N(0, c^2)
#     2  alpha = min{1, [L(theta*) p(theta*)] / [L(theta^(s)) p(theta^(s))]}
#     3  draw w uniform on the unit interval
#     4  set theta^(s+1) = theta* if w < alpha, otherwise theta^(s)
#   One parameter is estimated, rho, with the rest held at the values the
#   sample was simulated at. The acceptance ratio is the symmetric-proposal
#   form. There is no adaptation: the student tunes the step on a slider.
#   Whelan's part 10 stops at the posterior; the sampler is the deck's.

###### C_07_01: The Log Posterior for One Persistence Parameter ################
# Note: Log prior plus log likelihood. A proposal off the unit interval
#   returns -Inf before the model is solved.

C_07_01_logpost_fn <- function(rho_num, par, obs_lst, a_num, b_num) {

  if (!is.finite(rho_num) || rho_num <= 0 || rho_num >= 1) return(-Inf)

  prior_num <- stats::dbeta(rho_num, a_num, b_num, log = TRUE)
  if (!is.finite(prior_num)) return(-Inf)

  space_vec <- C_06_03_space_fn(modifyList(par, list(rho = rho_num)))
  if (is.null(space_vec)) return(-Inf)

  sigma2_num <- (obs_lst$sigma_z_num^2) * (1 - rho_num^2)
  C_06_05_loglik_fn(space_vec, obs_lst, sigma2_num) + prior_num
}

###### C_07_02: One Chain ######################################################
# Note: Returns every draw, including the warm-up. The uniform and normal
#   draws are taken in one call each.

C_07_02_chain_fn <- function(par, obs_lst, a_num, b_num,
                             step_num = 0.05, start_num = 0.30,
                             n_draw_int = 1000L, seed_int = 101L) {

  set.seed(seed_int)
  draw_vec  <- numeric(n_draw_int)
  step_vec  <- stats::rnorm(n_draw_int, 0, step_num)
  unif_vec  <- log(stats::runif(n_draw_int))
  accept_n  <- 0L

  now_num <- start_num
  now_lp  <- C_07_01_logpost_fn(now_num, par, obs_lst, a_num, b_num)

  for (s in seq_len(n_draw_int)) {
    cand_num <- now_num + step_vec[s]
    cand_lp  <- C_07_01_logpost_fn(cand_num, par, obs_lst, a_num, b_num)
    if (unif_vec[s] < cand_lp - now_lp) {
      now_num  <- cand_num
      now_lp   <- cand_lp
      accept_n <- accept_n + 1L
    }
    draw_vec[s] <- now_num
  }

  list(draw_vec = draw_vec, accept_num = accept_n / n_draw_int,
       start_num = start_num)
}

###### C_07_03: Two Chains, the Warm-Up and the Running Means ##################
# Note: Two chains started far apart, the cheapest convergence check. Returns
#   one block per chain with the running mean alongside the draw, plus the
#   kept draws, acceptance rates and posterior moments.

C_07_03_chains_fn <- function(par, obs_lst, a_num, b_num,
                              step_num    = 0.05,
                              start_vec   = c(0.30, 0.99),
                              n_draw_int  = 1000L,
                              warm_up_int = 250L,
                              seed_vec    = c(101L, 202L)) {

  if (is.null(obs_lst)) return(NULL)

  chain_lst <- lapply(seq_along(start_vec), function(i) {
    C_07_02_chain_fn(par, obs_lst, a_num, b_num, step_num,
                     start_vec[i], n_draw_int, seed_vec[i])
  })

  # the start distinguishes the two chains, so it is in the name
  name_vec <- sprintf("Chain %d (from %.2f)", seq_along(start_vec),
                      start_vec)

  draw_df <- do.call(rbind, lapply(seq_along(chain_lst), function(i) {
    one_vec <- chain_lst[[i]]$draw_vec
    data.frame(
      mod_draw_tm      = seq_len(n_draw_int),
      mod_rho_val      = one_vec,
      mod_running_val  = cumsum(one_vec) / seq_len(n_draw_int),
      mod_chain_cat    = name_vec[i],
      stringsAsFactors = FALSE)
  }))

  kept_vec <- unlist(lapply(chain_lst, function(one_lst) {
    one_lst$draw_vec[seq.int(warm_up_int + 1L, n_draw_int)]
  }), use.names = FALSE)

  list(
    draw_df       = draw_df,
    kept_vec      = kept_vec,
    accept_vec    = vapply(chain_lst, function(x) x$accept_num, numeric(1)),
    warm_up_int   = warm_up_int,
    n_draw_int    = n_draw_int,
    step_num      = step_num,
    post_mean_num = mean(kept_vec),
    post_sd_num   = stats::sd(kept_vec),
    truth_num     = obs_lst$truth_num,
    simulated_lgl = isTRUE(obs_lst$simulated_lgl)
  )
}

###### C_07_04: Prior and Posterior on One Grid ################################
# Note: The prior density and a kernel density of the kept draws. The kernel
#   is evaluated on the unit interval and renormalised there, since it leaks
#   past one when the posterior sits near the boundary.

C_07_04_density_fn <- function(kept_vec, a_num, b_num, n_point_int = 512L) {

  if (is.null(kept_vec) || length(kept_vec) < 10L) return(NULL)

  grid_vec <- seq(0, 1, length.out = n_point_int)
  kde_lst  <- stats::density(kept_vec, from = 0, to = 1, n = n_point_int)
  post_vec <- pmax(kde_lst$y, 0)
  area_num <- sum((post_vec[-1L] + post_vec[-n_point_int]) / 2) *
    (grid_vec[2L] - grid_vec[1L])
  if (area_num > 0) post_vec <- post_vec / area_num

  rbind(
    data.frame(mod_value_val   = grid_vec,
               mod_density_val = stats::dbeta(grid_vec, a_num, b_num),
               mod_curve_cat   = "Prior",
               stringsAsFactors = FALSE),
    data.frame(mod_value_val   = kde_lst$x,
               mod_density_val = post_vec,
               mod_curve_cat   = "Posterior",
               stringsAsFactors = FALSE))
}

#--------------------------------- Script End ---------------------------------#
