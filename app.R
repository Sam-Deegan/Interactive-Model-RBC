################################################################################
## Project: ECON42240 Advanced Macroeconomics                                 ##
## The Real Business Cycle Model: Interactive Shiny App                       ##
################################################################################

## Author:      Sam Deegan
## Affiliation: University College Dublin
## Email:       sam.deegan@ucdconnect.ie

## Usage:
##   Open app.R in RStudio and click Run App, or from this folder:
##     shiny::runApp()
##   Needs R 4.1 or later with shiny, bslib and ggplot2 installed. A hosted
##   copy runs in the browser at https://sam-deegan.com/toy-models/rbc/
##   The stage selector builds the model up one layer at a time:
##     1.1a  how long a shock lasts, at three rates of decay
##     1.5a  forward and backward solutions, and a rational bubble
##     1.5b  roots, stability and the Blanchard and Kahn count
##     2.1a  log-linearising the model, and the error it costs
##     2.1b  the steady state and the ratios it fixes
##     2.1c  simulate the model and hold it against the facts
##     2.4a  prior densities and their supports
##     2.4b  a random-walk Metropolis sampler on a SIMULATED sample
##   Stage 2.4b estimates one parameter of this model from a sample drawn
##   from this model. Nothing in the app is an estimate of a real economy.
##   All text (worked examples, prompts, equations, notation) lives in B_02
##   and B_04. Every figure carries Save PNG and Save PDF at deck size.
##
## Inputs:
##   R/model.R (the solver, simulator, filter and sampler) and R/toolkit.R
##   (shared layout and helpers), both sourced automatically by Shiny.
##
## Outputs:
##   None. Figures are exported from the app through the Save buttons.
##
## Packages:
##   shiny, bslib, ggplot2.
##
## Version:
##   B_03_11_version_chr; history in CHANGELOG.md; git tag vX.Y.Z.
##
## References:
##   Whelan, K. MA Advanced Macroeconomics. Part 6 (solving models with
##     rational expectations), part 7 (the real business cycle model), part
##     10 (estimating DSGE models). Cited in notes as [W6 nn], [W7 nn],
##     [W10 nn], with nn the slide.
##   Whelan's RATS program rbc.prg, for the calibration and the innovation
##     standard deviation sqrt(1 - rho^2).
##   Blanchard, O. and Kahn, C. (1980). The solution of linear difference
##     models under rational expectations. Econometrica 48(5).
##   Cogley, T. and Nason, J. (1995). Output dynamics in real-business-cycle
##     models. American Economic Review 85(3), as quoted on [W7 34].
##   Gali, J. (1999). Technology, employment, and the business cycle.
##     American Economic Review 89(1), as quoted on [W7 38].
##   ECON42240 lecture decks 1.1, 1.5, 2.1 and 2.4, for the stage order and
##     the estimation stage's notation; cited as "deck 2.4".

#-------------------------------- Script Begin --------------------------------#

################################################################################
## A: Table of Contents ########################################################
################################################################################
# Note: C (the model) is R/model.R and T (the toolkit) is R/toolkit.R. This
#   file holds B, D, E, F and G.
#
#   B: Constants
#     B_00  Loading
#     B_01  Parameters
#     B_02  Stages and examples
#     B_03  Presentation
#     B_04  Text
#     B_05  The tests
#   C: Model (R/model.R)
#   T: Toolkit (R/toolkit.R)
#   D: Drawing
#     D_01  Shared furniture
#     D_02  Stage 1.5 (root paths, the two branches)
#     D_03  Stage 2.1 (cycles, technology, growth, responses)
#     D_04  The layers before the simulation (bubble, approximation, ratios)
#     D_05  Stage 2.4 (priors, chains, posterior)
#   E: User interface
#   F: Server
#   G: Run

################################################################################
## B: Constants ################################################################
################################################################################
# Note: Everything a figure or a control needs that is not the model itself.

#### B_00: Loading #############################################################
# Note: Shiny sources R/ before app.R; the two guards cover sourcing app.R by
#   hand, which is how tests/verify_against_rats.R runs it.

###### B_00_01: Load Packages ##################################################
# Note: shiny for the app, bslib for the page frame, ggplot2 for the figures.

library(shiny)
library(bslib)
library(ggplot2)

###### B_00_02: Load the Model #################################################
# Note: R/model.R, guarded so a second source is skipped.

if (!exists("C_01_04_solve_fn")) {
  source(file.path("R", "model.R"))
}

###### B_00_03: Load the Toolkit ###############################################
# Note: R/toolkit.R. B_03_01 reads T_01_02_series_vec at load, so the toolkit
#   has to be in scope first.

if (!exists("T_01_01_palette_vec")) {
  source(file.path("R", "toolkit.R"))
}

###### B_00_04: Global Options #################################################
# Note: No scientific notation; three digits in the console; a fixed seed.

options(scipen = 999, digits = 3)
set.seed(42)

#### B_01: Parameters ##########################################################
# Note: The calibration, the slider ranges and the estimation settings.

###### B_01_01: Defaults #######################################################
# Note: Whelan's calibration (part 7, slide 29, and rbc.prg). Quarterly.

B_01_01_default_lst <- list(
  alpha = 1 / 3,   # capital share
  beta  = 0.99,    # discount factor
  delta = 0.015,   # depreciation, per quarter
  eta   = 1.0,     # CRRA curvature; 1/eta is the intertemporal elasticity
  rho   = 0.95     # persistence of technology
)

###### B_01_02: Slider Ranges ##################################################
# Note: Ranges are wide enough to break the model, deliberately. The
#   diagnostics panel says so when they do, which is the lesson.

B_01_02_range_lst <- list(
  alpha = c(0.20, 0.45, 0.01),
  beta  = c(0.95, 0.999, 0.001),
  delta = c(0.005, 0.040, 0.001),
  eta   = c(0.25, 4.00, 0.25),
  rho   = c(0.00, 1.05, 0.01)
)

###### B_01_03: The Persistence Figure's Roots #################################
# Note: Stage 1.1a's three roots, descending, which is the order the lines
#   read down the panel and take their colours. The names are plotmath, so
#   "rho == 0.9" reaches the panel as the Greek letter. See CONVENTIONS.md 6.

B_01_03_persistence_vec <- c("rho == 0.9" = 0.9,
                             "rho == 0.6" = 0.6,
                             "rho == 0.3" = 0.3)

###### B_01_04: The Estimation Settings ########################################
# Note: Stage 2.4b's four sliders (prior mean and width, proposal step, sample
#   length) and the fixed chain settings. Two chains of 800 draws over 120
#   quarters is 1,600 posterior evaluations, each a solve plus a Kalman pass,
#   which keeps a redraw near a second; the chain length is therefore not a
#   slider. The two starts sit at opposite ends of the unit interval so the
#   burn-in is visible on the trace. The sample is SIMULATED from this model
#   at the sidebar's rho; it is not data.

B_01_04_estimation_lst <- list(
  prior_mean = 0.70,   # prior mean for rho
  prior_sd   = 0.15,   # prior standard deviation for rho
  step       = 0.05,   # proposal standard deviation, the c of deck 2.4
  n_obs      = 120     # quarters in the SIMULATED sample
)

B_01_04_estimation_range_lst <- list(
  prior_mean = c(0.05, 0.95, 0.05),
  prior_sd   = c(0.005, 0.300, 0.005),
  step       = c(0.005, 0.400, 0.005),
  n_obs      = c(40, 200, 20)
)

B_01_04_chain_lst <- list(
  n_draw_int  = 800L,
  warm_up_int = 200L,
  start_vec   = c(0.30, 0.99),
  seed_vec    = c(101L, 202L),
  me_frac_num = 0.50,          # measurement error, as a share of each series
  sigma_z_num = 1,             # standard deviation of technology, per cent
  data_seed   = 20240L
)

###### B_01_05: The Impulse-Response Horizon ###################################
# Note: 50 quarters after impact, the horizon of Whelan's figures [W7 36] and
#   [W7 37], over which investment and hours turn negative.

B_01_05_irf_horizon_int <- 51L

#### B_02: Stages and Examples #################################################
# Note: The stage radio and the worked examples that belong to each stage.

###### B_02_01: Stages #########################################################
# Note: Whelan's own layers: part 7 builds the model in five labelled parts
#   and simulates only in the last, and part 10 rewrites the solved model in
#   state-space form for estimation. The code is the lecture the stage
#   belongs to plus a letter; 2.4a is drawn, 2.4b runs the sampler.

B_02_01_stage_vec <- c(
  "Stage 1: How Long a Shock Lasts"          = "1.1a",
  "Stage 2: Forward, Backward, Bubbles"      = "1.5a",
  "Stage 3: Roots and Stability"             = "1.5b",
  "Stage 4: Log-Linearising the Model"       = "2.1a",
  "Stage 5: The Steady State and Its Ratios" = "2.1b",
  "Stage 6: Simulate, and Check It"          = "2.1c",
  "Stage 7: Priors and Their Supports"       = "2.4a",
  "Stage 8: Sampling the Posterior"          = "2.4b")

###### B_02_01a: Where a Stage Sits in the Order ###############################
# Note: Stage codes are strings, so their order is their position in B_02_01.
#   Returns NA for a code the vector does not carry.

B_02_01a_rank_fn <- function(stage_chr) {
  match(stage_chr, unname(B_02_01_stage_vec))
}

###### B_02_02: Worked Examples ################################################
# Note: One card per stage in the main window. A preset sets the sliders and
#   cannot move the stage; each stage opens on its first example. See
#   CONVENTIONS.md 1 to 3.

B_02_02_example_lst <- list(

  persistence = list(
    label = "Three Rates of Decay",
    stage = "1.1a",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "One surprise hits a series once and nothing hits it again, so the ",
      "path afterwards is y<sub>t</sub> = &rho;y<sub>t&minus;1</sub> and ",
      "nothing else. Persistence ",
      "(&rho;) is the only difference between the three lines, so all three ",
      "start at the same height on the vertical axis - the impact is one ",
      "whatever &rho; is - and they separate only as the horizontal axis runs ",
      "on. What is left of the shock after n periods is &rho;<sup>n</sup>, ",
      "which is a third of it after ten periods at &rho; = 0.9 and ",
      "effectively none of it ",
      "after three at &rho; = 0.3."
    ),
    prompt = paste0(
      "A third of the shock is still there after ten periods at &rho; = 0.9. ",
      "How many periods does &rho; = 0.6 need to fall that far? Read the line ",
      "under the figure before you guess."
    )
  ),

  bubble = list(
    label = "A Bubble on Top",
    stage = "1.5a",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "Take y<sub>t</sub> = aE<sub>t</sub>y<sub>t+1</sub> + x<sub>t</sub> ",
      "and substitute forward. The answer is ",
      "the discounted sum of expected future x, and it is the solution ",
      "economists usually take. It is not the only one. Add any term (b) ",
      "that the market expects to grow at 1/a and the equation is still ",
      "satisfied exactly - that is a rational bubble. The panel draws the ",
      "fundamental and the fundamental plus a bubble from the same x. ",
      "Nothing in the equation prefers one."
    ),
    prompt = paste0(
      "The two paths start almost together. How many periods before the ",
      "difference is obvious? That gap is the whole problem with testing ",
      "for bubbles."
    )
  ),

  approx = list(
    label = "How Good Is the Approximation?",
    stage = "2.1a",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "Log-linearisation rests on one Taylor approximation, exp(x) is about ",
      "1 + x, applied around the steady state. Every equation in the model ",
      "goes through it. The slide says the economy stays close enough for ",
      "that to be safe and does not say how close is close enough, so the ",
      "panel answers it: the error is a tenth of a per cent at a five per ",
      "cent deviation and about nine per cent at fifty."
    ),
    prompt = paste0(
      "Recessions move output a few per cent. Read the error there. Now ",
      "read it at a deviation the size of the Great Depression."
    )
  ),

  ratios = list(
    label = "The Steady State Is Not Free",
    stage = "2.1b",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "The return in steady state has to equal 1/&beta;, and that one ",
      "requirement fixes the ratios: output to capital, investment to ",
      "capital, and investment to output (s<sub>i</sub>). None of them is ",
      "a slider. ",
      "Move depreciation (&delta;) and the investment share moves with it, ",
      "whether or not that was the intention - which is why a calibration ",
      "is a set of joint choices rather than a list of independent ones."
    ),
    prompt = paste0(
      "Whelan's &delta; of 0.015 implies an investment share near 0.20. What ",
      "&delta; would you need for a share of 0.25, and is that a plausible ",
      "rate of depreciation?"
    )
  ),

  baseline = list(
    label = "Whelan's Calibration",
    stage = "2.1c",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "Technology (a) is the only shock, and it is persistent: &rho; is 0.95, ",
      "so a surprise this quarter is still most of the way there a year ",
      "later. A positive surprise raises the marginal product of capital, ",
      "which raises the return (r) and pulls investment (i) up far harder ",
      "than consumption (c), because households smooth. Output (y) follows ",
      "technology almost exactly. Watch how little of the cycle survives ",
      "when you take &rho; down."
    ),
    prompt = paste0(
      "Read the relative standard deviations under the figure. How many ",
      "times more volatile is investment than consumption? Then set &rho; to ",
      "0.5 and read them again."
    )
  ),

  iid = list(
    label = "The Propagation Test",
    stage = "2.1c",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.00),
    story = paste0(
      "Technology is now iid - a surprise this quarter is gone by the next. ",
      "The early claim for these models was that capital accumulation and ",
      "the labour response would turn even iid shocks into cycles, so the ",
      "model would add dynamics of its own. Set &rho; to zero and read the ",
      "correlation of output with technology under the figure. It is one. ",
      "There is no propagation mechanism to find."
    ),
    prompt = paste0(
      "If output is technology whatever &rho; is, what is the model adding? ",
      "Check the Tests tab before you answer."
    )
  ),

  weak = list(
    label = "Weak Propagation",
    stage = "2.1c",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.50),
    story = paste0(
      "The same model with technology far less persistent. The innovation ",
      "standard deviation is &radic;(1 &minus; &rho;<sup>2</sup>), so the ",
      "VARIANCE of technology ",
      "is held at one and only its shape changes: this is not a smaller ",
      "shock, it is a shorter one. The cycles nearly disappear. That is the ",
      "weak-propagation finding - the model has almost no machinery of its ",
      "own, so output dynamics are the dynamics of the driving process."
    ),
    prompt = paste0(
      "Compare the output series with the technology series. Does the model ",
      "add anything to the shock, or does it pass it through?"
    )
  ),

  # eta is the wealth effect in n_t = y_t - eta c_t; with V(N) = aN [W7 14]
  # labour supply is flat at every eta
  elastic = list(
    label = "A Weaker Wealth Effect",
    stage = "2.1c",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 0.5, rho = 0.95),
    story = paste0(
      "Work carries a constant disutility, so labour supply is flat at a ",
      "wage that rises with consumption, and hours (n) obey n<sub>t</sub> ",
      "= y<sub>t</sub> &minus; &eta;",
      "c<sub>t</sub>: they rise with output (y) and fall with ",
      "consumption (c). ",
      "Utility curvature (&eta;), whose inverse is the intertemporal ",
      "elasticity, sets that wealth effect. Lowering it lets the same ",
      "technology surprise draw out ",
      "more hours, so output moves further. Lowering &eta; also raises the ",
      "intertemporal elasticity 1/&eta;, so households also shift consumption ",
      "across time more readily. Estimates of the intertemporal ",
      "elasticity are low, which means &eta; is high: the wrong direction for ",
      "the model."
    ),
    prompt = paste0(
      "Set &eta; to 4 and then to 0.25. Which way does the volatility of ",
      "hours go? Then, at &eta; = 4, move &rho; to 0.99 and read the sign ",
      "of hours on ",
      "impact on the Tests tab."
    )
  ),

  stable = list(
    label = "A Root Inside the Circle",
    stage = "1.5b",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "A single difference equation, y<sub>t</sub> = ",
      "&lambda;y<sub>t&minus;1</sub> + &epsilon;<sub>t</sub>, with ",
      "one surprise and nothing after it. With the root (&lambda;) inside the ",
      "unit circle the effect decays and the series returns to where it ",
      "started. ",
      "With it outside, the same surprise never stops growing. The panel ",
      "draws both from the same shock. Underneath, the root count for the ",
      "full model on the next tab: two predetermined variables, one jump, ",
      "so one root outside is what a unique stable solution needs."
    ),
    prompt = paste0(
      "The explosive path is drawn at 1.05, barely outside the circle. How ",
      "many quarters before it is visibly leaving? Now move &rho; to 1.02 and ",
      "read what the diagnostics say about the full model."
    )
  ),

  supports = list(
    label = "Three Supports, Three Shapes",
    stage = "2.4a",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95),
    story = paste0(
      "A prior has to put probability only where the parameter can go, and ",
      "the support is the first decision rather than the last. A response ",
      "of policy to inflation is free to take either sign, so its prior is ",
      "normal. A standard deviation cannot be negative, so its prior is ",
      "gamma. Persistence (&rho;) has to sit below one or the model has no ",
      "solution at all, so its prior is beta. The beta panel carries two ",
      "curves because one family gives several shapes: the same support, ",
      "two quite different beliefs."
    ),
    prompt = paste0(
      "Put a normal prior on a persistence parameter instead. What share of ",
      "its probability sits on models this app cannot solve? Stage 3 is ",
      "where you can watch one of them."
    )
  ),

  speak = list(
    label = "A Prior the Data Can Move",
    stage = "2.4b",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95,
                  prior_mean = 0.70, prior_sd = 0.15,
                  step = 0.05, n_obs = 120),
    story = paste0(
      "The sample is SIMULATED from this model at the persistence (&rho;) on ",
      "the sidebar, then observed with error, so there is a known answer to ",
      "find. The prior is centred well away from it at 0.70 and is wide. ",
      "Two chains start at opposite ends of the unit interval and walk ",
      "uphill until they meet, which is the burn-in the shaded band covers. ",
      "The posterior that comes out sits near the truth and is far narrower ",
      "than the prior: that narrowing is the data speaking."
    ),
    prompt = paste0(
      "Move the prior mean to 0.40 and back to 0.90. How far does the ",
      "posterior follow it? Then shorten the sample and ask the same ",
      "question again."
    )
  ),

  stuck = list(
    label = "A Prior That Will Not Move",
    stage = "2.4b",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95,
                  prior_mean = 0.60, prior_sd = 0.010,
                  step = 0.015, n_obs = 120),
    story = paste0(
      "The same model and the same sample, with the prior's standard ",
      "deviation cut to 0.01. The posterior now sits on top of the prior ",
      "and the data have changed almost nothing. Nothing failed: a proper ",
      "prior always gives a proper posterior, which is why this case is ",
      "dangerous rather than obvious. A posterior mean of 0.60 would be ",
      "reported as though the sample had produced it, and the sample ",
      "prefers 0.95."
    ),
    prompt = paste0(
      "Read the prior and posterior means under the figure. Now widen the ",
      "prior one step at a time and find the width at which the data start ",
      "to win."
    )
  ),

  coarse = list(
    label = "Steps Too Large",
    stage = "2.4b",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95,
                  prior_mean = 0.70, prior_sd = 0.15,
                  step = 0.40, n_obs = 120),
    story = paste0(
      "The proposal's step size is the one number the algorithm does not ",
      "choose for itself. Set it far too large and almost every candidate ",
      "is worse than where the chain already is, so almost every candidate ",
      "is rejected. The trace becomes a flat line with occasional jumps, ",
      "the acceptance rate under the figure collapses, and the kept draws ",
      "repeat themselves instead of exploring."
    ),
    prompt = paste0(
      "Read the acceptance rate. A common target is roughly one accepted ",
      "proposal in four. How far below it are you, and what does the ",
      "posterior density look like as a result?"
    )
  ),

  crawl = list(
    label = "Steps Too Small",
    stage = "2.4b",
    values = list(alpha = 1 / 3, beta = 0.99, delta = 0.015,
                  eta = 1.0, rho = 0.95,
                  prior_mean = 0.70, prior_sd = 0.15,
                  step = 0.005, n_obs = 120),
    story = paste0(
      "The opposite failure. With a tiny step almost every candidate is ",
      "close enough to be accepted, so the acceptance rate goes up and the ",
      "chain looks healthy - and it crawls. Neither chain reaches the ",
      "posterior inside the warm-up, the two traces never meet, and the ",
      "running means are still drifting at the last iteration. A high ",
      "acceptance rate is not evidence of convergence."
    ),
    prompt = paste0(
      "Compare the two running means at the last iteration. Would you call ",
      "this chain converged if you had only been shown the acceptance rate?"
    )
  )
)

#### B_03: Presentation ########################################################
# Note: Palette, export shapes, file names, controls, help text and headers.

###### B_03_01: Palette ########################################################
# Note: The toolkit's series colours keyed by role: main, second and third
#   for up to three series, compare for a comparator, band for a shaded
#   support, annot for a shaded span, muted for resting-point rules.

B_03_01_palette_vec <- c(
  main    = T_01_02_series_vec[["main"]],
  second  = T_01_02_series_vec[["second"]],
  third   = T_01_02_series_vec[["third"]],
  compare = T_01_02_series_vec[["compare"]],
  band    = T_01_02_series_vec[["band"]],
  annot   = T_01_02_series_vec[["annot"]],
  muted   = T_01_01_palette_vec[["muted"]],
  zero    = T_01_01_zero_chr,
  ground  = T_01_01_palette_vec[["ground"]]
)

###### B_03_02: Export Shapes ##################################################
# Note: The pixel sizes a Save PNG is written at: wide is a full-width slide
#   figure, pair one half of a figure pair. Both are T_02_03c_export_fn's.

B_03_02_shape_lst <- list(
  wide = list(width_px = 1600L, height_px = 800L),
  pair = list(width_px = 1440L, height_px = 720L)
)

###### B_03_03: The Floor on the Aspect Ratio ##################################
# Note: Width over height, at least this on every exported PNG. A panel
#   squarer than 3:2 is height-limited on a 16:9 slide and shrinks.

B_03_03_min_ratio_num <- 1.5

###### B_03_03a: Enforce the Floor at Load #####################################
# Note: Stops the app at load if an export shape is squarer than the floor.

B_03_03a_check_fn <- function(shape_lst, floor_num) {
  for (nm_chr in names(shape_lst)) {
    ratio_num <- shape_lst[[nm_chr]]$width_px / shape_lst[[nm_chr]]$height_px
    if (ratio_num < floor_num) {
      stop(sprintf("Export shape '%s' is %.2f:1, squarer than the %.2f:1 ",
                   nm_chr, ratio_num, floor_num), "floor.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

B_03_03a_ratios_ok_lgl <- B_03_03a_check_fn(B_03_02_shape_lst,
                                            B_03_03_min_ratio_num)

###### B_03_05: Figure File Names ##############################################
# Note: One entry per figure, keyed by plotOutput id: the slug a download
#   takes and the shape it is written at. D_01_03 adds the app and the stage.

B_03_05_figfile_lst <- list(
  plot_persistence = list(name_chr = "persistence",   shape_chr = "wide"),
  plot_roots       = list(name_chr = "root-paths",    shape_chr = "wide"),
  plot_bubble      = list(name_chr = "bubble",        shape_chr = "wide"),
  plot_branch      = list(name_chr = "branches",      shape_chr = "wide"),
  plot_approx      = list(name_chr = "approximation", shape_chr = "wide"),
  plot_ratios      = list(name_chr = "ratios",        shape_chr = "wide"),
  plot_cycles      = list(name_chr = "cycles",        shape_chr = "wide"),
  plot_tech        = list(name_chr = "technology",    shape_chr = "wide"),
  plot_growth      = list(name_chr = "growth",        shape_chr = "wide"),
  plot_irf         = list(name_chr = "responses",     shape_chr = "pair"),
  plot_irf_tech    = list(name_chr = "responses-technology",
                          shape_chr = "pair"),
  plot_priors      = list(name_chr = "priors",        shape_chr = "wide"),
  plot_trace       = list(name_chr = "chains",        shape_chr = "wide"),
  plot_posterior   = list(name_chr = "posterior",     shape_chr = "wide")
)

###### B_03_06: The App's Own Name #############################################
# Note: The first part of every exported file name.

B_03_06_app_chr <- "rbc"

###### B_03_07: Controls #######################################################
# Note: Label, range and step for T_03_01_control_fn, built from B_01_02 and
#   B_01_04. Labels are HTML with the symbol as an entity.

B_03_07_range_fn <- function(label_chr, range_vec) {
  list(label = label_chr, min = range_vec[1L], max = range_vec[2L],
       step = range_vec[3L])
}

B_03_07_controls_lst <- list(
  alpha = B_03_07_range_fn("Capital Share (&alpha;)",
                           B_01_02_range_lst$alpha),
  beta  = B_03_07_range_fn("Discount Factor (&beta;)",
                           B_01_02_range_lst$beta),
  delta = B_03_07_range_fn("Depreciation (&delta;)",
                           B_01_02_range_lst$delta),
  eta   = B_03_07_range_fn("Utility Curvature (&eta;)",
                           B_01_02_range_lst$eta),
  rho   = B_03_07_range_fn("Technology Persistence (&rho;)",
                           B_01_02_range_lst$rho),
  prior_mean = B_03_07_range_fn("Prior Mean for &rho;",
                                B_01_04_estimation_range_lst$prior_mean),
  prior_sd   = B_03_07_range_fn("Prior Standard Deviation for &rho;",
                                B_01_04_estimation_range_lst$prior_sd),
  step       = B_03_07_range_fn("Proposal Step Size (c)",
                                B_01_04_estimation_range_lst$step),
  n_obs      = B_03_07_range_fn("Quarters in the Simulated Sample (T)",
                                B_01_04_estimation_range_lst$n_obs)
)

B_03_07_defaults_lst <- c(B_01_01_default_lst, B_01_04_estimation_lst)

###### B_03_08: What Each Control Does #########################################
# Note: The grey note under each control. One sentence each.

B_03_08_help_lst <- list(
  alpha = "Share of output paid to capital. Whelan sets a third.",
  beta  = paste("Quarterly. Sets the steady-state return, 1/&beta;, and",
                "with it every steady-state ratio."),
  delta = paste("Share of capital lost each quarter. Moves the investment",
                "share, s<sub>i</sub>, with it."),
  eta   = paste("Its inverse is the intertemporal elasticity. It sets how",
                "far higher consumption holds hours back."),
  rho   = paste("How long a technology shock lasts. At one or above the",
                "model has no stable solution."),
  prior_mean = "Where the Beta prior on &rho; is centred.",
  prior_sd   = "How wide it is. Very narrow and the data cannot move it.",
  step       = paste("Standard deviation of each proposed move. Too large",
                     "and almost nothing is accepted; too small and the",
                     "chain crawls."),
  n_obs      = "Length of the SIMULATED sample the sampler is given."
)

###### B_03_09: Figure Card Headers ############################################
# Note: What each card is called on screen. Title Case.

B_03_09_header_lst <- list(
  plot_persistence = "How Long a Shock Lasts",
  plot_bubble      = "The Fundamental and a Bubble",
  plot_roots       = "One Shock, Two Roots",
  plot_branch      = "The Stable Branch and the Other Root",
  plot_approx      = "The Log-Linear Approximation",
  plot_ratios      = "Steady-State Ratios Against Depreciation",
  plot_cycles      = "Simulated Cycles",
  plot_tech        = "Output Against Technology",
  plot_growth      = "Output Growth",
  plot_irf         = "Impulse Responses to a Technology Shock",
  plot_irf_tech    = "Output and Technology After the Shock",
  plot_priors      = "Prior Densities and Their Supports",
  plot_trace       = "The Two Chains",
  plot_posterior   = "Prior Against Posterior"
)

###### B_03_10: QR Code Source #################################################
# Note: The sidebar QR code, from the toolkit.

B_03_10_qr_src_chr <- T_07_04_qr_fn()

###### B_03_11: Version ########################################################
# Note: Semantic version, shown in the footer; CHANGELOG.md has the history.

B_03_11_version_chr <- "1.0.5"

###### B_03_12: Source Repository ##############################################
# Note: The GitHub repo, linked from the footer.

B_03_12_repo_chr <- paste0("https://github.com/Sam-Deegan/",
                        "Interactive-Model-RBC")

#### B_04: Text ################################################################
# Note: Equations, notation, stage guidance and the two scope paragraphs.

###### B_04_01: Equations ######################################################
# Note: The equations panel. One item per equation; "versions" maps the stage
#   a form first applies to its LaTeX and "notes" says what it adds. The
#   panel shows the latest version at or before the stage on screen and flags
#   what is new or changed there. Keys are stage codes, ordered by
#   B_02_01a_rank_fn. Sources: [W6 9] to [W6 14] for the forward equation,
#   the bubble and transversality; [W7 14], [W7 19], [W7 22] to [W7 28] for
#   utility, the approximation, the seven log-linear equations and the
#   steady-state ratios; [W10 11] and [W10 12] for the decision rules and the
#   state-space form; deck 2.4 for the sampler.

B_04_01_equation_lst <- list(

  # --- The model's equations --------------------------------------------------
  list(
    group = "model",
    label = "Autoregressive Process of Order One (AR(1))",
    versions = list(
      "1.1a" = "y_t = \\rho y_{t-1} + \\epsilon_t",
      "2.1a" = "a_t = \\rho a_{t-1} + \\epsilon_t",
      "2.4a" = "z_t = \\rho z_{t-1} + \\epsilon_t"),
    notes = list(
      "1.1a" = paste("One surprise, then nothing, and what survives each",
                     "period is &rho; times what was there before."),
      "2.1a" = paste("The same process is now the model's technology, and",
                     "its only shock."),
      "2.4a" = paste("Lecture 2.4 writes the technology state z where 1.5",
                     "and 2.1 write a. Same process, the deck's letter."))
  ),
  list(
    group = "model", label = "Impulse Response Function (IRF)",
    versions = list("1.1a" = "y_{t+h} = \\rho^{h}"),
    notes = list(
      "1.1a" = paste("The extra y in each period after a one-unit shock,",
                     "with no later shocks. Every path starts at one."))
  ),
  list(
    group = "model", label = "First-Order Stochastic Difference Equation",
    versions = list("1.5a" = "y_t = a\\,E_t y_{t+1} + x_t"),
    notes = list(
      "1.5a" = paste("Today's value depends on what today expects of",
                     "tomorrow, which is what makes the solution forward."))
  ),
  list(
    group = "model", label = "Rational Bubble",
    versions = list("1.5a" = "b_t = a\\,E_t b_{t+1}"),
    notes = list(
      "1.5a" = paste("Any term the market expects to grow at 1/a satisfies",
                     "the equation exactly, on top of the fundamental."))
  ),
  list(
    group = "model", label = "Backward Difference Equation",
    versions = list("1.5b" = "y_t = \\lambda y_{t-1} + \\epsilon_t"),
    notes = list(
      "1.5b" = paste("The same arithmetic run the other way, so the root",
                     "&lambda; decides whether one shock dies or grows."))
  ),
  list(
    group = "model", label = "Resource Constraint",
    versions = list(
      "2.1a" = "y_t = (1 - s_i)\\,c_t + s_i\\,i_t",
      "2.1b" = paste0("y_t = \\Big(1 - \\tfrac{\\alpha\\delta}",
                      "{\\beta^{-1} + \\delta - 1}\\Big) c_t",
                      " + \\tfrac{\\alpha\\delta}",
                      "{\\beta^{-1} + \\delta - 1}\\,i_t")),
    notes = list(
      "2.1a" = paste("Output is split between consumption and investment,",
                     "in the shares the steady state implies."),
      "2.1b" = paste("The investment share is substituted out, so the",
                     "equation is in the parameters alone."))
  ),
  list(
    group = "model", label = "Capital Accumulation",
    versions = list("2.1a" = "k_t = (1 - \\delta)\\,k_{t-1} + \\delta\\,i_t"),
    notes = list(
      "2.1a" = paste("The first of the two predetermined variables: capital",
                     "this period was decided last period."))
  ),
  list(
    group = "model", label = "Production Function (Cobb-Douglas)",
    versions = list(
      "2.1a" = "y_t = a_t + \\alpha\\,k_{t-1} + (1 - \\alpha)\\,n_t"),
    notes = list(
      "2.1a" = paste("Cobb-Douglas in logs, so the shares are the",
                     "coefficients and technology enters additively."))
  ),
  list(
    # wage equals the marginal rate of substitution a C^eta [W7 14]; logs
    # give the line the solver uses
    group = "model", label = "Labour Condition (Intratemporal)",
    versions = list("2.1a" = paste0(
      "(1 - \\alpha)\\tfrac{Y_t}{N_t} = a\\,C_t^{\\eta}",
      " \\implies n_t = y_t - \\eta\\,c_t")),
    notes = list(
      "2.1a" = paste("The wage equals the marginal rate of substitution.",
                     "With V(N) linear the right side has no hours in it,",
                     "so given consumption the wage is fixed."))
  ),
  list(
    group = "model", label = "Euler Equation",
    versions = list(
      "2.1a" = "c_t = E_t c_{t+1} - \\tfrac{1}{\\eta}\\,E_t r_{t+1}"),
    notes = list(
      "2.1a" = paste("The only forward-looking equation, and the reason",
                     "consumption is the jump variable."))
  ),
  list(
    group = "model", label = "Return on Capital",
    versions = list(
      "2.1a" = "r_t = \\mu\\,(y_t - k_{t-1})",
      "2.1b" = "r_t = \\big(1 - \\beta(1 - \\delta)\\big)(y_t - k_{t-1})"),
    notes = list(
      "2.1a" = paste("The return rises with output per unit of capital,",
                     "which is what pulls investment after a shock."),
      "2.1b" = paste("The return coefficient is substituted out, leaving",
                     "the discount factor and depreciation."))
  ),

  # --- Assumptions ------------------------------------------------------------
  list(
    group = "assumption", label = "Transversality Condition",
    versions = list("1.5a" = "\\lim_{N \\to \\infty} a^{N} E_t y_{t+N} = 0"),
    notes = list(
      "1.5a" = paste("Imposing it is what picks the fundamental out of the",
                     "family of solutions. It is not algebra."))
  ),
  list(
    group = "assumption", label = "Stability Condition",
    versions = list(
      "1.5a" = "|a| < 1",
      "1.5b" = "|\\lambda| < 1"),
    notes = list(
      "1.5a" = paste("Stated for the purely forward-looking case. Lecture",
                     "1.5 warns that this is not the general condition."),
      "1.5b" = paste("With a lag in the equation the root that matters is",
                     "&lambda;, and the same threshold decides."))
  ),
  list(
    group = "assumption", label = "Blanchard and Kahn Condition",
    versions = list("1.5b" = "m = q"),
    notes = list(
      "1.5b" = paste("As many roots outside the unit circle as jump",
                     "variables gives one stable solution, and no more."))
  ),
  list(
    group = "assumption", label = "Period Utility (CRRA, Separable Labour)",
    versions = list("2.1a" = paste0(
      "U(C_t) - V(N_t) = \\tfrac{C_t^{1-\\eta}}{1-\\eta} - a\\,N_t")),
    notes = list(
      "2.1a" = paste("Whelan's [W7 14]. &eta; is the curvature of utility in",
                     "consumption, and 1/&eta; the intertemporal elasticity",
                     "in the Euler equation."))
  ),
  list(
    group = "assumption", label = "Log-Linear Approximation",
    versions = list("2.1a" = "e^{x} \\approx 1 + x"),
    notes = list(
      "2.1a" = paste("Every equation in the model goes through this one",
                     "Taylor approximation around the steady state, and x",
                     "is each variable's log deviation from it."))
  ),
  list(
    group = "assumption", label = "Prior Support",
    versions = list(
      "2.4a" = paste0("0 < \\rho < 1 \\implies \\rho \\sim ",
                      "\\mathrm{Beta}(s_1, s_2)")),
    notes = list(
      "2.4a" = paste("The prior puts no mass where the model has no",
                     "solution, so the support is the first decision."))
  ),

  # --- Solved forms -----------------------------------------------------------
  list(
    group = "solved", label = "Forward Solution (the Fundamental)",
    versions = list(
      "1.5a" = "y_t = \\sum_{k=0}^{\\infty} a^{k}\\,E_t x_{t+k}"),
    notes = list(
      "1.5a" = paste("Substituting forward without limit leaves the",
                     "discounted sum of what the driving variable is",
                     "expected to be."))
  ),
  list(
    group = "solved", label = "Solution With a Bubble",
    versions = list(
      "1.5a" = "y_t = \\sum_{k=0}^{\\infty} a^{k}\\,E_t x_{t+k} + b_t"),
    notes = list(
      "1.5a" = paste("The fundamental plus any bubble solves the equation",
                     "too. Nothing in the equation prefers one."))
  ),
  list(
    group = "solved", label = "Capital's Decision Rule",
    versions = list("2.4a" = "k_t = a_{kk}\\,k_{t-1} + a_{kz}\\,z_t"),
    notes = list(
      "2.4a" = paste("The solver returns this. The two coefficients are",
                     "functions of the structural parameters, not free",
                     "numbers."))
  ),
  list(
    group = "solved", label = "Consumption's Decision Rule",
    versions = list("2.4a" = "c_t = a_{ck}\\,k_{t-1} + a_{cz}\\,z_t"),
    notes = list(
      "2.4a" = paste("Consumption follows the same two states, with its",
                     "own coefficients out of the same solution."))
  ),
  list(
    group = "model", label = "Observation Equations",
    versions = list(
      "2.4a" = paste0("k^{*}_t = k_t + u^{k}_t, \\quad ",
                      "c^{*}_t = c_t + u^{c}_t")),
    notes = list(
      "2.4a" = paste("What the data show is the model's own variable plus",
                     "an error, and the two errors are separate."))
  ),
  list(
    group = "solved", label = "Transition Equation",
    versions = list(
      "2.4a" = paste0("\\begin{pmatrix} k_{t-1} \\\\ z_t \\end{pmatrix}",
                      " = \\begin{pmatrix} a_{kk} & a_{kz} \\\\ 0 & \\rho",
                      " \\end{pmatrix}",
                      "\\begin{pmatrix} k_{t-2} \\\\ z_{t-1}",
                      " \\end{pmatrix}",
                      " + \\begin{pmatrix} 0 \\\\ \\epsilon_t",
                      " \\end{pmatrix}")),
    notes = list(
      "2.4a" = paste("Only the second row carries a shock, so the state",
                     "shock covariance has rank one. The filter runs",
                     "anyway."))
  ),
  list(
    group = "solved", label = "Measurement Equation",
    versions = list(
      "2.4a" = paste0("\\begin{pmatrix} k^{*}_{t-1} \\\\ c^{*}_t",
                      " \\end{pmatrix}",
                      " = \\begin{pmatrix} 1 & 0 \\\\ a_{ck} & a_{cz}",
                      " \\end{pmatrix}",
                      "\\begin{pmatrix} k_{t-1} \\\\ z_t \\end{pmatrix}",
                      " + \\begin{pmatrix} u^{k}_{t-1} \\\\ u^{c}_t",
                      " \\end{pmatrix}")),
    notes = list(
      "2.4a" = paste("Two measurement errors and one shock give three",
                     "random terms for two observables, which is what the",
                     "counting rule asks for."))
  ),
  list(
    group = "solved", label = "Log-Likelihood",
    versions = list(
      "2.4a" = paste0("\\log L(\\theta) = -\\tfrac{1}{2}",
                      "\\sum_{t=1}^{T}\\big[\\log|\\Omega_t|",
                      " + \\varepsilon_t' \\Omega_t^{-1}\\varepsilon_t",
                      "\\big] + \\mathrm{constant}")),
    notes = list(
      "2.4a" = paste("Lecture 1.4's filter, returning one forecast error",
                     "and its variance each period. What is new is that",
                     "&theta; is now structural."))
  ),
  list(
    group = "solved", label = "Posterior",
    versions = list(
      "2.4b" = paste0("p(\\theta \\mid Z^{o}) \\propto ",
                      "L(Z^{o} \\mid \\theta)\\,p(\\theta)")),
    notes = list(
      "2.4b" = paste("The prior turns the likelihood into a distribution",
                     "over the parameters. Its constant is unknown and",
                     "never needed."))
  ),
  list(
    group = "solved", label = "Random-Walk Proposal",
    versions = list(
      "2.4b" = paste0("\\theta^{*} = \\theta^{(s)} + \\nu, \\quad ",
                      "\\nu \\sim \\mathcal{N}(0, c^{2})")),
    notes = list(
      "2.4b" = paste("The step size c is the one number the algorithm does",
                     "not choose for itself, which is why it is a slider."))
  ),
  list(
    group = "solved", label = "Acceptance Rule",
    versions = list(
      "2.4b" = paste0("\\alpha = \\min\\Big\\{1, ",
                      "\\tfrac{L(Z^{o} \\mid \\theta^{*})\\,",
                      "p(\\theta^{*})}",
                      "{L(Z^{o} \\mid \\theta^{(s)})\\,",
                      "p(\\theta^{(s)})}\\Big\\}")),
    notes = list(
      "2.4b" = paste("The symmetric-proposal form. The posterior's unknown",
                     "constant cancels because it appears above and below."))
  ),

  # --- Thresholds and simplifications -----------------------------------------
  list(
    group = "descriptor", label = "Investment Share",
    versions = list(
      "2.1b" = paste0("s_i = \\tfrac{\\alpha\\delta}",
                      "{\\delta - 1 + \\beta^{-1}}")),
    notes = list(
      "2.1b" = paste("Implied by the others rather than chosen, so moving",
                     "depreciation moves the investment share with it."))
  ),
  list(
    group = "descriptor", label = "Return Coefficient",
    versions = list("2.1b" = "\\mu = 1 - \\beta(1 - \\delta)"),
    notes = list(
      "2.1b" = paste("What the steady-state return requirement leaves in",
                     "front of output per unit of capital."))
  ),
  list(
    group = "descriptor", label = "Innovation Standard Deviation",
    versions = list("2.1c" = "\\sigma_{\\epsilon} = \\sqrt{1 - \\rho^{2}}"),
    notes = list(
      "2.1c" = paste("Held so that the variance of technology is one at",
                     "every persistence. Lowering &rho; gives a shorter",
                     "shock, not a smaller one."))
  )
)

###### B_04_02: Equation Group Titles ##########################################
# Note: The four headed blocks the equations panel sorts into.

B_04_02_groups_vec <- c(
  model      = "Model Equations",
  assumption = "Assumptions",
  solved     = "Solved Forms",
  descriptor = "Thresholds and Simplifications"
)

###### B_04_03: Notation #######################################################
# Note: One entry per symbol: "grp" (var, par, tgt, shk) picks the column,
#   "from" is the stage the symbol first appears, "txt" the gloss. Two
#   letters carry two meanings and are listed twice: alpha (capital share,
#   then the acceptance probability) and a (the forward weight, then
#   technology). The tests check the key and the equations close over each
#   other.

B_04_03_notation_lst <- list(

  # Lecture 1.1
  list(grp = "var", sym = "y_t",
       txt = "the series a difference equation determines", from = "1.1a"),
  list(grp = "shk", sym = "\\epsilon_t",
       txt = "the shock, later the technology innovation", from = "1.1a"),
  list(grp = "par", sym = "\\rho", txt = "the persistence parameter",
       from = "1.1a"),
  list(grp = "var", sym = "h", txt = "periods elapsed since the shock",
       from = "1.1a"),

  # Lecture 1.5
  list(grp = "var", sym = "x_t", txt = "the driving variable, taken as given",
       from = "1.5a"),
  list(grp = "par", sym = "a",
       txt = "the forward weight; in 2.1 the disutility of an hour",
       from = "1.5a"),
  list(grp = "var", sym = "E_t",
       txt = "the expectation using all information at t", from = "1.5a"),
  list(grp = "var", sym = "b_t", txt = "the bubble component", from = "1.5a"),
  list(grp = "var", sym = "N", txt = "the horizon the substitution is taken to",
       from = "1.5a"),
  list(grp = "var", sym = "k", txt = "the index summed over", from = "1.5a"),
  list(grp = "par", sym = "\\lambda",
       txt = "a root of the characteristic equation", from = "1.5b"),
  list(grp = "tgt", sym = "m", txt = "roots outside the unit circle, counted",
       from = "1.5b"),
  list(grp = "tgt", sym = "q", txt = "jump variables, counted", from = "1.5b"),

  # Lecture 2.1
  list(grp = "var", sym = "x",
       txt = "a log deviation from steady state (per cent)", from = "2.1a"),
  list(grp = "var", sym = "a_t", txt = "technology", from = "2.1a"),
  list(grp = "var", sym = "c_t", txt = "consumption", from = "2.1a"),
  list(grp = "var", sym = "i_t", txt = "investment", from = "2.1a"),
  list(grp = "var", sym = "k_t",
       txt = "the capital stock at the end of the period", from = "2.1a"),
  list(grp = "var", sym = "n_t", txt = "hours worked", from = "2.1a"),
  list(grp = "var", sym = "r_t", txt = "the return on capital", from = "2.1a"),
  list(grp = "par", sym = "\\alpha",
       txt = "the share of output paid to capital", from = "2.1a"),
  list(grp = "par", sym = "\\beta", txt = "the household's discount factor",
       from = "2.1a"),
  list(grp = "par", sym = "\\delta",
       txt = "the fraction of capital lost each period", from = "2.1a"),
  list(grp = "par", sym = "\\eta",
       txt = paste("the curvature of utility in consumption (CRRA); 1/&eta; is",
                   "the intertemporal elasticity of substitution"),
       from = "2.1a"),
  list(grp = "var", sym = "Y_t", txt = "output, in levels", from = "2.1a"),
  list(grp = "var", sym = "C_t", txt = "consumption, in levels",
       from = "2.1a"),
  list(grp = "var", sym = "N_t", txt = "hours worked, in levels",
       from = "2.1a"),
  list(grp = "var", sym = "U", txt = "utility from consumption",
       from = "2.1a"),
  list(grp = "var", sym = "V", txt = "disutility of hours", from = "2.1a"),
  list(grp = "tgt", sym = "s_i", txt = "the investment share of output",
       from = "2.1a"),
  list(grp = "tgt", sym = "\\mu", txt = "the return coefficient",
       from = "2.1a"),
  list(grp = "par", sym = "\\sigma_{\\epsilon}",
       txt = "the standard deviation of the technology innovation",
       from = "2.1c"),

  # Lecture 2.4
  list(grp = "var", sym = "z_t",
       txt = "the technology state, lecture 2.4's letter for a", from = "2.4a"),
  list(grp = "par", sym = "a_{kk}",
       txt = "capital's loading on its own lag", from = "2.4a"),
  list(grp = "par", sym = "a_{kz}",
       txt = "capital's loading on the technology state", from = "2.4a"),
  list(grp = "par", sym = "a_{ck}",
       txt = "consumption's loading on lagged capital", from = "2.4a"),
  list(grp = "par", sym = "a_{cz}",
       txt = "consumption's loading on the technology state", from = "2.4a"),
  list(grp = "var", sym = "k^{*}_t", txt = "capital as measured",
       from = "2.4a"),
  list(grp = "var", sym = "c^{*}_t", txt = "consumption as measured",
       from = "2.4a"),
  list(grp = "shk", sym = "u^{k}_t",
       txt = "the measurement error in capital", from = "2.4a"),
  list(grp = "shk", sym = "u^{c}_t",
       txt = "the measurement error in consumption", from = "2.4a"),
  list(grp = "var", sym = "Z^{o}",
       txt = "the observed part of the model's variables", from = "2.4a"),
  list(grp = "var", sym = "L", txt = "the likelihood of the sample",
       from = "2.4a"),
  list(grp = "var", sym = "\\varepsilon_t",
       txt = "the filter's forecast error", from = "2.4a"),
  list(grp = "var", sym = "\\Omega_t",
       txt = "the variance of the forecast error", from = "2.4a"),
  list(grp = "par", sym = "\\theta", txt = "the structural parameters",
       from = "2.4a"),
  list(grp = "par", sym = "T", txt = "quarters in the simulated sample",
       from = "2.4a"),
  list(grp = "par", sym = "s_1",
       txt = "the beta prior's first shape parameter", from = "2.4a"),
  list(grp = "par", sym = "s_2",
       txt = "the beta prior's second shape parameter", from = "2.4a"),
  list(grp = "var", sym = "p(\\theta)", txt = "the prior density",
       from = "2.4a"),
  list(grp = "var", sym = "\\theta^{*}",
       txt = "the proposed parameter vector", from = "2.4b"),
  list(grp = "var", sym = "\\theta^{(s)}",
       txt = "the chain's parameter vector at draw s", from = "2.4b"),
  list(grp = "shk", sym = "\\nu", txt = "the proposal innovation",
       from = "2.4b"),
  list(grp = "par", sym = "c", txt = "the proposal step size", from = "2.4b"),
  list(grp = "tgt", sym = "\\alpha",
       txt = "the acceptance probability, not 2.1's capital share",
       from = "2.4b")
)

###### B_04_04: Stage Guidance #################################################
# Note: The Note line under the prompt, one per stage.

B_04_04_guidance_lst <- list(
  "1.1a" = paste0(
    "The same &rho; is the technology persistence slider on the sidebar. It ",
    "does not move this figure, which is drawn at three fixed values so it ",
    "matches the slide; stage 6 is where &rho; drives a whole model."
  ),
  "1.5a" = paste0(
    "The bubble grows at 1/a, so a below one means it explodes. This is the ",
    "root question in its simplest setting, met before any matrix."
  ),
  "2.1a" = paste0(
    "Every equation in the model is put through this one approximation. ",
    "The error is what it costs."
  ),
  "2.1b" = paste0(
    "Move &delta; on the sidebar and watch the investment share move with it. ",
    "s<sub>i</sub> is reported on the Diagnostics tab."
  ),
  "1.5b" = paste0(
    "Load 'A root inside the circle'. The two paths come from one shock and ",
    "differ only in &lambda;. Then move &rho; above one on the sidebar and ",
    "watch ",
    "the diagnostics below refuse the model: two jump-side roots where there ",
    "should be one."
  ),
  "2.1c" = paste0(
    "Load 'Whelan's calibration' and read the three simulated figures in ",
    "order: the cycles, output against technology, and growth. Then set &rho; ",
    "to 0.5 and read them again. The pair at the bottom is Whelan's impulse ",
    "responses. The faded curve is where each series was before you moved ",
    "anything."
  ),
  "2.4a" = paste0(
    "Nothing on this panel is estimated and nothing on it is data. The ",
    "shading under each curve is the range the parameter is allowed to ",
    "take, and the curve reaches the axis at the edge of it on its own."
  ),
  "2.4b" = paste0(
    "The sidebar's &rho; is the TRUTH the sample is drawn at, not a guess at ",
    "it. The four controls below the rule are the estimation's own: the ",
    "prior's mean and width, the proposal's step size, and how many ",
    "quarters the sampler is given. Each redraw reruns both chains."
  )
)

###### B_04_05: Scope ##########################################################
# Note: The sidebar note saying what the app is not: the full DSGE of part 11
#   and the accelerator of part 12 are not built, and the absent frictions
#   are the stage 2.1c argument.

B_04_05_scope_chr <- paste0(
  "Whelan's Part 7, and Part 10's estimation of it. One technology shock, ",
  "no money, no frictions, markets clear. The full DSGE - Smets and ",
  "Wouters, seven shocks and the frictions this model leaves out - is Part ",
  "11, and Part 12 builds the financial accelerator on it. Neither is in ",
  "this app. What the missing frictions cost you is the whole point of the ",
  "scorecard."
)

###### B_04_06: What Stage 2.4 Is, and Is Not ##################################
# Note: Printed on stage 2.4b: the sample is simulated from this model and
#   nothing is an estimate of a real economy. The figure subtitles say the
#   same, so an exported PNG carries it too.

B_04_06_simulated_chr <- paste0(
  "The sample is SIMULATED from this model at the parameters on the ",
  "sidebar and then observed with error. It is not data, and nothing here ",
  "is an estimate of any real economy. The point of a known truth is that ",
  "the sampler can be held against it, which no real estimation allows."
)

#### B_05: The Tests ###########################################################
# Note: The Tests tab's scorecard.

###### B_05_01: What the Model Is Held Against #################################
# Note: Whelan runs the simulated model against four facts and it fails three
#   ([W7 31] to [W7 38]). The data column is his number where he gives one:
#   Cogley and Nason's 0.34 as quoted on [W7 34], Gali's hours sign as quoted
#   on [W7 38], and the iid propagation claim of [W7 32]. The model side is
#   computed live.

B_05_01_test_lst <- list(
  list(name = "Growth autocorrelation",
       data = "+0.34",
       src  = "Cogley and Nason (1995), [W7 34]",
       pass = function(m_lst, irf_df) m_lst$growth_acf > 0.15),
  list(name = "Output against technology",
       data = "cycles are more than the shock",
       src  = "[W7 32], the propagation claim",
       pass = function(m_lst, irf_df) m_lst$tech_corr < 0.95),
  list(name = "Hours after a positive technology shock",
       data = "fall",
       src  = "Gali (1999), [W7 38]",
       pass = function(m_lst, irf_df) irf_df$mod_hours_val[1] < 0),
  list(name = "Investment against consumption",
       data = "investment far more volatile",
       src  = "the one it gets right",
       pass = function(m_lst, irf_df) {
         m_lst$sd_rel_vec[["investment"]] >
           2 * m_lst$sd_rel_vec[["consumption"]]
       })
)

################################################################################
## D: Drawing ##################################################################
################################################################################
# Note: Every builder takes a solved model and a parameter list and returns a
#   ggplot; none touches a Shiny input. All use T_02_01_theme_fn.

#### D_01: Shared Furniture ####################################################
# Note: Export, ghost, end-of-line names, file names and the on-screen form.

###### D_01_01: Export a Figure ################################################
# Note: Writes a figure at deck size through T_02_03c_export_fn, which strips
#   the title and subtitle, drops any legend and paints a white ground. A
#   "pair" shape is one half of a figure pair; format is "png" or "pdf".

D_01_01_export_fn <- function(plot_obj, path_chr, shape_chr = "wide",
                              format_chr = "png") {
  T_02_03c_export_fn(path_chr, plot_obj,
                     pair = identical(shape_chr, "pair"),
                     format = format_chr)
  invisible(path_chr)
}

###### D_01_02: Ghost Layer ####################################################
# Note: The loaded example's figure drawn first at the toolkit's opacity, and
#   skipped when the sliders agree with it. The cycles figure takes none:
#   three ghosted series on one axis read as noise. See CONVENTIONS.md 8.

D_01_02_ghost_alpha_num <- T_01_04_ghost_alpha_num

D_01_02_ghost_off_fn <- function(par_lst, ref_lst) {
  if (is.null(ref_lst)) return(TRUE)
  isTRUE(all.equal(par_lst[names(B_01_01_default_lst)],
                   ref_lst[names(B_01_01_default_lst)],
                   tolerance = 1e-12))
}

###### D_01_02a: Spreading the Names at the End of the Lines ###################
# Note: Pushes end-of-line names apart from the bottom up, each lifted only
#   far enough to clear the one below by gap_num. See CONVENTIONS.md 6.

D_01_02a_spread_fn <- function(y_vec, gap_num) {
  ord_int <- order(y_vec)
  run_vec <- y_vec[ord_int]
  for (i in seq_along(run_vec)[-1L]) {
    run_vec[i] <- max(run_vec[i], run_vec[i - 1L] + gap_num)
  }
  out_vec <- y_vec
  out_vec[ord_int] <- run_vec
  out_vec
}

###### D_01_02b: The Name at the End of Each Line ##############################
# Note: Each series' last point, spread, as a geom_text layer. The caller adds
#   the x expansion that makes room for the names.

D_01_02b_endname_fn <- function(plot_df, x_chr, y_chr, series_chr,
                                gap_frac_num = 0.075, size_num = 3.5,
                                parse_lgl = FALSE) {
  end_df <- do.call(rbind, lapply(unique(plot_df[[series_chr]]),
    function(one_chr) {
      one_df <- plot_df[plot_df[[series_chr]] == one_chr, , drop = FALSE]
      one_df[nrow(one_df), , drop = FALSE]
    }))
  end_df$mod_name_val <- D_01_02a_spread_fn(
    end_df[[y_chr]], gap_frac_num * diff(range(plot_df[[y_chr]])))
  geom_text(
    data = end_df,
    aes(x = .data[[x_chr]], y = .data$mod_name_val,
        label = .data[[series_chr]], colour = .data[[series_chr]]),
    parse = parse_lgl, hjust = 0, vjust = 0.5, size = size_num,
    nudge_x = 0.015 * diff(range(plot_df[[x_chr]])),
    show.legend = FALSE, inherit.aes = FALSE)
}

###### D_01_03: The Name a Download Takes ######################################
# Note: {app}-{stage}-{figure}.{ext}, with the stage read at download time.

D_01_03_figfile_fn <- function(id_chr, stage_chr, ext_chr = "png") {
  paste0(D_01_03a_figstem_fn(id_chr, stage_chr), ".", ext_chr)
}

###### D_01_03a: The Stem of That Name #########################################
# Note: The same name without its extension, for T_07_07h_exports_fn.

D_01_03a_figstem_fn <- function(id_chr, stage_chr) {
  paste0(B_03_06_app_chr, "-", stage_chr, "-",
         B_03_05_figfile_lst[[id_chr]]$name_chr)
}

###### D_01_04: A Figure as a Half-Width Card Shows It #########################
# Note: On screen the plot title is dropped (the card header names the
#   figure) and the subtitle becomes the caption, which T_02_01c_draw_fn
#   lifts into the note under the card. A title carrying a live verdict is
#   kept as the caption's first sentence. Drawn at 96 dpi.

D_01_04_res_int <- 96L

D_01_04_screen_fn <- function(plot_obj, title_to_caption = FALSE) {
  if (!inherits(plot_obj, "ggplot")) return(plot_obj)
  cap_vec <- character(0)
  ttl_chr <- plot_obj$labels$title
  if (title_to_caption && is.character(ttl_chr) && nzchar(ttl_chr)) {
    cap_vec <- c(cap_vec, paste0(ttl_chr, "."))
  }
  sub_chr <- plot_obj$labels$subtitle
  if (is.character(sub_chr) && length(sub_chr) == 1L && nzchar(sub_chr)) {
    cap_vec <- c(cap_vec, sub_chr)
  }
  plot_obj$labels$title    <- NULL
  plot_obj$labels$subtitle <- NULL
  if (length(cap_vec) > 0L) {
    plot_obj$labels$caption <- paste(cap_vec, collapse = " ")
  }
  T_02_01c_draw_fn(plot_obj)
}

#### D_02: Stage 1.5 ###########################################################
# Note: The scalar root paths and the two branches of the full model.

###### D_02_01: Root Paths #####################################################
# Note: One builder for stages 1.1a and 1.5b: the path of a scalar shock, one
#   root per line, up to three roots. Colours follow the order the roots are
#   given; each root is named at the end of its line as plotmath, with no
#   legend. See CONVENTIONS.md 6.

D_02_01_line_vec <- c(B_03_01_palette_vec[["main"]],
                      B_03_01_palette_vec[["second"]],
                      B_03_01_palette_vec[["third"]])

D_02_01_roots_fn <- function(
    lambda_vec = c("lambda == 0.85" = 0.85,
                   "lambda == 1.05" = 1.05),
    n_periods_int    = 40L,
    first_period_int = 0L,
    shock_period_int = 0L,
    title_chr        = "One Shock, Two Roots",
    x_title          = expression(bold("Quarters Since the Shock (" * h *
                                         ")")),
    y_title          = expression(bold("Deviation (" * y[t] * ")")),
    linetype_vec     = "solid",
    name_on_line_lgl = TRUE,
    lambda_stable_num    = NULL,
    lambda_explosive_num = NULL) {

  # a stable-and-explosive pair given by name becomes the vector
  if (!is.null(lambda_stable_num) || !is.null(lambda_explosive_num)) {
    lambda_vec <- c(
      if (is.null(lambda_stable_num)) 0.85 else lambda_stable_num,
      if (is.null(lambda_explosive_num)) 1.05 else lambda_explosive_num)
    names(lambda_vec) <- sprintf("lambda == %s", format(lambda_vec))
  }

  path_df <- C_03_01_roots_fn(lambda_vec, n_periods_int,
                              first_period_int, shock_period_int)
  lab_vec <- unique(path_df$mod_case_cat)
  n_line  <- length(lab_vec)

  colour_vec <- stats::setNames(
    rep_len(D_02_01_line_vec, n_line), lab_vec)
  line_vec   <- stats::setNames(
    rep_len(linetype_vec, n_line), lab_vec)

  # names spread from the bottom up, since decayed paths end together; the
  # zero rule is clipped to the data so it does not run under the names
  name_lyr <- NULL
  x_scale  <- NULL
  zero_lyr <- geom_hline(yintercept = 0, linetype = "dashed",
                         linewidth = 0.35,
                         colour = B_03_01_palette_vec[["zero"]])
  if (isTRUE(name_on_line_lgl)) {
    end_df <- do.call(rbind, lapply(lab_vec, function(one_chr) {
      one_df <- path_df[path_df$mod_case_cat == one_chr, , drop = FALSE]
      one_df[nrow(one_df), , drop = FALSE]
    }))
    end_df$mod_name_val <- D_01_02a_spread_fn(
      end_df$mod_path_val, 0.095 * diff(range(path_df$mod_path_val)))

    name_lyr <- geom_text(
      data = end_df,
      aes(x = mod_period_tm, y = mod_name_val, label = mod_case_cat,
          colour = mod_case_cat),
      parse = TRUE, hjust = 0, vjust = 0.5, size = 3.5,
      nudge_x = 0.02 * diff(range(path_df$mod_period_tm)),
      show.legend = FALSE, inherit.aes = FALSE)
    x_scale <- scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.14)))
    zero_lyr <- annotate(
      "segment", x = min(path_df$mod_period_tm),
      xend = max(path_df$mod_period_tm), y = 0, yend = 0,
      linetype = "dashed", linewidth = 0.35,
      colour = B_03_01_palette_vec[["zero"]])
  }

  ggplot(path_df, aes(x = mod_period_tm, y = mod_path_val,
                      colour = mod_case_cat, linetype = mod_case_cat)) +
    zero_lyr +
    geom_line(linewidth = 0.7) +
    name_lyr +
    x_scale +
    scale_colour_manual(values = colour_vec) +
    scale_linetype_manual(values = line_vec) +
    (if (isTRUE(name_on_line_lgl)) {
      guides(colour = "none", linetype = "none")
    }) +
    labs( x = x_title, y = y_title) +
    T_02_01_theme_fn(grid = "h")
}

###### D_02_02: The Two Branches ###############################################
# Note: The model solved twice from one shock, imposing transversality and
#   taking the other root. Both satisfy the structural equations; the y axis
#   is logarithmic because the divergent path would flatten the stable one.

D_02_02_branch_fn <- function(par, n_horizon_int = 40L) {

  one_fn <- function(branch_chr, label_chr) {
    sol_lst <- C_01_04_solve_fn(par, branch_chr = branch_chr)
    if (!isTRUE(sol_lst$converged_lgl)) return(NULL)
    irf_df <- C_02_02_irf_fn(sol_lst, n_horizon_int = n_horizon_int)
    data.frame(mod_horizon_tm = irf_df$mod_horizon_tm,
               mod_value_val  = abs(irf_df$mod_output_val),
               mod_branch_cat = label_chr)
  }

  plot_df <- rbind(one_fn("stable",   "Transversality"),
                   one_fn("unstable", "Other root"))
  if (is.null(plot_df) || nrow(plot_df) == 0L) return(NULL)

  ggplot(plot_df, aes(x = mod_horizon_tm, y = mod_value_val,
                      colour = mod_branch_cat, linetype = mod_branch_cat)) +
    geom_line(linewidth = 0.7) +
    # no spread: the gap is in data units and this axis is logarithmic
    D_01_02b_endname_fn(plot_df, "mod_horizon_tm", "mod_value_val",
                        "mod_branch_cat", gap_frac_num = 0) +
    scale_y_log10() +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.36))) +
    # the other root is the comparator, solid in the compare colour
    scale_colour_manual(values = c(
      "Transversality" = B_03_01_palette_vec[["main"]],
      "Other root"     = B_03_01_palette_vec[["compare"]])) +
    scale_linetype_manual(values = c(
      "Transversality" = "solid",
      "Other root"     = "solid")) +
    guides(colour = "none", linetype = "none") +
    labs(
      x = expression(bold("Quarters After the Shock (" * h * ")")),
      y = expression(bold("Output (" * group("|", y[t], "|") * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

#### D_03: Stage 2.1 ###########################################################
# Note: The simulated series and the impulse responses, Whelan's own figures.

###### D_03_01: Simulated Cycles ###############################################
# Note: Whelan's simulated output, consumption and investment [W7 30] on one
#   panel, each named at the end of its own line.

D_03_01_cycles_fn <- function(sim_df, n_show_int = 200L) {

  if (is.null(sim_df)) return(NULL)
  keep_df <- sim_df[sim_df$mod_period_tm <= n_show_int, , drop = FALSE]

  long_fn <- function(one_df) {
    do.call(rbind, lapply(
      c("output", "consumption", "investment"), function(v_chr) {
        data.frame(mod_period_tm = one_df$mod_period_tm,
                   mod_value_val = one_df[[paste0("mod_", v_chr, "_val")]],
                   mod_series_cat = tools::toTitleCase(v_chr))
      }))
  }
  plot_df <- long_fn(keep_df)

  # the zero rule is clipped to the data, clear of the names in the margin
  ggplot(plot_df, aes(x = mod_period_tm, y = mod_value_val,
                      colour = mod_series_cat)) +
    annotate("segment", x = min(plot_df$mod_period_tm),
             xend = max(plot_df$mod_period_tm), y = 0, yend = 0,
             linetype = "dashed", linewidth = 0.35,
             colour = B_03_01_palette_vec[["zero"]]) +
    geom_line(linewidth = 0.6) +
    # a wider gap and margin than the default, for a half-width card
    D_01_02b_endname_fn(plot_df, "mod_period_tm", "mod_value_val",
                        "mod_series_cat", gap_frac_num = 0.14) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.36))) +
    scale_colour_manual(values = c(
      Output      = B_03_01_palette_vec[["main"]],
      Consumption = B_03_01_palette_vec[["second"]],
      Investment  = B_03_01_palette_vec[["third"]])) +
    guides(colour = "none") +
    labs(
      x = expression(bold("Quarter (" * t * ")")),
      y = expression(bold("Per Cent Deviation"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

###### D_03_02: Output Against Technology ######################################
# Note: Whelan's output against technology [W7 32], with a ghost.

D_03_02_technology_fn <- function(sim_df, ref_df = NULL, n_show_int = 200L) {

  if (is.null(sim_df)) return(NULL)
  keep_df <- sim_df[sim_df$mod_period_tm <= n_show_int, , drop = FALSE]

  long_fn <- function(one_df) {
    rbind(
      data.frame(mod_period_tm = one_df$mod_period_tm,
                 mod_value_val = one_df$mod_output_val,
                 mod_series_cat = "Output"),
      data.frame(mod_period_tm = one_df$mod_period_tm,
                 mod_value_val = one_df$mod_technology_val,
                 mod_series_cat = "Technology"))
  }
  plot_df <- long_fn(keep_df)

  ghost_lyr <- NULL
  if (!is.null(ref_df)) {
    ghost_df  <- long_fn(ref_df[ref_df$mod_period_tm <= n_show_int, ,
                                drop = FALSE])
    ghost_lyr <- geom_line(
      data = ghost_df,
      aes(x = mod_period_tm, y = mod_value_val, colour = mod_series_cat),
      alpha = D_01_02_ghost_alpha_num, linewidth = 0.6, inherit.aes = FALSE
    )
  }

  ggplot(plot_df, aes(x = mod_period_tm, y = mod_value_val,
                      colour = mod_series_cat)) +
    annotate("segment", x = min(plot_df$mod_period_tm),
             xend = max(plot_df$mod_period_tm), y = 0, yend = 0,
             linetype = "dashed", linewidth = 0.35,
             colour = B_03_01_palette_vec[["zero"]]) +
    ghost_lyr +
    geom_line(linewidth = 0.6) +
    D_01_02b_endname_fn(plot_df, "mod_period_tm", "mod_value_val",
                        "mod_series_cat", gap_frac_num = 0.14) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    scale_colour_manual(values = c(
      Output     = B_03_01_palette_vec[["main"]],
      Technology = B_03_01_palette_vec[["second"]])) +
    guides(colour = "none") +
    labs(
      x = expression(bold("Quarter (" * t * ")")),
      y = expression(bold("Per Cent Deviation"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

###### D_03_03: Growth #########################################################
# Note: Whelan's output growth [W7 34]. The autocorrelation is printed under
#   the figure, not inside it.

D_03_03_growth_fn <- function(sim_df, n_show_int = 200L) {

  if (is.null(sim_df)) return(NULL)
  keep_df <- sim_df[sim_df$mod_period_tm <= n_show_int &
                      !is.na(sim_df$mod_growth_val), , drop = FALSE]

  ggplot(keep_df, aes(x = mod_period_tm, y = mod_growth_val)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35,
               colour = B_03_01_palette_vec[["zero"]]) +
    geom_line(linewidth = 0.5, colour = B_03_01_palette_vec[["main"]]) +
    labs(
      x = expression(bold("Quarter (" * t * ")")),
      y = expression(bold("Output Growth (" * Delta * y[t] * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

###### D_03_04: Impulse Responses, as a Pair ###################################
# Note: [W7 36] and [W7 37] as two panels on one shared y scale: output,
#   consumption and hours on the left, output against technology on the
#   right. Investment is left to the cycles figure, where its scale does not
#   flatten the others. At Whelan's rho the output response peaks on impact.

D_03_04_irf_long_fn <- function(one_df, var_vec) {
  do.call(rbind, lapply(var_vec, function(v_chr) {
    data.frame(mod_horizon_tm = one_df$mod_horizon_tm,
               mod_value_val  = one_df[[paste0("mod_", v_chr, "_val")]],
               mod_series_cat = tools::toTitleCase(v_chr),
               stringsAsFactors = FALSE)
  }))
}

# the shared y range of the pair, live and ghost
D_03_04_irf_range_fn <- function(irf_df, ref_df = NULL) {
  var_vec <- c("output", "consumption", "hours", "technology")
  val_vec <- unlist(lapply(list(irf_df, ref_df), function(one_df) {
    if (is.null(one_df)) return(NULL)
    unlist(one_df[paste0("mod_", var_vec, "_val")], use.names = FALSE)
  }), use.names = FALSE)
  range(c(0, val_vec))
}

D_03_04_irf_panel_fn <- function(irf_df, ref_df, var_vec, colour_vec,
                                 title_chr) {

  if (is.null(irf_df)) return(NULL)

  plot_df <- D_03_04_irf_long_fn(irf_df, var_vec)
  ylim_vec <- D_03_04_irf_range_fn(irf_df, ref_df)

  ghost_lyr <- NULL
  if (!is.null(ref_df)) {
    ghost_lyr <- geom_line(
      data = D_03_04_irf_long_fn(ref_df, var_vec),
      aes(x = mod_horizon_tm, y = mod_value_val, colour = mod_series_cat),
      alpha = D_01_02_ghost_alpha_num, linewidth = 0.6, inherit.aes = FALSE
    )
  }

  # the gap is a share of the shared range, so both panels spread alike
  name_lyr <- D_01_02b_endname_fn(plot_df, "mod_horizon_tm",
                                  "mod_value_val", "mod_series_cat",
                                  gap_frac_num = 0.12 *
                                    diff(ylim_vec) /
                                    max(diff(range(plot_df$mod_value_val)),
                                        1e-09))

  ggplot(plot_df, aes(x = mod_horizon_tm, y = mod_value_val,
                      colour = mod_series_cat)) +
    annotate("segment", x = min(plot_df$mod_horizon_tm),
             xend = max(plot_df$mod_horizon_tm), y = 0, yend = 0,
             linetype = "dashed", linewidth = 0.35,
             colour = B_03_01_palette_vec[["zero"]]) +
    ghost_lyr +
    geom_line(linewidth = 0.7) +
    name_lyr +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.30)),
                       breaks = pretty(range(plot_df$mod_horizon_tm), 5L)) +
    scale_y_continuous(limits = ylim_vec,
                       expand = expansion(mult = c(0.06, 0.06))) +
    scale_colour_manual(values = colour_vec) +
    guides(colour = "none") +
    labs(
      x = expression(bold("Quarters After the Shock (" * h * ")")),
      y = expression(bold("Per Cent Deviation"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

D_03_04_irf_fn <- function(irf_df, ref_df = NULL) {
  D_03_04_irf_panel_fn(
    irf_df, ref_df, c("output", "consumption", "hours"),
    c(Output      = B_03_01_palette_vec[["main"]],
      Consumption = B_03_01_palette_vec[["second"]],
      Hours       = B_03_01_palette_vec[["third"]]),
    "Responses to a Technology Shock")
}

###### D_03_05: Output Against Technology, the Pair's Right Panel ##############
# Note: [W7 37]. The gap between the two responses is all the propagation the
#   model adds.

D_03_05_irf_tech_fn <- function(irf_df, ref_df = NULL) {
  D_03_04_irf_panel_fn(
    irf_df, ref_df, c("output", "technology"),
    c(Output     = B_03_01_palette_vec[["main"]],
      Technology = B_03_01_palette_vec[["second"]]),
    "Output Against Technology")
}

#### D_04: The Layers Before the Simulation ####################################
# Note: One figure per layer Whelan builds before he simulates; each is drawn,
#   not estimated.

###### D_04_01: Forward, Backward and Bubbles ##################################
# Note: [W6 13]. Two paths satisfying one equation: the fundamental is flat
#   because x is constant, and everything above it is bubble.

D_04_01_bubble_fn <- function(a_num = 0.95, n_periods_int = 40L) {

  plot_df <- C_05_03_bubble_fn(a_num = a_num, n_periods_int = n_periods_int)

  ggplot(plot_df, aes(x = mod_period_tm, y = mod_value_val,
                      colour = mod_path_cat, linetype = mod_path_cat)) +
    geom_line(linewidth = 0.7) +
    D_01_02b_endname_fn(plot_df, "mod_period_tm", "mod_value_val",
                        "mod_path_cat") +
    # a wide margin for long names; breaks from the data so the margin
    # earns no tick
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.38)),
                       breaks = pretty(range(plot_df$mod_period_tm), 4L)) +
    scale_colour_manual(values = c(
      "Fundamental"   = B_03_01_palette_vec[["main"]],
      "With a bubble" = B_03_01_palette_vec[["second"]])) +
    scale_linetype_manual(values = c(
      "Fundamental"   = "solid",
      "With a bubble" = "solid")) +
    guides(colour = "none", linetype = "none") +
    labs(
      x = expression(bold("Period (" * t * ")")),
      y = expression(bold("Price (" * y[t] * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

###### D_04_02: The Log-Linearisation Error ####################################
# Note: [W7 19]. The exact curve and the linear approximation on one panel;
#   the error is the gap between them.

D_04_02_approx_fn <- function(max_dev_num = 0.6) {

  raw_df  <- C_05_01_approx_fn(max_dev_num)
  plot_df <- rbind(
    data.frame(mod_deviation_val = raw_df$mod_deviation_val,
               mod_value_val = raw_df$mod_exact_val,
               mod_series_cat = "Exact"),
    data.frame(mod_deviation_val = raw_df$mod_deviation_val,
               mod_value_val = raw_df$mod_linear_val,
               mod_series_cat = "Linear"))

  # the axis is in per cent, and the end-of-line name needs a data column
  plot_df$mod_pct_val <- 100 * plot_df$mod_deviation_val

  ggplot(plot_df, aes(x = mod_pct_val, y = mod_value_val,
                      colour = mod_series_cat, linetype = mod_series_cat)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35,
               colour = B_03_01_palette_vec[["zero"]]) +
    geom_line(linewidth = 0.7) +
    D_01_02b_endname_fn(plot_df, "mod_pct_val", "mod_value_val",
                        "mod_series_cat") +
    # breaks from the data, so the margin earns no tick
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.34)),
                       breaks = pretty(range(plot_df$mod_pct_val), 5L)) +
    scale_colour_manual(values = c(
      "Exact"  = B_03_01_palette_vec[["main"]],
      "Linear" = B_03_01_palette_vec[["second"]])) +
    scale_linetype_manual(values = c(
      "Exact"  = "solid",
      "Linear" = "solid")) +
    guides(colour = "none", linetype = "none") +
    labs(
      x = expression(bold("Log Deviation, Per Cent (" * 100 * x * ")")),
      y = expression(bold("Level (" * e^x * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

###### D_04_03: Steady-State Ratios ############################################
# Note: [W7 26] and [W7 27]. The investment share against the swept
#   parameter, with the current setting marked on the curve.

D_04_03_ratios_fn <- function(par, sweep_chr = "delta") {

  plot_df <- C_05_02_ratios_fn(par, sweep_chr)
  here_num <- par[[sweep_chr]]
  now_num  <- C_01_01_derived_fn(par)$s_i
  x_title  <- if (identical(sweep_chr, "delta")) {
    expression(bold("Depreciation (" * delta * ")"))
  } else {
    expression(bold("Discount factor (" * beta * ")"))
  }

  ggplot(plot_df, aes(x = mod_sweep_val, y = mod_i_over_y_val)) +
    geom_line(linewidth = 0.7, colour = B_03_01_palette_vec[["main"]]) +
    geom_vline(xintercept = here_num, linetype = "dotted", linewidth = 0.45,
               colour = B_03_01_palette_vec[["muted"]]) +
    # the marker sits on the one series, so it takes that series' colour
    geom_point(x = here_num, y = now_num, size = 2.2,
               colour = B_03_01_palette_vec[["main"]]) +
    # the current share is named as a tick on the opposite axis
    scale_y_continuous(sec.axis = dup_axis(
      breaks = now_num, labels = sprintf("%.3f", now_num), name = NULL)) +
    labs(
      x = x_title,
      y = expression(bold("Investment Share (" * s[i] * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

#### D_05: Stage 2.4 ###########################################################
# Note: The prior densities, the two chains and the posterior. The chains and
#   the posterior rest on a sample SIMULATED from this model, and both say
#   so in their subtitle so an exported PNG carries it.

###### D_05_01: Prior Densities ################################################
# Note: Whelan's three prior families ([W10 17] to [W10 19]), one facet each
#   on its own scales. The range drawn is wider than the support, so each
#   curve lands on the axis at its own boundary; the support is shaded.

D_05_01_priors_fn <- function(n_point_int = 601L) {

  plot_df  <- C_06_02_priors_fn(n_point_int)
  inside_df <- plot_df[plot_df$mod_inside_lgl, , drop = FALSE]

  name_df <- plot_df[!duplicated(plot_df$mod_curve_cat),
                     c("mod_curve_cat", "mod_family_cat",
                       "mod_mode_val", "mod_peak_val")]

  # Beta(2, 2)'s name runs leftwards from its peak, clear of Beta(5, 2)'s
  # flank; Gamma's runs rightwards down its tail
  name_df$mod_hjust_val <- ifelse(
    name_df$mod_curve_cat == "Beta(2, 2)", 1, 0.5)
  name_df$mod_hjust_val[name_df$mod_curve_cat == "Gamma(2, 4)"] <- 0.3
  name_df$mod_label_chr <- ifelse(name_df$mod_curve_cat == "Beta(2, 2)",
                                  "Beta\n(2, 2)", name_df$mod_curve_cat)

  one_df <- data.frame(mod_value_val = 1,
                       mod_family_cat = factor("Beta",
                         levels = levels(plot_df$mod_family_cat)))

  # one series per facet takes the first colour; the beta facet's second
  # curve takes the second
  colour_vec <- c("Normal(0, 1)" = B_03_01_palette_vec[["main"]],
                  "Gamma(2, 4)"  = B_03_01_palette_vec[["main"]],
                  "Beta(2, 2)"   = B_03_01_palette_vec[["main"]],
                  "Beta(5, 2)"   = B_03_01_palette_vec[["second"]])

  ggplot(plot_df, aes(x = mod_value_val, y = mod_density_val,
                      colour = mod_curve_cat, group = mod_curve_cat)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35,
               colour = B_03_01_palette_vec[["zero"]]) +
    geom_vline(data = one_df, aes(xintercept = mod_value_val),
               linetype = "dotted", linewidth = 0.45,
               colour = B_03_01_palette_vec[["muted"]]) +
    # position = "identity": geom_area stacks by default
    geom_area(data = inside_df,
              aes(x = mod_value_val, y = mod_density_val,
                  group = mod_curve_cat),
              position = "identity",
              fill = B_03_01_palette_vec[["band"]], alpha = 0.30,
              inherit.aes = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_text(data = name_df,
              aes(x = mod_mode_val, y = mod_peak_val,
                  label = mod_label_chr, colour = mod_curve_cat,
                  hjust = mod_hjust_val),
              vjust = -0.4, size = 3.4, lineheight = 0.9,
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~ mod_family_cat, nrow = 1L, scales = "free") +
    scale_colour_manual(values = colour_vec) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.20))) +
    guides(colour = "none") +
    labs(
      x = expression(bold("Parameter Value (" * theta * ")")),
      y = expression(bold("Prior Density (" * p(theta) * ")"))
    ) +
    T_02_01_theme_fn(grid = "none")
}

###### D_05_02: The Chain's Path Through One Parameter #########################
# Note: Two chains against iteration, with the warm-up shaded, beside their
#   running means on the same rho axis. Each chain is named at the end of
#   its line with its start in brackets. The title is computed: it says the
#   chains have met only when their post-warm-up means agree.

D_05_02_trace_fn <- function(chain_lst) {

  if (is.null(chain_lst)) return(NULL)

  panel_vec <- c("The Draws", "The Running Mean")
  long_fn <- function(col_chr, panel_chr) {
    data.frame(mod_draw_tm   = chain_lst$draw_df$mod_draw_tm,
               mod_value_val = chain_lst$draw_df[[col_chr]],
               mod_chain_cat = chain_lst$draw_df$mod_chain_cat,
               mod_panel_cat = factor(panel_chr, levels = panel_vec),
               stringsAsFactors = FALSE)
  }
  plot_df <- rbind(long_fn("mod_rho_val", panel_vec[1L]),
                   long_fn("mod_running_val", panel_vec[2L]))

  # names spread from the bottom up within each panel; the gap suits a
  # two-line label on a half-width card
  gap_num <- 0.30 * diff(range(plot_df$mod_value_val))
  end_df  <- do.call(rbind, lapply(panel_vec, function(p_chr) {
    one_df <- plot_df[plot_df$mod_panel_cat == p_chr, , drop = FALSE]
    last_df <- do.call(rbind, lapply(unique(one_df$mod_chain_cat),
      function(c_chr) {
        two_df <- one_df[one_df$mod_chain_cat == c_chr, , drop = FALSE]
        two_df[nrow(two_df), , drop = FALSE]
      }))
    last_df$mod_name_val <- D_01_02a_spread_fn(last_df$mod_value_val, gap_num)
    last_df
  }))

  end_df$mod_label_chr <- sub(" \\(", "\n(", end_df$mod_chain_cat)

  band_df <- data.frame(
    mod_panel_cat = factor(panel_vec[1L], levels = panel_vec))

  colour_vec <- stats::setNames(
    c(B_03_01_palette_vec[["main"]], B_03_01_palette_vec[["second"]]),
    unique(chain_lst$draw_df$mod_chain_cat))

  # post-warm-up means held against the larger of the two standard
  # deviations, a cheap between-against-within comparison
  kept_df  <- chain_lst$draw_df[
    chain_lst$draw_df$mod_draw_tm > chain_lst$warm_up_int, , drop = FALSE]
  mean_vec <- tapply(kept_df$mod_rho_val, kept_df$mod_chain_cat, mean)
  sd_vec   <- tapply(kept_df$mod_rho_val, kept_df$mod_chain_cat, stats::sd)
  met_lgl  <- diff(range(mean_vec)) < 0.75 * max(sd_vec, 1e-09)

  title_chr <- if (isTRUE(met_lgl)) {
    "Two Chains, Started Far Apart, Settling on One Level"
  } else {
    "Two Chains, Started Far Apart, That Have Not Met"
  }

  ggplot(plot_df, aes(x = mod_draw_tm, y = mod_value_val,
                      colour = mod_chain_cat)) +
    annotate("rect", xmin = 0.5, xmax = chain_lst$warm_up_int + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = B_03_01_palette_vec[["annot"]], alpha = 0.12) +
    # the truth rule is clipped to the draws, clear of the names
    annotate("segment", x = 1, xend = chain_lst$n_draw_int,
             y = chain_lst$truth_num, yend = chain_lst$truth_num,
             linetype = "dotted", linewidth = 0.45,
             colour = B_03_01_palette_vec[["muted"]]) +
    geom_line(linewidth = 0.45) +
    geom_text(data = band_df,
              aes(x = chain_lst$warm_up_int * 0.04, y = -Inf),
              label = "Warm-up", vjust = -0.8, hjust = 0, size = 3.2,
              colour = B_03_01_palette_vec[["muted"]],
              inherit.aes = FALSE) +
    geom_text(data = end_df,
              aes(x = mod_draw_tm, y = mod_name_val, label = mod_label_chr,
                  colour = mod_chain_cat),
              hjust = 0, vjust = 0.5, size = 3.0, lineheight = 0.95,
              nudge_x = 0.03 * chain_lst$n_draw_int,
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~ mod_panel_cat, nrow = 1L) +
    # breaks from the chain length, so the margin earns no tick
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.46)),
                       breaks = pretty(c(0, chain_lst$n_draw_int), 4L)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    scale_colour_manual(values = colour_vec) +
    guides(colour = "none") +
    labs(
      x = expression(bold("Iteration (s)")),
      y = expression(bold("Persistence (" * rho * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

###### D_05_03: Prior Against Posterior ########################################
# Note: The prior as drawn and a kernel density of the kept draws, each named
#   at its peak. The two means are dotted rules from the axis to the curve,
#   named as breaks on the primary axis; a default break near a mean is
#   dropped. See CONVENTIONS.md 6.

D_05_03_posterior_fn <- function(dens_df, prior_lst, chain_lst) {

  if (is.null(dens_df) || is.null(chain_lst)) return(NULL)

  # prior drawn first, so the posterior sits on top where they coincide
  dens_df$mod_curve_cat <- factor(dens_df$mod_curve_cat,
                                  levels = c("Prior", "Posterior"))

  prior_mean_num <- prior_lst$mean_num
  post_mean_num  <- chain_lst$post_mean_num

  peak_fn <- function(one_chr) {
    one_df <- dens_df[dens_df$mod_curve_cat == one_chr, , drop = FALSE]
    one_df[which.max(one_df$mod_density_val), , drop = FALSE]
  }
  name_df <- rbind(peak_fn("Prior"), peak_fn("Posterior"))

  # the names are spread apart, since a tight prior's two curves coincide
  name_df$mod_name_val <- D_01_02a_spread_fn(
    name_df$mod_density_val, 0.09 * max(dens_df$mod_density_val))

  at_fn <- function(one_chr, x_num) {
    one_df <- dens_df[dens_df$mod_curve_cat == one_chr, , drop = FALSE]
    one_df$mod_density_val[which.min(abs(one_df$mod_value_val - x_num))]
  }
  rule_df <- data.frame(
    mod_value_val   = c(prior_mean_num, post_mean_num),
    mod_density_val = c(at_fn("Prior", prior_mean_num),
                        at_fn("Posterior", post_mean_num)))

  # named breaks on the primary axis; two means that coincide become one
  # break
  same_lgl <- abs(prior_mean_num - post_mean_num) < 0.035
  mean_vec <- if (same_lgl) prior_mean_num else
    c(prior_mean_num, post_mean_num)
  mean_lab <- if (same_lgl) {
    sprintf("%.2f\nprior, posterior", prior_mean_num)
  } else {
    c(sprintf("%.2f\nprior", prior_mean_num),
      sprintf("%.2f\nposterior", post_mean_num))
  }
  base_vec <- c(0.25, 0.50, 0.75)
  keep_lgl <- vapply(base_vec, function(b_num) {
    min(abs(b_num - mean_vec)) > 0.12
  }, logical(1))
  # an end break is dropped when a mean lands within 0.1 of it
  end_vec   <- c(0, 1)[vapply(c(0, 1), function(e_num) {
    min(abs(e_num - mean_vec)) > 0.10
  }, logical(1))]
  break_vec <- c(end_vec, base_vec[keep_lgl], mean_vec)
  label_vec <- c(sprintf("%.2f", end_vec), sprintf("%.2f", base_vec[keep_lgl]),
                 mean_lab)
  ord_int   <- order(break_vec)

  ggplot(dens_df, aes(x = mod_value_val, y = mod_density_val,
                      colour = mod_curve_cat, linetype = mod_curve_cat)) +
    geom_segment(data = rule_df,
                 aes(x = mod_value_val, xend = mod_value_val,
                     y = 0, yend = mod_density_val),
                 linetype = "dotted", linewidth = 0.45,
                 colour = B_03_01_palette_vec[["muted"]],
                 inherit.aes = FALSE) +
    geom_line(linewidth = 0.7) +
    geom_text(data = name_df,
              aes(x = mod_value_val, y = mod_name_val,
                  label = mod_curve_cat, colour = mod_curve_cat),
              hjust = 0.5, vjust = -0.5, size = 3.5,
              show.legend = FALSE, inherit.aes = FALSE) +
    scale_x_continuous(breaks = break_vec[ord_int],
                       labels = label_vec[ord_int],
                       limits = c(0, 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
    # the prior is the comparator, solid in the compare colour
    scale_colour_manual(values = c(
      Prior     = B_03_01_palette_vec[["compare"]],
      Posterior = B_03_01_palette_vec[["main"]])) +
    scale_linetype_manual(values = c(Prior = "solid",
                                     Posterior = "solid")) +
    guides(colour = "none", linetype = "none") +
    labs(
      x = expression(bold("Persistence (" * rho * ")")),
      y = expression(bold("Density (" * p(rho) * ")"))
    ) +
    T_02_01_theme_fn(grid = "h")
}

################################################################################
## E: User Interface ###########################################################
################################################################################
# Note: The toolkit's page frame: nav bar, title, sidebar of slider-plus-box
#   controls, equations card, worked-example card, prompt, tiles, figure
#   cards and footer. The sidebar chooses the model; the main window chooses
#   what to run in it. See CONVENTIONS.md 1.

#### E_01: Sidebar #############################################################
# Note: The stage radio, the scope note, the controls, Reset and the QR code.

###### E_01_01: Control Shorthand ##############################################
# Note: One slider-plus-box control from B_03_07 and B_03_08.

E_01_01_ctl_fn <- function(id_chr) {
  T_03_01_control_fn(id_chr, B_03_07_controls_lst, B_03_08_help_lst,
                     B_03_07_defaults_lst)
}

###### E_01_02: Sidebar ########################################################
# Note: The five model controls, then stage 2.4b's four estimation controls,
#   shown only where they do something. On 2.4b the persistence slider above
#   them becomes the truth the sample is drawn at.

E_01_02_sidebar_lst <- sidebar(
  width = 380,
  radioButtons("stage", "Stage of the Model", choices = B_02_01_stage_vec,
               selected = "2.1c"),
  T_03_05_note_fn(B_04_05_scope_chr),
  accordion(
    open = c("The Model", "The Estimation"),
    accordion_panel(
      "The Model",
      E_01_01_ctl_fn("alpha"),
      E_01_01_ctl_fn("beta"),
      E_01_01_ctl_fn("delta"),
      E_01_01_ctl_fn("eta"),
      E_01_01_ctl_fn("rho")
    ),
    accordion_panel(
      "The Estimation",
      conditionalPanel(
        "input.stage == '2.4b'",
        E_01_01_ctl_fn("prior_mean"),
        E_01_01_ctl_fn("prior_sd"),
        E_01_01_ctl_fn("step"),
        E_01_01_ctl_fn("n_obs")
      ),
      conditionalPanel(
        "input.stage != '2.4b'",
        tags$p(class = "stat-caption",
               "These four controls act at stage 8, Sampling the Posterior.")
      )
    )
  ),
  actionButton("reset", "Reset Everything",
               class = "btn-outline-secondary btn-sm w-100"),
  T_07_10b_sidebarqr_fn(B_03_10_qr_src_chr)
)

#### E_02: Main Panel ##########################################################
# Note: The figure cards, the rows they sit in and the page itself.

###### E_02_01: One Figure Card ################################################
# Note: The toolkit's card, titled from B_03_09.

E_02_01_fig_fn <- function(id_chr) {
  T_07_07f_figcard_fn(id_chr, HTML(B_03_09_header_lst[[id_chr]]))
}

###### E_02_02: A Row of Figures for One Stage #################################
# Note: Two cards to a row, stacked on a narrow screen, shown only at their
#   stage. A lone figure takes the left column so it is the same size as a
#   paired one.

E_02_02_row_fn <- function(...) {
  card_lst <- list(...)
  if (length(card_lst) == 1L) card_lst <- c(card_lst, list(tags$div()))
  do.call(T_07_07g_pair_fn, card_lst)
}

E_02_02_stage_fn <- function(stage_chr, ...) {
  conditionalPanel(sprintf("input.stage == '%s'", stage_chr), ...)
}

###### E_02_03: The Figures for Each Stage #####################################
# Note: The rows each stage shows; stage 2.4b carries the simulated-sample
#   paragraph above its pair.

E_02_03_figures_lst <- tagList(
  E_02_02_stage_fn("1.1a", E_02_02_row_fn(E_02_01_fig_fn("plot_persistence"))),
  E_02_02_stage_fn("1.5a", E_02_02_row_fn(E_02_01_fig_fn("plot_bubble"))),
  E_02_02_stage_fn("1.5b", E_02_02_row_fn(E_02_01_fig_fn("plot_roots"),
                                          E_02_01_fig_fn("plot_branch"))),
  E_02_02_stage_fn("2.1a", E_02_02_row_fn(E_02_01_fig_fn("plot_approx"))),
  E_02_02_stage_fn("2.1b", E_02_02_row_fn(E_02_01_fig_fn("plot_ratios"))),
  E_02_02_stage_fn(
    "2.1c",
    E_02_02_row_fn(E_02_01_fig_fn("plot_cycles"),
                   E_02_01_fig_fn("plot_tech")),
    E_02_02_row_fn(E_02_01_fig_fn("plot_irf"),
                   E_02_01_fig_fn("plot_irf_tech")),
    E_02_02_row_fn(E_02_01_fig_fn("plot_growth"))
  ),
  E_02_02_stage_fn("2.4a", E_02_02_row_fn(E_02_01_fig_fn("plot_priors"))),
  E_02_02_stage_fn(
    "2.4b",
    tags$div(class = "narrative",
             tags$div(class = "nar-head", "A Simulated Sample, Not Data"),
             tags$p(B_04_06_simulated_chr)),
    E_02_02_row_fn(E_02_01_fig_fn("plot_trace"),
                   E_02_01_fig_fn("plot_posterior"))
  )
)

###### E_02_04: Page ###########################################################
# Note: The whole page: nav bar, sidebar, equation tabs, presets, prompt,
#   problems, tiles, figures and footer.

E_02_04_ui <- tagList(
  T_07_08b_nav_fn(),
  page_sidebar(
    title        = T_07_09_title_fn("The RBC Model, Solved and Simulated",
                                    B_03_10_qr_src_chr),
    window_title = paste("The RBC Model -", T_07_01_author_chr),
    fillable     = FALSE,
    theme        = T_07_05_theme_fn(),
    sidebar      = E_01_02_sidebar_lst,
    T_07_08_head_fn(),
    tags$head(
      tags$style(HTML(T_05_07_preset_css_chr)),
      tags$script(HTML(T_05_05_preset_js_chr))
    ),
    T_07_07j_eqtabs_fn(
      title = textOutput("eq_title", inline = TRUE),
      nav_panel("Equations", uiOutput("equation_ui")),
      nav_panel("Notation", uiOutput("notation_ui")),
      nav_panel("In Words", uiOutput("explain_ui")),
      nav_panel("Tests", uiOutput("test_ui")),
      nav_panel("Diagnostics", uiOutput("diagnostic_ui"))
    ),
    T_05_04_presets_fn(B_02_02_example_lst, B_02_01_stage_vec,
                       stage_word = ""),
    uiOutput("prompt"),
    uiOutput("problems"),
    uiOutput("tiles"),
    E_02_03_figures_lst,
    T_07_11_footer_fn(paste(
      "Notation follows Whelan's Parts 6, 7 and 10 and lectures 1.1 to 2.4.",
      "The sample in stages 7 and 8 is simulated.",
      paste0("Version ", B_03_11_version_chr, ".")), repo = B_03_12_repo_chr)
  )
)

################################################################################
## F: Server ###################################################################
################################################################################
# Note: One function. Presets and sliders (F_01), the solved model live and
#   ghost (F_02), the text outputs (F_03), the figures and their downloads
#   (F_04), the Tests and Diagnostics tabs (F_05), the equation tabs (F_06).

#### F_01: The Server Function #################################################
# Note: Everything below runs inside it.

F_01_01_server <- function(input, output, session) {

  ###### F_01_00: Captions and Controls ########################################
  # Note: The lifted-caption outputs, then the slider-plus-box pairs.

  T_07_07d_cap_fn(output)
  T_03_02_sync_fn(input, session, B_03_07_controls_lst)

  val <- function(id_chr) T_03_04_val_fn(input, id_chr)

  set_values_fn <- function(values_lst) {
    for (nm in names(values_lst)) {
      if (!is.null(B_03_07_controls_lst[[nm]])) {
        T_03_03_set_fn(session, B_03_07_controls_lst, nm, values_lst[[nm]])
      }
    }
    invisible(NULL)
  }

  ###### F_01_01: Assemble the Parameter List ##################################
  # Note: Called for the live run and for the ghost.

  F_01_01_assemble_fn <- function(values_lst) {
    list(alpha = values_lst$alpha, beta = values_lst$beta,
         delta = values_lst$delta, eta = values_lst$eta, rho = values_lst$rho)
  }

  F_01_02_ref_rv      <- reactiveVal(NULL)
  F_01_02_scenario_rv <- reactiveVal(NULL)

  ###### F_01_03: Presets ######################################################
  # Note: Each stage opens on its first example; a preset sets the sliders and
  #   never moves the stage. The loaded button is marked through dgPreset.

  F_01_03_stage_examples_fn <- function(stage_chr) {
    keep_lgl <- vapply(B_02_02_example_lst,
                       function(e) identical(e$stage, stage_chr), logical(1))
    B_02_02_example_lst[keep_lgl]
  }

  F_01_04_load_fn <- function(key_chr) {
    one_lst <- B_02_02_example_lst[[key_chr]]
    if (is.null(one_lst)) return(invisible(NULL))
    set_values_fn(one_lst$values)
    F_01_02_ref_rv(one_lst$values)
    F_01_02_scenario_rv(key_chr)
    session$sendCustomMessage("dgPreset", key_chr)
    invisible(NULL)
  }

  observeEvent(input$stage, {
    ex_lst <- F_01_03_stage_examples_fn(input$stage)
    if (length(ex_lst) > 0L) F_01_04_load_fn(names(ex_lst)[1L])
  }, ignoreInit = FALSE)

  lapply(names(B_02_02_example_lst), function(k_chr) {
    observeEvent(input[[paste0("preset_", k_chr)]], {
      F_01_04_load_fn(k_chr)
    }, ignoreInit = TRUE)
  })

  observeEvent(input$reset, {
    # the four estimation controls reset with the rest
    set_values_fn(B_03_07_defaults_lst)
    F_01_02_ref_rv(NULL)
    F_01_02_scenario_rv(NULL)
    session$sendCustomMessage("dgPreset", "")
  })

  F_01_05_scn_r <- reactive({
    k_chr <- F_01_02_scenario_rv()
    if (is.null(k_chr)) return(NULL)
    B_02_02_example_lst[[k_chr]]
  })

  ###### F_02_01: Solve, Live and Ghost ########################################
  # Note: The parameter list, the solution, the diagnostics, the simulation
  #   and the impulse responses, each at the sliders and at the ghost.

  F_02_01_par_r <- reactive({
    req(!is.null(input$alpha))
    F_01_01_assemble_fn(list(alpha = val("alpha"), beta = val("beta"),
                             delta = val("delta"), eta = val("eta"),
                             rho = val("rho")))
  })

  F_02_02_sol_r <- reactive(C_01_04_solve_fn(F_02_01_par_r()))
  F_02_03_dia_r <- reactive(C_04_03_diagnostics_fn(F_02_02_sol_r(),
                                                   F_02_01_par_r()))

  F_02_04_ghost_par_r <- reactive({
    ref_lst <- F_01_02_ref_rv()
    if (is.null(ref_lst)) return(B_01_01_default_lst)
    F_01_01_assemble_fn(ref_lst)
  })

  F_02_05_ghost_on_r <- reactive({
    !D_01_02_ghost_off_fn(F_02_01_par_r(), F_02_04_ghost_par_r())
  })

  F_02_06_sim_r <- reactive({
    sol_lst <- F_02_02_sol_r()
    if (!isTRUE(sol_lst$converged_lgl)) return(NULL)
    C_02_01_simulate_fn(sol_lst, F_02_01_par_r(), n_periods_int = 400L)
  })

  F_02_07_ghost_sim_r <- reactive({
    if (!F_02_05_ghost_on_r()) return(NULL)
    g_par <- F_02_04_ghost_par_r()
    g_sol <- C_01_04_solve_fn(g_par)
    if (!isTRUE(g_sol$converged_lgl)) return(NULL)
    C_02_01_simulate_fn(g_sol, g_par, n_periods_int = 400L)
  })

  F_02_08_irf_r <- reactive({
    sol_lst <- F_02_02_sol_r()
    if (!isTRUE(sol_lst$converged_lgl)) return(NULL)
    C_02_02_irf_fn(sol_lst, n_horizon_int = B_01_05_irf_horizon_int)
  })

  F_02_09_ghost_irf_r <- reactive({
    if (!F_02_05_ghost_on_r()) return(NULL)
    g_sol <- C_01_04_solve_fn(F_02_04_ghost_par_r())
    if (!isTRUE(g_sol$converged_lgl)) return(NULL)
    C_02_02_irf_fn(g_sol, n_horizon_int = B_01_05_irf_horizon_int)
  })

  ###### F_02_10: Stage 2.4, the Estimation Chain ##############################
  # Note: Evaluated only when a stage 2.4b output asks. The chain is one
  #   reactive read by both figures, so a redraw runs the sampler once.

  F_02_10_prior_r <- reactive({
    C_06_01_beta_fn(val("prior_mean"), val("prior_sd"))
  })

  F_02_11_obs_r <- reactive({
    C_06_04_sample_fn(F_02_01_par_r(),
                      n_obs_int   = as.integer(val("n_obs")),
                      me_frac_num = B_01_04_chain_lst$me_frac_num,
                      sigma_z_num = B_01_04_chain_lst$sigma_z_num,
                      seed_int    = B_01_04_chain_lst$data_seed)
  })

  F_02_12_chain_r <- reactive({
    obs_lst <- F_02_11_obs_r()
    if (is.null(obs_lst)) return(NULL)
    pri_lst <- F_02_10_prior_r()
    C_07_03_chains_fn(F_02_01_par_r(), obs_lst,
                      a_num       = pri_lst$a_num,
                      b_num       = pri_lst$b_num,
                      step_num    = val("step"),
                      start_vec   = B_01_04_chain_lst$start_vec,
                      n_draw_int  = B_01_04_chain_lst$n_draw_int,
                      warm_up_int = B_01_04_chain_lst$warm_up_int,
                      seed_vec    = B_01_04_chain_lst$seed_vec)
  })

  F_02_13_dens_r <- reactive({
    ch_lst <- F_02_12_chain_r()
    if (is.null(ch_lst)) return(NULL)
    pri_lst <- F_02_10_prior_r()
    C_07_04_density_fn(ch_lst$kept_vec, pri_lst$a_num, pri_lst$b_num)
  })

  ###### F_03_01: Worked Example, Prompt, Problems #############################
  # Note: The card title and story, the loaded example's prompt with the
  #   stage's guidance under it, and the problems panel.

  output$preset_title <- renderUI({
    T_05_06_preset_title_fn(F_01_05_scn_r(), input$stage, B_02_01_stage_vec,
                            stage_word = "")
  })

  output$scenario_story <- renderUI({
    T_05_02_story_fn(F_01_05_scn_r(), B_03_07_controls_lst, B_03_08_help_lst)
  })

  output$prompt <- renderUI({
    scn_lst <- F_01_05_scn_r()
    if (!is.null(scn_lst)) scn_lst$prompt <- HTML(scn_lst$prompt)
    guide_lst <- lapply(B_04_04_guidance_lst, HTML)
    tagList(
      if (!is.null(scn_lst) && identical(scn_lst$stage, input$stage)) {
        div(class = "mb-2",
            T_07_12_prompt_fn(scn_lst, input$stage, guide_lst))
      },
      div(class = "mb-2", T_07_12_prompt_fn(NULL, input$stage, guide_lst))
    )
  })

  # parameter values with no stable solution are named above the figures
  output$problems <- renderUI({
    T_07_13_problems_fn(lapply(F_02_03_dia_r()$problems, HTML))
  })

  ###### F_03_02: Tiles ########################################################
  # Note: Live readouts on the stages whose figures move with the sliders.

  output$tiles <- renderUI({
    if (!input$stage %in% c("1.5b", "2.1b", "2.1c", "2.4b")) return(NULL)
    dia_lst <- F_02_03_dia_r()
    rt      <- dia_lst$roots_lst
    if (is.null(rt)) return(NULL)
    ok_lgl  <- rt$outside_n == rt$jump_n && rt$on_circle_n == 0L
    res_num <- max(dia_lst$residual_lst$state_num,
                   dia_lst$residual_lst$shock_num)
    T_04_03_row_fn(
      T_04_01_tile_fn(
        "Roots outside the unit circle",
        sprintf("%d of %d", rt$outside_n, rt$jump_n),
        if (ok_lgl) "One per jump variable: a unique stable solution"
        else "Not one per jump variable: no unique stable solution",
        class = if (ok_lgl) "good" else "bad"),
      T_04_01_tile_fn(
        "Investment share, s<sub>i</sub>", T_02_05_num_fn(dia_lst$s_i, 3),
        "Implied by &alpha;, &beta; and &delta;, not chosen"),
      T_04_01_tile_fn(
        "Return coefficient, &mu;", T_02_05_num_fn(dia_lst$mu, 4),
        "1 &minus; &beta;(1 &minus; &delta;)"),
      T_04_01_tile_fn(
        # scientific, since scipen = 999 would write out eighteen digits
        "Structural residual", formatC(res_num, format = "e", digits = 1),
        "Anything above 10<sup>&minus;6</sup> is a bug",
        class = if (res_num <= 1e-06) "good" else "bad")
    )
  })

  ###### F_04_00: One Function Per Figure ######################################
  # Note: The only place a builder is called with its arguments. F_04_01
  #   renders it and F_04_02 writes it to PNG and PDF, so the download cannot
  #   drift from the panel.

  F_04_00_build_lst <- list(

    # stage 1.1a is the root-path builder with three stable roots, twenty-one
    # periods from the shock itself
    plot_persistence = function() D_02_01_roots_fn(
      lambda_vec       = B_01_03_persistence_vec,
      n_periods_int    = 21L,
      first_period_int = 0L,
      shock_period_int = 0L,
      title_chr        = "Same Impact, Three Rates of Decay",
      x_title          = expression(bold("Periods Since the Shock (" * h *
                                           ")")),
      y_title          = expression(bold("Response (" * y[t + h] * ")"))),

    plot_roots     = function() D_02_01_roots_fn(),
    plot_bubble    = function() D_04_01_bubble_fn(),
    plot_approx    = function() D_04_02_approx_fn(),
    plot_ratios    = function() D_04_03_ratios_fn(F_02_01_par_r(), "delta"),
    plot_branch    = function() D_02_02_branch_fn(F_02_01_par_r()),
    plot_cycles    = function() D_03_01_cycles_fn(F_02_06_sim_r()),
    plot_tech      = function() D_03_02_technology_fn(F_02_06_sim_r(),
                                                      F_02_07_ghost_sim_r()),
    plot_growth    = function() D_03_03_growth_fn(F_02_06_sim_r()),
    plot_irf       = function() D_03_04_irf_fn(F_02_08_irf_r(),
                                               F_02_09_ghost_irf_r()),
    plot_irf_tech  = function() D_03_05_irf_tech_fn(F_02_08_irf_r(),
                                                    F_02_09_ghost_irf_r()),
    plot_priors    = function() D_05_01_priors_fn(),
    plot_trace     = function() D_05_02_trace_fn(F_02_12_chain_r()),
    plot_posterior = function() D_05_03_posterior_fn(F_02_13_dens_r(),
                                                     F_02_10_prior_r(),
                                                     F_02_12_chain_r())
  )

  # a builder returns NULL with no stable solution; the card shows the
  # placeholder and a Save writes it
  F_04_00_plot_fn <- function(id_chr) {
    p <- F_04_00_build_lst[[id_chr]]()
    if (is.null(p)) {
      T_02_02_placeholder_fn("No stable solution at these settings.")
    } else {
      p
    }
  }

  ###### F_04_01: Figures ######################################################
  # Note: Registered in a loop over the same list; local() fixes the id.

  for (id_chr in names(F_04_00_build_lst)) {
    local({
      one_chr <- id_chr
      output[[one_chr]] <- renderPlot(
        D_01_04_screen_fn(F_04_00_plot_fn(one_chr),
                          title_to_caption = identical(one_chr, "plot_trace")),
        res = D_01_04_res_int)
    })
  }

  ###### F_04_02: Save PNG and Save PDF Under Every Figure #####################
  # Note: T_07_07h_exports_fn at deck size, named rbc-<stage>-<figure>; the
  #   two impulse-response panels are written at the pair size.

  for (id_chr in names(F_04_00_build_lst)) {
    local({
      one_chr <- id_chr
      T_07_07h_exports_fn(
        output, one_chr,
        plot_fn = function() F_04_00_plot_fn(one_chr),
        stem    = function() D_01_03a_figstem_fn(one_chr, input$stage),
        pair    = identical(B_03_05_figfile_lst[[one_chr]]$shape_chr, "pair"))
    })
  }

  ###### F_04_03: The Notes Under the Figures ##################################
  # Note: Each note goes in its card's caption slot, after the figure's own
  #   lifted caption, and is computed from the objects the panel draws.

  F_04_03_note_fn <- function(id_chr, ...) {
    store <- T_02_01d_capstore_fn()
    cap_chr <- if (is.null(store)) NULL else store[[id_chr]]
    tagList(
      if (!is.null(cap_chr) && nzchar(cap_chr)) tags$p(cap_chr),
      tags$p(HTML(paste0(...))))
  }

  output$plot_branch__cap <- renderUI({
    par_lst <- F_02_01_par_r()
    res_fn  <- function(branch_chr) {
      s_lst <- C_01_04_solve_fn(par_lst, branch_chr = branch_chr)
      if (!isTRUE(s_lst$converged_lgl)) return(NA_real_)
      r_lst <- C_04_02_residual_fn(s_lst)
      max(r_lst$state_num, r_lst$shock_num)
    }
    F_04_03_note_fn("plot_branch",
      "Structural residual, transversality imposed: ",
      formatC(res_fn("stable"), format = "e", digits = 1),
      " &nbsp;&middot;&nbsp; the other root: ",
      formatC(res_fn("unstable"), format = "e", digits = 1),
      "<br>Log scale, absolute value of output. Both are exact ",
      "solutions of the model's equations. The ",
      "condition that separates them is not algebra.")
  })

  output$plot_persistence__cap <- renderUI({
    path_df <- C_03_01_roots_fn(B_01_03_persistence_vec, n_periods_int = 21L,
                                first_period_int = 0L, shock_period_int = 0L)
    at_fn <- function(n_int) {
      keep_df <- path_df[path_df$mod_period_tm == n_int, , drop = FALSE]
      paste(sprintf("&rho; = %.1f: %.3f", keep_df$mod_root_val,
                    keep_df$mod_path_val),
            collapse = " &nbsp;&middot;&nbsp; ")
    }
    F_04_03_note_fn("plot_persistence",
      "Share of the shock left after 5 periods &nbsp; ", at_fn(5L),
      "<br>after 10 periods &nbsp; ", at_fn(10L),
      "<br>Every path starts at one. Persistence changes how fast the ",
      "shock goes, not how big it is on arrival.")
  })

  output$plot_priors__cap <- renderUI({
    F_04_03_note_fn("plot_priors",
      "Normal(0, 1) on the whole line &nbsp;&middot;&nbsp; ",
      "Gamma(shape 2, rate 4) on the positive half &nbsp;&middot;&nbsp; ",
      "Beta(2, 2) and Beta(5, 2) on the unit interval",
      "<br>Drawn, not estimated. The parameters are illustrative and ",
      "none of them is load-bearing.")
  })

  output$plot_trace__cap <- renderUI({
    ch_lst <- F_02_12_chain_r()
    if (is.null(ch_lst)) return(NULL)
    F_04_03_note_fn("plot_trace",
      "Acceptance rate &nbsp; chain 1 ",
      sprintf("%.2f", ch_lst$accept_vec[1L]),
      " &nbsp;&middot;&nbsp; chain 2 ",
      sprintf("%.2f", ch_lst$accept_vec[2L]),
      " &nbsp;&middot;&nbsp; step size ",
      sprintf("%.3f", ch_lst$step_num),
      "<br>", ch_lst$n_draw_int, " draws per chain, the first ",
      ch_lst$warm_up_int, " discarded. A common target is roughly one ",
      "accepted proposal in four; it is a convention, not a theorem.",
      "<br>The sample is SIMULATED from this model at &rho; = ",
      sprintf("%.2f", ch_lst$truth_num),
      ", which is the dotted rule. It is not data.")
  })

  output$plot_posterior__cap <- renderUI({
    ch_lst  <- F_02_12_chain_r()
    pri_lst <- F_02_10_prior_r()
    if (is.null(ch_lst)) return(NULL)
    F_04_03_note_fn("plot_posterior",
      "Prior &nbsp; Beta(", sprintf("%.2f", pri_lst$a_num), ", ",
      sprintf("%.2f", pri_lst$b_num), "), mean ",
      sprintf("%.3f", pri_lst$mean_num), ", standard deviation ",
      sprintf("%.3f", pri_lst$sd_num),
      "<br>Posterior &nbsp; mean ",
      sprintf("%.3f", ch_lst$post_mean_num),
      ", standard deviation ", sprintf("%.3f", ch_lst$post_sd_num),
      ", from ", length(ch_lst$kept_vec), " kept draws",
      "<br>The SIMULATED sample was drawn at &rho; = ",
      sprintf("%.2f", ch_lst$truth_num), ". A posterior sitting on top ",
      "of the prior means the sample changed nothing.")
  })

  output$plot_approx__cap <- renderUI({
    raw_df <- C_05_01_approx_fn()
    at_fn  <- function(d_num) {
      i <- which.min(abs(raw_df$mod_deviation_val - d_num))
      sprintf("%.2f%%", raw_df$mod_error_pct_val[i])
    }
    F_04_03_note_fn("plot_approx",
      "Approximation error at a deviation of &nbsp;",
      "5%: ", at_fn(0.05), " &nbsp;&middot;&nbsp; ",
      "10%: ", at_fn(0.10), " &nbsp;&middot;&nbsp; ",
      "25%: ", at_fn(0.25), " &nbsp;&middot;&nbsp; ",
      "50%: ", at_fn(0.50),
      "<br>Exact: e<sup>x</sup> &nbsp;&middot;&nbsp; linear: ",
      "1 + x, the first-order Taylor approximation at zero.",
      "<br>Quarterly output rarely moves more than a few per cent from ",
      "trend, which is the defence of the method.")
  })

  # the pair's note; at Whelan's calibration output 2.13 and hours 1.70 on
  # impact, output 0.16 against technology 0.08 at quarter 50 [W7 36, 37]
  output$plot_irf_tech__cap <- renderUI({
    irf_df <- F_02_08_irf_r()
    if (is.null(irf_df)) return(NULL)
    last_int <- nrow(irf_df)
    F_04_03_note_fn("plot_irf_tech",
      "On impact &nbsp; output ", sprintf("%.2f", irf_df$mod_output_val[1]),
      " &nbsp;&middot;&nbsp; consumption ",
      sprintf("%.2f", irf_df$mod_consumption_val[1]),
      " &nbsp;&middot;&nbsp; hours ",
      sprintf("%.2f", irf_df$mod_hours_val[1]),
      " &nbsp;&middot;&nbsp; investment ",
      sprintf("%.2f", irf_df$mod_investment_val[1]),
      "<br>After ", irf_df$mod_horizon_tm[last_int], " quarters &nbsp; ",
      "output ", sprintf("%.2f", irf_df$mod_output_val[last_int]),
      " against technology ",
      sprintf("%.2f", irf_df$mod_technology_val[last_int]),
      ". Investment is on the cycles figure, where its scale does not ",
      "flatten the others.")
  })

  output$plot_cycles__cap <- renderUI({
    mom_lst <- C_02_03_moments_fn(F_02_06_sim_r())
    if (is.null(mom_lst)) return(NULL)
    rel_vec <- mom_lst$sd_rel_vec
    F_04_03_note_fn("plot_cycles",
      "Relative standard deviations, output = 1 &nbsp;&middot;&nbsp; ",
      "consumption ", round(rel_vec[["consumption"]], 2),
      " &nbsp;&middot;&nbsp; investment ",
      round(rel_vec[["investment"]], 2),
      " &nbsp;&middot;&nbsp; hours ", round(rel_vec[["hours"]], 2),
      "<br>Output growth, first-order autocorrelation ",
      round(mom_lst$growth_acf, 3),
      " &nbsp;&middot;&nbsp; correlation of output with technology ",
      round(mom_lst$tech_corr, 3))
  })

  ###### F_05_01: The Tests Tab ################################################
  # Note: B_05_01's four tests with the model side live.

  output$test_ui <- renderUI({
    sim_df <- F_02_06_sim_r(); irf_df <- F_02_08_irf_r()
    if (is.null(sim_df) || is.null(irf_df)) {
      return(tags$p("No stable solution at these settings."))
    }
    mom_lst <- C_02_03_moments_fn(sim_df)
    now_fn <- function(i_int) {
      switch(i_int,
             sprintf("%+.3f", mom_lst$growth_acf),
             sprintf("correlation %.3f", mom_lst$tech_corr),
             sprintf("%+.2f on impact", irf_df$mod_hours_val[1]),
             sprintf("%.1f times output, against %.1f",
                     mom_lst$sd_rel_vec[["investment"]],
                     mom_lst$sd_rel_vec[["consumption"]]))
    }
    rows_lst <- lapply(seq_along(B_05_01_test_lst), function(i_int) {
      one_lst <- B_05_01_test_lst[[i_int]]
      ok_lgl  <- isTRUE(one_lst$pass(mom_lst, irf_df))
      tags$tr(
        tags$td(one_lst$name),
        tags$td(one_lst$data),
        tags$td(now_fn(i_int)),
        tags$td(if (ok_lgl) "matches" else "does not"),
        tags$td(class = "chg-note", one_lst$src))
    })
    tagList(
      div(class = "eq-group-title", "What the Model Is Held Against"),
      tags$table(class = "table table-sm",
        tags$thead(tags$tr(tags$th("Test"), tags$th("The Data"),
                           tags$th("This Model"), tags$th(""),
                           tags$th("Source"))),
        tags$tbody(rows_lst)),
      # no corner of the sliders passes all three; the tests sweep them
      tags$p(class = "stat-caption", HTML(paste0(
        "Three of the four are the criticisms the lecture ends on, and at ",
        "Whelan's calibration the model fails all three. A corner of the ",
        "sliders can flip one of them - at &eta; = 4 and &rho; = 0.99 the ",
        "wealth effect makes hours fall on impact - but no setting passes ",
        "all three at once."))))
  })

  ###### F_05_02: The Diagnostics Tab ##########################################
  # Note: The root count, the residual and the implied parameters, live on
  #   every stage.

  output$diagnostic_ui <- renderUI({
    dia_lst <- F_02_03_dia_r()
    sol_lst <- F_02_02_sol_r()
    if (is.null(dia_lst$roots_lst)) {
      return(tags$p("No stable solution at these parameters."))
    }
    rt <- dia_lst$roots_lst
    tagList(
      div(class = "eq-group-title", "Blanchard and Kahn, Counted"),
      tags$p(HTML(paste0(
        "The reduced system has three roots: ",
        paste(round(rt$modulus_vec, 4), collapse = ", "), ".<br>",
        "Inside the unit circle: <b>", rt$inside_n, "</b> against <b>",
        rt$predetermined_n,
        "</b> predetermined variables (capital and technology).<br>",
        "Outside: <b>", rt$outside_n, "</b> against <b>", rt$jump_n,
        "</b> jump variable (consumption).<br>",
        "On the circle: <b>", rt$on_circle_n, "</b>."))),
      div(class = "eq-group-title", "Does the Solution Satisfy the Model?"),
      tags$p(HTML(paste0(
        "Substituting the solution back into the structural equations ",
        "leaves a residual of ",
        formatC(max(dia_lst$residual_lst$state_num,
                    dia_lst$residual_lst$shock_num), format = "e", digits = 1),
        ". Anything above 10<sup>&minus;6</sup> is a bug."))),
      div(class = "eq-group-title", "Implied Parameters"),
      tags$p(HTML(paste0(
        "Investment share s<sub>i</sub> = ", round(dia_lst$s_i, 4),
        " &nbsp;&middot;&nbsp; return coefficient &mu; = ",
        round(dia_lst$mu, 5), "<br>",
        "Consumption policy: c<sub>t</sub> = ",
        round(sol_lst$policy_k, 4), " k<sub>t&minus;1</sub> + ",
        round(sol_lst$policy_a, 4), " a<sub>t</sub>")))
    )
  })

  ###### F_06_01: The Equations in Force at This Stage #########################
  # Note: An item is in force once its earliest version is at or before the
  #   stage; the latest version at or before it is shown, flagged "new" or
  #   "changed" when it arrives here. Stage codes are ordered by
  #   B_02_01a_rank_fn, which is why this is not T_06_03 to T_06_06.

  output$eq_title <- renderText({
    trimws(names(B_02_01_stage_vec)[match(input$stage, B_02_01_stage_vec)])
  })

  F_06_01_items_r <- reactive({
    now_int <- B_02_01a_rank_fn(input$stage)
    if (is.na(now_int)) return(list())
    keep_lst <- Filter(function(one_lst) {
      min(B_02_01a_rank_fn(names(one_lst$versions))) <= now_int
    }, B_04_01_equation_lst)
    lapply(keep_lst, function(one_lst) {
      rank_vec <- B_02_01a_rank_fn(names(one_lst$versions))
      cur_int  <- max(rank_vec[rank_vec <= now_int])
      cur_chr  <- names(one_lst$versions)[rank_vec == cur_int]
      status_chr <- if (cur_int != now_int) {
        ""
      } else if (cur_int == min(rank_vec)) {
        "new"
      } else {
        "changed"
      }
      was_chr <- if (identical(status_chr, "changed")) {
        one_lst$versions[[
          names(one_lst$versions)[rank_vec == max(rank_vec[rank_vec < cur_int])]
        ]]
      }
      list(group = one_lst$group, label = one_lst$label,
           tex = one_lst$versions[[cur_chr]], status = status_chr,
           was = was_chr, note = one_lst$notes[[cur_chr]])
    })
  })

  ###### F_06_02: The Two Small Renderers ######################################
  # Note: The LaTeX is escaped before it is marked as HTML, since a bare < or
  #   & in an equation opens a tag; the browser decodes them for MathJax.

  F_06_02_mj_fn <- function(tex_chr) {
    safe_chr <- gsub("&", "&amp;", tex_chr, fixed = TRUE)
    safe_chr <- gsub("<", "&lt;", safe_chr, fixed = TRUE)
    safe_chr <- gsub(">", "&gt;", safe_chr, fixed = TRUE)
    HTML(paste0("\\(", safe_chr, "\\)"))
  }

  F_06_03_flag_fn <- function(status_chr) {
    if (nzchar(status_chr)) {
      tags$span(class = paste0("eq-flag eq-", status_chr), status_chr)
    }
  }

  ###### F_06_04: Panel One, the Equations #####################################
  # Note: Four headed blocks; an empty block says so.

  output$equation_ui <- renderUI({
    item_lst <- F_06_01_items_r()
    block_fn <- function(grp_chr) {
      row_lst <- lapply(
        Filter(function(x_lst) identical(x_lst$group, grp_chr), item_lst),
        function(x_lst) {
          tags$tr(tags$td(class = "eq-label", HTML(x_lst$label),
                          F_06_03_flag_fn(x_lst$status)),
                  tags$td(class = "eq-math", F_06_02_mj_fn(x_lst$tex)))
        })
      div(class = "eq-group",
          div(class = "eq-group-title", B_04_02_groups_vec[[grp_chr]]),
          if (length(row_lst) == 0L) {
            div(class = "chg-note text-muted",
                "Nothing here yet at this stage.")
          } else {
            tags$table(class = "eq-table", do.call(tagList, row_lst))
          })
    }
    withMathJax(tagList(
      do.call(layout_columns, c(
        list(col_widths = breakpoints(sm = 12, md = c(6, 6, 6, 6),
                                      xl = c(7, 5, 7, 5))),
        lapply(names(B_04_02_groups_vec), block_fn))),
      div(class = "eq-legend",
          tags$span(class = "eq-flag eq-new", "new"), " and ",
          tags$span(class = "eq-flag eq-changed", "changed"),
          " mark what this stage adds to the one before. The In Words ",
          "tab says what each one does.")
    ))
  })

  ###### F_06_05: Panel Two, the Notation ######################################
  # Note: One column per group, only the symbols reached; a symbol arriving
  #   at this stage is flagged.

  output$notation_ui <- renderUI({
    now_int  <- B_02_01a_rank_fn(input$stage)
    item_lst <- Filter(function(x_lst) B_02_01a_rank_fn(x_lst$from) <= now_int,
                       B_04_03_notation_lst)
    col_fn <- function(grp_vec, title_chr) {
      keep_lst <- Filter(function(x_lst) x_lst$grp %in% grp_vec, item_lst)
      div(div(class = "eq-group-title", title_chr),
          tags$table(class = "nota-table", lapply(keep_lst, function(x_lst) {
            tags$tr(
              tags$td(F_06_02_mj_fn(x_lst$sym)),
              tags$td(HTML(paste0(toupper(substr(x_lst$txt, 1L, 1L)),
                                  substring(x_lst$txt, 2L))),
                      if (identical(x_lst$from, input$stage) && now_int > 1L) {
                        tags$span(class = "eq-flag eq-new", "new")
                      }))
          })))
    }
    withMathJax(layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(4, 4, 4)),
      col_fn("var", "Variables"),
      col_fn("par", "Parameters"),
      col_fn(c("tgt", "shk"), "Targets, Thresholds and Shocks")
    ))
  })

  ###### F_06_06: Panel Three, In Words ########################################
  # Note: Every equation in force with its note, and what a changed item was.

  output$explain_ui <- renderUI({
    item_lst  <- F_06_01_items_r()
    block_lst <- lapply(names(B_04_02_groups_vec), function(grp_chr) {
      keep_lst <- Filter(function(x_lst) identical(x_lst$group, grp_chr),
                         item_lst)
      if (length(keep_lst) == 0L) return(NULL)
      tagList(
        tags$tr(tags$td(colspan = "3", class = "eq-group-title",
                        B_04_02_groups_vec[[grp_chr]])),
        lapply(keep_lst, function(x_lst) {
          tags$tr(
            tags$td(class = "eq-label", HTML(x_lst$label),
                    F_06_03_flag_fn(x_lst$status)),
            tags$td(div(F_06_02_mj_fn(x_lst$tex)),
                    if (!is.null(x_lst$was)) {
                      div(class = "chg-was", "was ", F_06_02_mj_fn(x_lst$was))
                    }),
            tags$td(class = "chg-note", HTML(x_lst$note)))
        })
      )
    })
    withMathJax(tags$table(class = "eq-table eq-explain", block_lst))
  })
}

################################################################################
## G: Run ######################################################################
################################################################################
# Note: Launch.

#### G_01: Launch ##############################################################
# Note: UI from E, server from F.

shinyApp(E_02_04_ui, F_01_01_server)

#--------------------------------- Script End ---------------------------------#
