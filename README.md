# Interactive Model: RBC

A Shiny app for teaching the real business cycle (RBC) model: how it is
built, solved and simulated, what it gets right and wrong, and how one of its
parameters can be estimated. Built by [Sam Deegan](https://sam-deegan.com)
for ECON42240 Advanced Macroeconomics, University College Dublin.

**Try it in the browser (nothing to install):**
https://sam-deegan.com/toy-models/rbc/

Current version: **1.0.8** (see [CHANGELOG.md](CHANGELOG.md)). The version
is shown in the app footer; releases are tagged `vX.Y.Z`.

## What it does

The stage selector builds the model up one layer at a time, following the
order of Whelan's notes rather than opening on a solved model:

| Stage | What is added |
|---|---|
| 1 How Long a Shock Lasts | One shock to an AR(1) series at three rates of decay; what persistence does to the impulse response |
| 2 Forward, Backward, Bubbles | The first-order expectational equation, its forward solution, and a rational bubble that solves it too |
| 3 Roots and Stability | A root inside against a root outside the unit circle; the full model solved on the stable branch and on the other root; the Blanchard and Kahn count on the Diagnostics tab |
| 4 Log-Linearising the Model | The one Taylor approximation every equation goes through, and the error it costs at each deviation from steady state |
| 5 The Steady State and Its Ratios | The return 1/β and the ratios it fixes; the investment share against depreciation |
| 6 Simulate, and Check It | Simulated cycles, output against technology, output growth, and Whelan's impulse-response pair; a Tests tab that holds the model against four facts |
| 7 Priors and Their Supports | Normal, gamma and beta priors drawn over ranges wider than their supports |
| 8 Sampling the Posterior | A random-walk Metropolis sampler for ρ, run on a sample **simulated from this model** and observed with error |

Each stage opens on a worked example (three rates of decay, a bubble on top,
a root inside the circle, how good the approximation is, the steady state is
not free, Whelan's calibration, the propagation test, weak propagation, a
weaker wealth effect, three supports, a prior the data can move, a prior that
will not move, steps too large, steps too small). Every slider has a box
beside it for an exact value, and the loaded example's figures stay on screen
as faded ghosts while the sliders move. The Equations, Notation and In Words
tabs show the model as it stands at the chosen stage and flag what that stage
changed; the Tests and Diagnostics tabs report the scorecard and the root
count live. Every figure carries Save PNG at slide size. Periods
are quarters.

**Stage 8 is estimation on simulated data, not on any economy.** The sample
the sampler is given is drawn from this model at the persistence on the
sidebar and then observed with error, so there is a known answer to hold the
posterior against. The app says so on the stage, in the figure subtitles and
in the notes under the figures. Nothing in the app is an estimate of a real
economy.

## Run it locally

1. Install [R](https://cran.r-project.org/) (4.1 or later) and, ideally,
   [RStudio](https://posit.co/download/rstudio-desktop/).
2. Install the three packages once:

   ```r
   install.packages(c("shiny", "bslib", "ggplot2"))
   ```

3. Open `app.R` in RStudio and click **Run App**, or from R in this folder:

   ```r
   shiny::runApp()
   ```

Equations are typeset with MathJax from a CDN, so they need an internet
connection; everything else runs offline. Stage 8 runs two chains of 800
draws on every redraw, which takes about a second.

## Files

```
app.R          the app: settings and text (section B), figures (D),
               interface (E), server (F)
R/model.R      the model: solver, simulator, impulse responses, moments,
               diagnostics, the layers before the simulation, the
               state-space form, the Kalman likelihood and the sampler.
               Sources on its own, so slides can reuse it.
R/toolkit.R    layout and helpers shared with the other toy-model apps
tests/         verify_against_rats.R, the checks described below
README.md      this file
CHANGELOG.md   version history
CONVENTIONS.md how the figures and worked examples are laid out
LICENSE        CC BY-NC-ND 4.0
```

All text on screen (worked examples, prompts, equations, notation) is in
sections `B_02` and `B_04` of `app.R`, so it can be edited without touching
the rest.

## Checks

From the repo root:

```
Rscript tests/verify_against_rats.R
```

The script needs `shiny`, `bslib`, `ggplot2` and `htmltools`, prints a pass
or fail line per check and exits with a non-zero status if any fails. It
takes about three minutes, most of it the sampler.

Sections `V_02` to `V_05` compare the solver against Whelan's RATS program
`rbc.prg`, whose arithmetic is transcribed into the script: the derived
parameters, the shock loading computed his way, the structural residual, and
his brute-force iteration run in R and held against the eigenvalue solution.
`V_06` and `V_07` check the findings Whelan states in part 7 (volatile
investment, output tracking technology, no positive growth autocorrelation,
hours rising against Gali, no propagation from iid technology) and his
impulse responses as read off his figures, then every equation the Equations
panel prints against the solver's own path at four calibrations. `V_08` is
not against RATS, which does not estimate: it checks the scalar Kalman filter
against a plain matrix one, the sampler's acceptance rates, recovery of a
known ρ from a simulated sample, a tight prior returning itself, and the
drawn prior densities against their exact masses. `V_09` and `V_10` check the
app around the model: the export shapes, the equation and notation lists
closing over each other, the page frame, every stage's panels through
`shiny::testServer`, and every Save PNG handler.

## The model

The real business cycle model is the original dynamic stochastic general
equilibrium model: a representative household chooses consumption, hours and
investment to maximise expected utility, markets clear, and the only shock is
to technology. It is the model of part 7 of Whelan's *MA Advanced
Macroeconomics* notes, and part 10 of the same notes rewrites its solution in
state-space form to show how such models are estimated. Every variable is a
log deviation from steady state, in per cent; `k_{t-1}` is the capital
carried into period `t`.

```
Utility:     U(C_t) - V(N_t) = C_t^{1-η}/(1-η) - a N_t
Resource:    y_t = (1 - s_i) c_t + s_i i_t
Capital:     k_t = (1 - δ) k_{t-1} + δ i_t
Production:  y_t = a_t + α k_{t-1} + (1 - α) n_t
Labour:      (1 - α) Y_t/N_t = a C_t^η   ⟹   n_t = y_t - η c_t
Euler:       c_t = E_t c_{t+1} - (1/η) E_t r_{t+1}
Return:      r_t = μ (y_t - k_{t-1}),   μ = 1 - β(1 - δ)
Technology:  a_t = ρ a_{t-1} + ε_t,     σ_ε = √(1 - ρ²)
Steady state: R* = 1/β,  I*/K* = δ,  s_i = I*/Y* = αδ / (β⁻¹ + δ - 1)
```

**Utility** is CRRA in consumption with a linear disutility of work. `η` is
the curvature of utility; its inverse is the intertemporal elasticity of
substitution. `a` is the disutility of an hour.

**The resource constraint** splits output between consumption and investment
in the shares the steady state implies. **Capital accumulation** carries
`1 - δ` of last period's stock forward and adds investment. **Production** is
Cobb-Douglas in logs, so the shares are the coefficients and technology enters
additively; `α` is capital's share.

**The labour condition** sets the wage equal to the marginal rate of
substitution. With `V(N)` linear the right side has no hours in it, so given
consumption the wage is fixed: labour supply is flat at every `η`, and what
`η` sets is the wealth effect, how far higher consumption holds hours back.

**The Euler equation** is the only forward-looking equation, which makes
consumption the model's one jump variable. **The return on capital** rises
with output per unit of capital; its coefficient `μ` and the investment share
`s_i` are not chosen but implied by `α`, `β` and `δ` through the steady-state
requirement `R* = 1/β`.

**Technology** is an AR(1) with persistence `ρ`. The innovation standard
deviation is `√(1 - ρ²)`, Whelan's own convention, so the variance of
technology is one at every `ρ`: moving the slider changes how long a shock
lasts, not how big it is.

**How it is solved.** The four static equations pin output, hours, the return
and investment given `k_{t-1}`, `a_t` and `c_t`, so the model reduces to a
three-equation core in two predetermined variables (capital, technology) and
one jump (consumption). The app writes the core as
`M0 E_t v_{t+1} = M1 v_t`, takes the eigenvalues of `M0⁻¹ M1`, and imposes
the Blanchard and Kahn condition: with one jump variable exactly one root
must lie outside the unit circle, and zeroing that root's coordinate is the
transversality condition that pins consumption. Stage 3 can take the other
root instead, which solves the same equations and walks away from the steady
state. The seven-variable transition and shock loading are rebuilt from the
core, and the solution is substituted back into the structural system so the
residual can be shown on the Diagnostics tab.

**How it is estimated (stage 8).** The solved model is written in state-space
form with state `(k_{t-1}, z_t)`, capital and consumption observed with
error (two measurement errors and one shock give three random terms for two
observables), and the Kalman filter returns the log-likelihood. A beta prior
on `ρ` is combined with it, and a random-walk Metropolis sampler with a
normal proposal of step `c` draws from the posterior; two chains are started
at opposite ends of the unit interval and the warm-up is discarded. Only `ρ`
is estimated; the other parameters and the measurement-error variances are
held at the values the sample was simulated at.

**What the stages show with it**

- *1* Persistence changes how fast a shock goes, not how big it is on
  arrival: every path starts at one and what survives after `h` periods is
  `ρ^h`.
- *2* The forward solution is the discounted sum of expected future `x`, and
  it is not the only solution: any bubble growing at `1/a` satisfies the
  equation too. Nothing in the equation prefers one.
- *3* A root inside the unit circle decays, a root outside explodes. The full
  model has two predetermined variables and one jump, so one root outside is
  what a unique stable solution needs; both branches solve the model exactly
  and only transversality separates them.
- *4* The approximation `e^x ≈ 1 + x` costs a tenth of a per cent at a
  five per cent deviation and about nine per cent at fifty.
- *5* `R* = 1/β` fixes the ratios: moving `δ` moves the investment share
  whether or not that was the intention. A calibration is a set of joint
  choices.
- *6* At Whelan's calibration the model gets one fact right (investment far
  more volatile than consumption) and three wrong: output growth is not
  positively autocorrelated (Cogley and Nason's 0.34 in the data), output
  simply follows technology so the model adds almost no propagation of its
  own, and hours rise after a positive technology shock where Gali finds they
  fall. The Tests tab reports the model side live; no corner of the sliders
  passes all three.
- *7* A prior has to live where the parameter lives: normal for a free sign,
  gamma for a standard deviation, beta for a persistence.
- *8* On a simulated sample, a wide prior centred away from the truth is
  pulled onto it and narrows; a very tight prior returns itself and nothing
  looks wrong; a step too large rejects almost everything, a step too small
  accepts almost everything and crawls, and a high acceptance rate is not
  evidence of convergence.

**Where it departs from the textbook.** Whelan's RATS program solves the
model by iterating on `C ← (I - BC)⁻¹ A` from `C = I`. That map solves a
matrix quadratic with more than one root and converges to whichever root its
arithmetic path leads to; run in R it reaches the explosive fixed point, so
the app selects the stable manifold by eigenvalue instead (the test script
runs both and records the difference). The shock loading is the same object
as the program's. The estimation stage goes beyond part 10, which stops at
the choice of prior: the Metropolis sampler, its step-size slider and the
two-chain trace follow the ECON42240 lecture deck, and the sampler is not
tuned during warm-up because the student is meant to tune it. The sample is
simulated, so the exercise has a known answer that no real estimation has.
The model has one shock, no money, no frictions and no policy; the full DSGE
of Smets and Wouters (Whelan's part 11) and the financial accelerator (part
12) are not in this app, and the missing frictions are what the stage 6
scorecard is about.

## References

- Whelan, K. *MA Advanced Macroeconomics*. Part 6 (solving models with
  rational expectations), part 7 (the real business cycle model), part 10
  (estimating DSGE models). https://www.karlwhelan.com/ma-advanced-macroeconomics/
- Whelan's RATS program `rbc.prg`, for the calibration and the innovation
  standard deviation.
- Blanchard, O. and Kahn, C. (1980). The solution of linear difference models
  under rational expectations. *Econometrica* 48(5).
- Cogley, T. and Nason, J. (1995). Output dynamics in real-business-cycle
  models. *American Economic Review* 85(3).
- Gali, J. (1999). Technology, employment, and the business cycle: do
  technology shocks explain aggregate fluctuations? *American Economic
  Review* 89(1).

## Licence

© Sam Deegan. Released under
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/):
free to use and share for teaching with attribution; not for commercial use
or redistribution in modified form.
