# Changelog

All notable changes to this app. Versions follow [Semantic Versioning](https://semver.org/):
MAJOR for a change to the model or its notation, MINOR for new features
(a stage, a worked example, a figure), PATCH for fixes and wording.
Each release is tagged in git as `vX.Y.Z` and shown in the app footer.

## [1.0.8] - 2026-09-28

### App
- Figure titles are back, on the page: each figure's title (with its
  live numbers) is the card header and its subtitle sits above the image,
  so the PNG itself stays bare for slides with their own captions.
- The stage name and the Equations tabs sit on one line with no rule
  under them. The card opens folded, with the stage name drawn as the
  selected tab, so the worked examples and figures sit high on the page;
  a tab opens it and the stage name folds it again.
- Every figure box is 3:2 at any width, and the theme holds the image at
  3:2 inside it.
- Text sizes and the figure card follow the macro apps: Save PNG sits in
  the card header, exports are 3:2 (1500 x 1000, or 1440 x 960 for a pair
  panel) and the tests check that.
- The title is "Real Business Cycle (RBC) Model".

## [1.0.7] - 2026-09-28

### App
- The notes under the figures are rewritten as short plain prose: no bold
  lead-in sentences, one point per paragraph.

## [1.0.6] - 2026-09-28

### App
- Card headers in the blue used for headings, not body grey.
- Save PNG only under each figure; the PDF button is gone.

## [1.0.5] - 2026-09-28

### App
- Cards have no border or header rule: figures, equations and stories sit
  on the page separated by whitespace alone.

## [1.0.4] - 2026-09-28

### App
- Cards, panels, tiles and buttons are square with no shadow: they organise
  the page rather than decorate it.

## [1.0.3] - 2026-09-28

### App
- The QR code returns to the foot of the sidebar, with the name and site
  address, alongside the small one in the title bar.

## [1.0.2] - 2026-09-28

### App
- No figure carries a title or subtitle inside the image; the card header
  and the caption under it name and explain the figure (CONVENTIONS.md 6).
- Figures are drawn on a white ground, so the image sits flat in its card
  instead of showing as a tinted tile.

## [1.0.1] - 2026-09-28

### App
- The In Words tab lays out its three columns at fixed widths, so an
  equation no longer collapses to one term per line beside its note.
- The preset card no longer doubles the word "Stage" in front of a stage
  name that already carries it.

## [1.0.0] - 2026-09-28

First public release as a standalone repository.

### Model
- The seven-equation log-linear RBC model of Whelan's MA Advanced
  Macroeconomics part 7, at his calibration, with the innovation standard
  deviation sqrt(1 - rho^2) of his RATS program.
- Solved by reducing to a three-equation core and selecting the stable
  manifold by eigenvalue (Blanchard and Kahn), with the other root available
  on request; the solution is substituted back into the structural system
  and the residual reported.
- Simulation, impulse responses and moments; scalar root paths for the
  stage 1 and 3 figures; the log-linear approximation, the steady-state
  ratios and the forward solution with a bubble.
- The solved model in state-space form (part 10), a sample SIMULATED from it
  and observed with error, a scalar Kalman log-likelihood, and a random-walk
  Metropolis-Hastings sampler for rho with two chains and a warm-up.

### App
- Eight stages that build the model up one layer at a time, from a single
  AR(1) to estimation on a simulated sample.
- Fourteen worked examples, each belonging to one stage and shown only
  there.
- Equations, Notation and In Words tabs that track the model at each stage;
  Tests and Diagnostics tabs that hold the model against four facts and show
  the root count, the residual and the implied parameters live.
- Fourteen figures with Save PNG and Save PDF at slide size: root paths, the
  two branches, the bubble, the approximation, the investment share,
  simulated cycles, output against technology, growth, the impulse-response
  pair, the prior densities, the two chains and prior against posterior.
- Ghost curves showing the loaded worked example alongside the live sliders.
- Readout tiles for the root count, the investment share, the return
  coefficient and the structural residual.
- A test script, `tests/verify_against_rats.R`, that checks the solver
  against Whelan's RATS program, the findings of part 7, the estimation
  stage against its own definitions, and the app's panels and downloads.
