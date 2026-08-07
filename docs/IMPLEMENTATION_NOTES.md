# Implementation Notes

Date: 2026-08-07

## Cleanup Scope

This `open_source_release` folder is a cleaned public-code view of the original research working directory. The original numerical update formulas were intentionally kept close to the source code in `myLeapfrog/myLeapfrog.c`.

Changes made for open-source readiness:

- Replaced hard-coded absolute paths with command-line options:
  - `--nt-points`
  - `--points`
  - `--output-dir`
- Renamed the solver header from the ambiguous `complex.h` usage to `leapfrog.h`.
- Renamed `restrict` to `restrict_grid` to avoid the C99 keyword conflict.
- Added file-open checks and readable error messages.
- Added output-directory creation.
- Added optional `--no-phi` mode to avoid writing very large `phi*.m` outputs during quick checks.
- Clamped point-neighborhood initialization loops to the allocated grid range to avoid out-of-bounds writes near domain boundaries.
- Made `cube_max`, `error`, and `error3` respect their `nxt`, `nyt`, and `nzt` arguments.
- Removed damaged comments and unused local variables from the public C files.
- Rewrote MATLAB scripts as path-parameterized functions.
- Rewrote the Python preprocessing script as a command-line utility.

No large datasets or generated `phi*.m` result files are included.

## Alignment With The Paper

The official article preview describes the method as:

- point-cloud-driven 3D volume reconstruction;
- a phase-field Allen-Cahn-type model with a non-negative control function;
- finite-difference spatial discretization;
- a linear, second-order, energy-stable Leap-frog time-stepping scheme;
- numerical tests on objects including horse, teapot, and costa.

The public solver matches these implementation-level features:

- `initialization()` reads scattered point-cloud samples and constructs the control-related field `gf`.
- `Heat3d2()` advances the phase field using current and previous states, matching a Leap-frog-type update structure.
- `relax2()` solves a linearized finite-difference update.
- The solver records the discrete energy history in `Oene.m`.
- MATLAB scripts reconstruct the zero level-set surface with `isosurface`.

## Remaining Manual Verification

The ScienceDirect public preview does not expose every equation in the numerical scheme. Before final publication, verify these points against the full paper PDF or final manuscript:

- The Leap-frog startup strategy. The current code initializes `ooc = oc` and calls `Heat3d2()` from the first time step. There is a disabled `Heat3d1()` block in the original source. If the paper uses a distinct first-step method, restore and document it.
- The exact parameter values for each numerical example, especially `dt`, `T`, `gam`, `Sb`, domain extents, and grid dimensions.
- The exact object-specific coordinate scaling used for each dataset.
- Whether `gf = 1 - phi^2` is the final control-function implementation used for all paper figures.

These are not changed here because they are numerical-method decisions, not packaging bugs.

## Recommended Data Policy

Do not commit generated volume outputs to Git:

- `phi*.m`
- `Oene*.m`
- `remarks.m`

For reproducibility, publish one of the following separately:

- a small demo point cloud in `examples/`;
- full datasets in a GitHub Release;
- full datasets on Zenodo with a DOI.
