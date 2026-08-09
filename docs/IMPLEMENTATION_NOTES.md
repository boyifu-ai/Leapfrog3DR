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

