# ViscousFlow Legacy-to-Modern Migration Method

## Context

This repository now runs Whirl-style FSI workflow on `ViscousFlow v0.1.7` (legacy API),
with Julia 1.x compatibility patches in the local legacy checkout.

Legacy APIs used by this repo:

- `ViscousFlow.Systems.NavierStokes`
- `ViscousFlow.TimeMarching.IFHERK_sc2d` (ported compatibility layer)
- `ViscousFlow.TimeMarching.r₁`
- `ViscousFlow.TimeMarching.r₂`
- `ViscousFlow.TimeMarching.plan_constraints`

Modern API target (`ViscousFlow >= 0.6`) is based on:

- `viscousflow_system`
- `init_sol`, `init`, `step!`
- `ConstrainedODEFunction` internals

## Recommended Equivalence Test Cases

To migrate legacy FSI code to modern API with controlled risk, use three levels:

1. Operator-level equivalence
- Same geometry, grid, Reynolds number, and timestep.
- Compare one-step outputs:
  - vorticity norm (`L2`, `Linf`)
  - constraint force norm
- Acceptance: relative error <= `1e-6` for static-body cases.

2. Coupled short-horizon equivalence
- Run 10-50 steps in a fixed benchmark setup.
- Compare:
  - integrated force history (`fx`, `fy`)
  - body generalized coordinates/velocities
  - convergence iterations per coupled step
- Acceptance: relative trend agreement and bounded drift; document tolerances per case.

3. Long-horizon qualitative equivalence
- Full notebook or benchmark horizon.
- Compare:
  - wake topology milestones
  - force-period statistics (mean, RMS, dominant frequency)
- Acceptance: same qualitative regime, no numerical instability/regression.

## Practical Migration Pattern

1. Freeze a legacy baseline:
- Save `thist`, `fx`, `fy`, `solns`, `uhist` snapshots.

2. Build adapter layer for modern API:
- Keep old call-site signatures in repository code.
- Implement adapter functions that route to modern `viscousflow_system + init/step!`.

3. Migrate incrementally:
- Static-body notebooks first.
- Then moving prescribed kinematics.
- Finally fully coupled FSI loops.

4. Record deltas at each phase:
- API changes
- parameter changes
- numerical differences and accepted tolerance

## Pilot Notes From This Repo

- Legacy route (`ViscousFlow v0.1.7`) is a workable bridge to remove direct Whirl dependency.
- Full test suite passes (`Pkg.test()`).
- Notebook suite passes when executed with per-notebook working directory (`18/18`).

Reference log:
- `.logs/notebook_suite_after_legacy_vf_cwd.log`

