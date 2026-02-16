# Whirl -> ViscousFlow Migration Notes

## Confirmed Version Baseline

- `ViscousFlow.jl` (local path `.external/ViscousFlow.jl`) is `0.6.15`.
- `ViscousFlow.jl` declares `julia = "1.10"` in its `Project.toml`.
- Legacy bridge used in this repo migration: `.external/ViscousFlow_legacy` at `v0.1.7`.

## API Mapping Used By This Repo

The FSI template notebook currently uses the old Whirl APIs below:

- `Whirl.Systems.NavierStokes`
- `Whirl.IFHERK_sc2d`
- `Whirl.plan_constraints`
- `Whirl.r₁`
- `Whirl.r₂`
- `Whirl.TimeMarching.RK31`

Their modern ViscousFlow/ImmersedLayers counterparts are:

- `Whirl.Systems.NavierStokes` -> `ViscousFlow.viscousflow_system(...)`
- `Whirl.IFHERK_sc2d` -> `ImmersedLayers.init(u0,tspan,sys; alg=LiskaIFHERK())` + `step!`
- `Whirl.plan_constraints` -> now internalized in `ConstrainedODEFunction(sys)`
- `Whirl.r₁`, `Whirl.r₂` -> now internalized as ODE/boundary RHS assembly in system operators
- `Whirl.TimeMarching.RK31` -> no direct exposed equivalent; modern solver selection is via algorithm objects (`LiskaIFHERK`, `HETrapezoidalAB2`, etc.)

## Equivalence Assessment

- The mapping is **not 1:1 at call-site level**.
- Whirl exposed a stage-level, callback-style IFHERK interface.
- ViscousFlow now exposes a higher-level constrained ODE system API and hides stage internals.
- Therefore, direct textual replacement of Whirl calls in this notebook is not behavior-equivalent.

## Practical Migration Paths

1. High-level solver migration (recommended):
   - Rewrite notebook fluid time marching to `viscousflow_system + init/step!`.
   - Keep Dyn3d coupling at time-step level (not stage level).
   - Re-validate force/motion histories against Whirl baseline.

2. Compatibility adapter layer:
   - Implement an internal adapter that emulates old stage-level API on top of modern solver.
   - More work and higher maintenance cost, but closest to old coupling workflow.

## Notes From FSInteraction.jl Reference

`FSInteraction.jl` is helpful conceptually but not directly reusable with modern `ViscousFlow`:

- `FSInteraction.jl` pins `ViscousFlow = "0.1.7"` in `.external/FSInteraction.jl/Project.toml`.
- Its source imports legacy namespaces:
  - `ViscousFlow.Fields`
  - `ViscousFlow.RigidBodyMotions`
  - `ViscousFlow.Systems`
  - `ViscousFlow.TimeMarching`
- These namespaces do not exist in `ViscousFlow 0.6.15` (Julia 1.10 target).

Implication:

- We cannot copy/paste FSInteraction coupling code into this repository as-is.
- To use FSInteraction code style, we would need either:
  - a legacy ViscousFlow stack (not aligned with Julia 1.10 modernization), or
  - a substantial port of the coupling layer to the new `ImmersedLayers` ODE API.

## Current Status In This Repository

- Route 2 (legacy-compatible migration) has been implemented.
- Whirl dependency in the FSI template notebook was replaced with `ViscousFlow` legacy API usage.
- A compatibility `IFHERK_sc2d` layer was added in local legacy checkout to preserve old coupling call shape.
- Tests pass:
  - `Pkg.test()` -> pass
- Notebook execution pass:
  - `notebook/` + `test/` non-checkpoint notebooks: `18/18` pass
  - log: `.logs/notebook_suite_after_legacy_vf_cwd.log`
