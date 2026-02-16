# AGENTS.md — Build / Edit / Test / Verify Rules

## Tooling defaults
- When you need OpenAI/Codex/API details, ALWAYS consult the `openaiDeveloperDocs` MCP server.
- When you need third-party library/framework docs, consult `context7` before guessing.
- For UI/E2E validation or repro steps, prefer `playwright` MCP.

## Change discipline
- Prefer the smallest correct diff. Avoid refactors unless asked.
- Keep public interfaces stable unless explicitly required.
- Add/adjust tests for every behavior change or bug fix.

## Build/Test/Verify contract (mandatory)
1) Identify the current build/test command(s) from repo scripts/docs.
2) After changes, run:
   - fast checks (lint/typecheck/unit) first
   - then integration/e2e if applicable
3) If any command fails, fix and rerun until green.

## Output format you must follow
- Always provide:
  - What changed (bullet list)
  - Why it changed
  - How it was verified (exact commands + key results)
  - Risks/edge cases

## Resumable workflow (mandatory)

To make work resumable after interruptions (Stop / End session / Reload Window / VS Code restart),
the agent MUST maintain a progress log file:

- File: `CODEx_PROGRESS.md` (repo root)

Rules:
1) At the start of every task, create or update `CODEx_PROGRESS.md`.
2) Before running any build/test command, record the exact command in the log.
3) After each meaningful step, append:
   - what changed
   - commands executed
   - key output (errors or success summary)
   - current status (PASS/FAIL/PARTIAL)
   - next concrete action (1-3 steps)
4) Always keep section "Current State Snapshot (for resume)" up to date.
5) If the session is interrupted, the next session must resume from
   `CODEx_PROGRESS.md` without redoing completed work.

## Dyn3d FSI Debug Lessons (must-follow)

1) Keep `gap` semantics consistent across FSI preprocessing and coupling:
   - If bodies are connected in 2d plate setups, use:
     - `CutOut2d(bd,bgs; gap=0.0)`
     - `T₁ᵀ(...; gap=false, plane=[1,2])`
   - Do not mix `CutOut2d` default gap (`1.0`) with `T₁ᵀ(...; gap=false/true)` arbitrarily.
   - Mismatch can produce rank-deficient body Schur matrix (`Sbmat`) and trigger `SingularException` in coupled marching.

2) For IFHERK-coupled regressions, validate in this order:
   - first run short horizon (`tf = 2*Δt` or `0.02`) to verify each coupled step is executable;
   - then run target horizon (`tf = 0.5` or full configured tf) to confirm no long-step instability.
   This avoids expensive full-horizon retries while root-causing setup mistakes.

3) Notebook automation path safety:
   - notebook cells that infer `repo_root` from `pwd()` may fail under headless `jupyter execute`;
   - when running CI/non-interactive execution, ensure cwd is notebook directory or set absolute `repo_root` explicitly for reproducibility.

