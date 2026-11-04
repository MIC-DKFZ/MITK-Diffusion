# Handoff plan template

The sections a handoff plan has, in this order, and why each is there. Leave out a section only
when it would be empty; "Tests" is never left out. The plan is read by a fresh Claude session that knows nothing of the
conversation, so every fact it relies on must be in it, with a file and line it can check.

```markdown
# Plan: <issue title>

Handoff plan, carried out by a session started in the workspace root `<ws>` (started with
"carry out this plan: <this file>"). All work happens in the worktree `<wt>` =
`<repo>/.worktrees/<branch>`; see "Working in the worktree". `AGENTS.md` is loaded and its
rules apply throughout. Checked against `<base>` at `<base commit>`.

## Issue and branch (already created)

- Issue: #<nr> <title> - <link>
- Branch: `<branch>` from `<base>` at `<base commit>`, worktree `<wt>`
- Commit message (from the issue skill): `<message>`

## Working in the worktree

- Read, search and edit only under `<wt>`, always by absolute path. A `mitk-diffusion/...` path
  in `AGENTS.md`, a skill or this plan means `<wt>/...` for this task. The main checkout
  `<repo>` is not edited.
- Git only as `git -C <wt> ...`. At start and end, `git -C <wt> status --short --branch` shows
  `<branch>`, and `git -C <repo> status --short --branch` is the same as at the start.
- Skills that default to `mitk-diffusion/` (e.g. `review`) get paths under `<wt>` as scope.
- Builds and tests use `<tree>`, the build tree bound to `<wt>` (see "Environment").
- Other sessions may work in other worktrees at the same time; their worktrees and build trees
  are not touched.

## Decisions already made by the user

- <each decision, in the user's words - the implementer does not reopen these>

## Goal

1. <observable outcome>
2. ...

## Facts the design rests on

- <fact> (`path/to/file:line`)
- <fact that was not verified> - marked "not verified", with what would verify it

## Environment

- MITK: `<ws>/mitk`, read-only.
- Build and tests: <where, as the user decided - build tree `<tree>`, configure command,
  targets; or "no build" and why>.
- GitLab: `.env` is readable from `<ws>`, but nothing here creates GitLab content; issue text
  goes into a file for the user.
- <memory facts the work depends on, copied>

## Work, in order

### 1. <step>

<files with anchors, what changes, why, the edge cases to cover; "factor out, do not copy"
where an existing helper fits; the tests that check this step>

### 2. ...

## Docs

<each user manual or doc to update, current state only>

## Tests

- Module tests: <new and changed test cases per step, file and name; the existing tests that
  cover the changed code and must still pass>
- Command line apps: <for a changed app, its case in
  `<wt>/Modules/DiffusionCmdApps/Testing/mitkDiffusionCmdApps<Group>Test.cpp`, or "no app changes">
- UI tests: <for a change a view can show (its code, `.ui` file, or what it calls), the case
  added to or changed in the plugin's suite in
  `<wt>/Plugins/org.mitk.gui.qt.diffusionimaging.uitests/src/internal/`, with the widgets
  (`.ui` object names) and results it checks; what stays a manual check at a display
  (rendering, #45); or "no view affected" and why>
- Commands: <the exact build and ctest commands, run from where, with a `TMPDIR` per tree,
  e.g. `ctest -R mitkDiffusionViewsUiTest_<suite>`; which need the test data
  (`MITK_DIFFUSION_DATA_DIR`); a skip (exit code 77: no Qt platform plugin, no OpenGL) is
  reported as a skip>

## Open points for the user

- <what the implementer must not decide alone: checks on real data, values only the user has -
  state what is done in code meanwhile>

## Finish

- When implemented and tested, commit on `<branch>`. The user asked for this through the
  handoff.
  - Stage only this task's files by name, never `git add -A`.
  - Message: the issue skill's commit message as subject, in the shape of `mitk/CLAUDE.md`
    (blank line, body with what and why). No `Co-Authored-By` or other attribution lines.
- No push, merge request, merge, rebase or tag unless the user asks.
- Report in the user's format (Answer / Files changed / Potential side effects / Not checked /
  Open questions), with each build and test's actual result.
```

Why these sections:

- **Issue and branch** says what already exists, so the session neither creates a second issue
  nor branches again.
- **Working in the worktree** exists because the session starts in the workspace root, so the
  user needs no second window; without it, edits and git commands land in the main checkout.
- **Decisions** stop the session from re-asking or second-guessing what the user settled.
- **Facts with anchors** replace the investigation. The session checks an anchor instead of
  searching, and a fact marked unverified is not built on blindly.
- **Environment** carries what lives outside the worktree (MITK, the build tree, the memory
  facts), so the session does not have to rediscover it.
- **Work in order** is ordered so each step can be tested before the next builds on it.
- **Tests** are part of every handoff, so the implementer writes them instead of only running
  the existing ones; the UI tests are the only automatic check of a view, since nobody uses
  the Workbench at a display during the run.
- **Open points** keep decisions that belong to the user with the user.
- **Finish** is the explicit permission to commit, scoped to this branch, which AGENTS.md
  otherwise withholds.
