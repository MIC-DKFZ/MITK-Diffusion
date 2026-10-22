# Handoff plan template

The sections a handoff plan has, in this order, and why each is there. Leave out a section only
when it would be empty. The plan is read by a fresh Claude session that knows nothing of the
conversation, so every fact it relies on must be in it, with a file and line it can check.

```markdown
# Plan: <issue title>

Handoff plan for a new session in `<repo>/.worktrees/<branch>`. Read `AGENTS.md` first; its
rules apply throughout. Its paths are relative to `<ws>`, which is not this session's working
directory. Checked against `<base>` at `<base commit>`.

## Issue and branch (already created)

- Issue: #<nr> <title> - <link>
- Branch: `<branch>` from `<base>` at `<base commit>`, worktree `<repo>/.worktrees/<branch>`
- Commit message (from the issue skill): `<message>`

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
- Build and tests: <where, as the user decided - build tree, configure command, targets; or
  "no build" and why>.
- The worktree has no `.env` (git-ignored): no `GITLAB_TOKEN`.
  Nothing here may create GitLab content.
- <memory facts the work depends on, copied>

## Work, in order

### 1. <step>

<files with anchors, what changes, why, the edge cases to cover; "factor out, do not copy"
where an existing helper fits>

### 2. ...

## Docs

<each user manual or doc to update, current state only>

## Tests

<new and changed tests per step, then the exact build and ctest commands, run from where; which
need MITK-Data; a skip is reported as a skip>

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
- **Decisions** stop the session from re-asking or second-guessing what the user settled.
- **Facts with anchors** replace the investigation. The session checks an anchor instead of
  searching, and a fact marked unverified is not built on blindly.
- **Environment** carries what lives outside the worktree (MITK, the build tree, the memory,
  the missing `.env`), which a worktree session cannot discover.
- **Work in order** is ordered so each step can be tested before the next builds on it.
- **Open points** keep decisions that belong to the user with the user.
- **Finish** is the explicit permission to commit, scoped to this branch, which AGENTS.md
  otherwise withholds.
