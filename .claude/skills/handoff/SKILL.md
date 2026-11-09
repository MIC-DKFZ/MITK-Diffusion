---
name: handoff
description: Investigate a feature or fix read-only, settle the open design decisions with the user, then create the GitLab issue (through the feature-request, bug-report or improvement-request skill), a <nr>-branch in a free worktree slot (mitk-diffusion/.worktrees/slot-<n>, kept built in wt-build/slot-<n>), a self-contained handoff plan, and the start sentence ("carry out this plan: <path>") for a new Claude session in the workspace root that implements, tests and commits in that worktree. Use whenever the user wants work prepared for a separate session or worktree - "hand this off", "prepare a handoff", "plan this for a new session", "create issue, branch and plan", "set this up so another session can implement it" - also when an /investigate discussion should turn into issue + branch + plan, even if the word "handoff" is not used. Not for implementing in the current session, and not for issue text alone (that is feature-request / bug-report / improvement-request).
---

# Handoff: investigate, plan, create issue and branch, hand over

The argument is the feature or fix. With no argument, the source is the current conversation:
what the user asked, what was found and what they decided. Do not redo an investigation the
conversation already holds; re-check only the anchors the plan will cite.

Paths are relative to the workspace root (see `AGENTS.md`); `<ws>` below is its absolute path
and `<repo>` is `<ws>/mitk-diffusion`.

The result is a plan that a new session carries out cold from the **workspace root**: the user
opens a new chat there (a new Claude Code tab in the same VS Code window) and says
`carry out this plan: <ws>/gitlab_issues/<branch>.plan.md`. Nobody switches to the worktree or
opens a second window. That session loads `AGENTS.md`, the skills, the project memory and
`.env` like this one, but has none of this conversation, and its working directory is `<ws>`,
not the worktree. Consequences:

- The issue and the branch are created here, so the implementer neither files nor branches.
- The plan must carry every fact, decision and command the implementer needs.
- The plan pins every edit, git command, build and test to the worktree by absolute path;
  without that, work lands in the main checkout `<repo>`.
- The work runs in a worktree slot (`AGENTS.md`), whose Debug tree `<ws>/wt-build/slot-<n>` is
  already built and recompiles only what the branch changes; `debug/` and `release/` are bound to
  the main checkout.

## 1. Investigate - read-only

Follow the `investigate` rules until step 4: no edits to source, no git state changes, no
GitLab writes. Read the actual code and name file and line for every claim. Look for an
existing issue (`python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py search <key words>`)
and an existing `<nr>-…` branch (`git -C mitk-diffusion branch -a | grep <nr or key word>`).

Find the cause for a fix and the smallest design that fits the existing architecture for a
feature. Note what the design rests on, with anchors, as you go: that list becomes the plan's
"Facts" section and saves the next session from repeating the work.

Find the tests that cover the changed code and where new ones go: the module's `Testing/`
folder; for a command line app its group's smoke test in `Modules/DiffusionCmdApps/Testing/`;
for a view, its `.ui` file or anything a view calls, the plugin's UI test suite in
`Plugins/org.mitk.gui.qt.diffusionimaging.uitests/src/internal/` (`mitkDiffusionViewsUiTest_<suite>`).
Every plan gets tests; the UI tests whenever a view's behaviour can change.

## 2. Settle the decisions

Ask the user only what changes the work and cannot be read from the code: scope, trade-offs,
behaviour choices, and **where the implementing session builds and tests** (by default the
slot's Debug tree `<ws>/wt-build/slot-<n>`; a Release check or another build configuration needs a
tree of its own, see `AGENTS.md`). Ask one question at a time, with a recommendation. Record each answer
in the user's words; the plan states them as decisions, so the implementer does not reopen them.

Some points cannot be settled from here: checks on real data, values only the user has. Do not
block on them. They go into the plan as open points the implementer must not decide alone.

## 3. Confirm before creating anything

Show in a few lines the issue type and title, the branch name (with `<nr>`), the base branch, the
slot,
and the plan's goals, steps and tests (module, command line app and UI tests). Then ask once whether to create the issue and the branch in the
slot. Skip this only if the user's request already said to go ahead without asking.
Creating a GitLab issue is visible to others and is not undone by deleting a file.

## 4. Issue

Invoke the matching skill with `chat create` as its argument, so it writes the issue from this
conversation and creates it in the mitk-diffusion project:

| The work | Skill |
| --- | --- |
| new or changed behaviour | `feature-request` |
| something broken | `bug-report` |
| refactoring, porting, cleanup, docs, tests, build | `improvement-request` |

Keep what the skill returns: issue number, link, title, commit message, and the branch name
with the number filled in.

- **Existing issue** (the skill's search finds one): use its number and do not create another.
- **Creation fails** (no token, network): show the saved issue file's path, ask the user to post
  it and give the number, and wait.
- **Work only in MITK:** it does not belong here. Say so and stop.

## 5. Branch in a worktree slot

The main checkout stays as it is: no switching, no stashing. The base is `master` unless the
user named another. Record the base commit (`git -C <repo> rev-parse --short <base>`) for the
plan.

Take a free slot `<repo>/.worktrees/slot-<n>` (`AGENTS.md`): its HEAD is detached and its working
tree is clean. A slot with a branch checked out belongs to a handoff that is running, not yet
carried out or not yet merged; never take it, even when the branch shows as merged (a new branch
without commits does).

```bash
git -C <repo> worktree list                                    # slots: detached or on a branch
git -C <repo>/.worktrees/slot-<n> status --short               # empty when clean
git -C <repo>/.worktrees/slot-<n> switch -c <branch> <base>    # new branch
git -C <repo>/.worktrees/slot-<n> switch <branch>              # branch already exists
```

When no slot is free, ask the user whether to add one (one full build, about 20 minutes at
`-j8`, 14 GB), then:

```bash
git -C <repo> worktree add --detach <repo>/.worktrees/slot-<n> <base>
<repo>/.claude/scripts/configure-worktree-build.sh <repo>/.worktrees/slot-<n> <ws>/wt-build/slot-<n>
cmake --build <ws>/wt-build/slot-<n> -j8    # in the background
```

Use absolute paths: `git -C <repo>` resolves a relative worktree path inside the repository.
Afterwards, check that the main checkout's status is unchanged.

## 6. Plan file

Write the plan to `gitlab_issues/<branch>.plan.md` in the workspace root, beside the issue file
the skill saved. `gitlab_issues/` lies outside both repositories and survives a reboot. Follow
`references/plan-template.md`; it holds the sections and why each exists.

The plan is written for a fresh Claude session, not for the user: complete, anchored and terse.
Copy into it every memory fact the work depends on; memory changes, the plan must not.

## 7. Hand over

Reply with:

- issue number and link;
- branch, base commit and slot path;
- the plan file's absolute path;
- the start sentence, fenced so it copies in one click, with the absolute plan path:

```text
carry out this plan: <ws>/gitlab_issues/<branch>.plan.md
```

Add one line: paste it into a new chat started in the workspace root (a new Claude Code tab in
the same VS Code window); no worktree switch, no second window. Then list the plan's open
points for the user, if any. Nothing is committed by this skill; the
commit is the implementing session's job, as the plan says.
