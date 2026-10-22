---
name: handoff
description: Investigate a feature or fix read-only, settle the open design decisions with the user, then create the GitLab issue (through the feature-request, bug-report or improvement-request skill), a <nr>-branch in its own mitk-diffusion/.worktrees/ worktree, a self-contained handoff plan, and the one command that starts a new Claude session in that worktree to implement, test and commit on the branch. Use whenever the user wants work prepared for a separate session or worktree - "hand this off", "prepare a handoff", "plan this for a new session", "create issue, branch and plan", "set this up so another session can implement it" - also when an /investigate discussion should turn into issue + branch + plan, even if the word "handoff" is not used. Not for implementing in the current session, and not for issue text alone (that is feature-request / bug-report / improvement-request).
---

# Handoff: investigate, plan, create issue and branch, hand over

The argument is the feature or fix. With no argument, the source is the current conversation:
what the user asked, what was found and what they decided. Do not redo an investigation the
conversation already holds; re-check only the anchors the plan will cite.

Paths are relative to the workspace root (see `AGENTS.md`); `<ws>` below is its absolute path
and `<repo>` is `<ws>/mitk-diffusion`.

The result is a new session that can start cold in its own worktree. A session started in a
worktree has none of this conversation, may not load this session's memory (memory is kept per
project path), and has the versioned `.claude/` skills but no `.env`, because it is
git-ignored. That has three consequences:

- The issue and the branch are created here, where `GITLAB_TOKEN` exists.
- The plan must carry every fact, decision and command the implementer needs.
- A worktree is not built by `debug/` or `release/`, which are bound to the main checkout.

## 1. Investigate - read-only

Follow the `investigate` rules until step 4: no edits to source, no git state changes, no
GitLab writes. Read the actual code and name file and line for every claim. Look for an
existing issue (`python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py search <key words>`)
and an existing `<nr>-…` branch (`git -C mitk-diffusion branch -a | grep <nr or key word>`).

Find the cause for a fix and the smallest design that fits the existing architecture for a
feature. Note what the design rests on, with anchors, as you go: that list becomes the plan's
"Facts" section and saves the next session from repeating the work.

## 2. Settle the decisions

Ask the user only what changes the work and cannot be read from the code: scope, trade-offs,
behaviour choices, and **where the implementing session builds and tests** (the worktree has no
build tree of its own). Ask one question at a time, with a recommendation. Record each answer
in the user's words; the plan states them as decisions, so the implementer does not reopen them.

Some points cannot be settled from here: checks on real data, values only the user has. Do not
block on them. They go into the plan as open points the implementer must not decide alone.

## 3. Confirm before creating anything

Show in a few lines the issue type and title, the branch name (with `<nr>`), the base branch,
and the plan's goals and steps. Then ask once whether to create the issue, the branch and the
worktree. Skip this only if the user's request already said to go ahead without asking.
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

## 5. Branch and worktree

The main checkout stays as it is: no switching, no stashing. The base is `master` unless the
user named another. Record the base commit (`git -C mitk-diffusion rev-parse --short <base>`)
for the plan.

```bash
git -C mitk-diffusion worktree add .worktrees/<branch> -b <branch> master   # new branch
git -C mitk-diffusion worktree add .worktrees/<branch> <branch>             # branch already exists
```

Afterwards, check that the main checkout's status is unchanged.

## 6. Plan file

Write the plan to `gitlab_issues/<branch>.plan.md` in the workspace root, beside the issue file
the skill saved. `gitlab_issues/` lies outside both repositories, survives a reboot, and the
worktree session can read it by absolute path. Follow `references/plan-template.md`; it holds
the sections and why each exists.

The plan is written for a fresh Claude session, not for the user: complete, anchored and terse.
Copy into it every memory fact the work depends on, because the worktree session may not see
this session's memory.

## 7. Hand over

Reply with:

- issue number and link;
- branch, base commit and worktree path;
- the plan file's absolute path;
- the start command, fenced so it copies in one click, with absolute paths; `--add-dir` gives
  the session access to `mitk/` and the build trees beside the repo:

```bash
cd <repo>/.worktrees/<branch> && claude --add-dir <ws> "Read <ws>/gitlab_issues/<branch>.plan.md and carry it out."
```

Add one line for VS Code: open the worktree folder in a new window and paste the quoted prompt.
Then list the plan's open points for the user, if any. Nothing is committed by this skill; the
commit is the implementing session's job, as the plan says.
