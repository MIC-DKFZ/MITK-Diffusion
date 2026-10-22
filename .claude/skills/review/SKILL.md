---
name: review
description: Review exactly what the user names, against the criteria the user names - e.g. "only the user manuals of the tractography plugin" reviewed for "easy to read for a researcher, not a developer". Like the regular code review, but the scope and the yardstick come from the instruction instead of being the current diff and correctness. Use whenever the user asks to review, check, proofread or critique a named set of files, a directory, a module, a plugin, a branch or a diff with a stated focus, audience or quality goal.
---

# Guided review

`/review <what to review> [, <how to review it>]` - one free-form instruction. Read it as two
parts and say back, in one line, how you split it:

- **Scope** - which files. "only the user manuals of the tractography plugin", "the
  FiberTracking module", "my uncommitted changes", "branch 12-port-fibertracking".
- **Focus** - the yardstick. "easy to read for a Workbench user", "uses current MITK
  interfaces", "CMake wiring complete in every build configuration", "thread safety".

Missing scope means the working-tree changes in `mitk-diffusion/`. Missing focus means the
regular correctness review. Ask only if a guess would send the review at the wrong files;
otherwise state your reading and continue.

## Flags

Same flags as the regular code review, anywhere in the instruction; strip them before reading
the rest as scope and focus.

- `--fix` - after reporting, apply the findings to the working tree. Only the ones reported,
  nothing else, and never a commit.
- `--comment` - attach the findings as inline comments on the code where the host offers it
  (the editor's review comments), instead of listing them in the reply. GitLab stays
  read-only: never post there, print the text instead.
- `low` / `medium` / `high` / `max` - effort. Low and medium report few findings you are sure
  of; high and max widen the sweep and may include uncertain ones, marked as such. With no
  level given, use the one from the last review in this session, else medium.

## 1. Resolve the scope to a file list

Turn the words into actual paths (`git -C mitk-diffusion status --short`,
`git -C mitk-diffusion diff --name-only master...<branch>`, `ls`, `grep -rl`) and print the
list before reviewing.

Review every file on that list and nothing outside it. If the list looks wrong - empty, huge,
or clearly not what was meant - say so and fix the scope before reading on.

## 2. Review against the focus, not against habit

The focus **replaces** the default criteria; it does not stack on top of them. A readability
review of a user manual does not report missing tests, and an audience review does not report
style nits the audience will never notice.

Read each file whole - never judge from a filename or a diff hunk alone. Judge against what
the focus asks for, in the words of the focus. For a target-audience focus, read as that
audience: what does a researcher using the Workbench already know, what do they have in front
of them, where does the text assume developer context (class names, CMake options) that they
cannot act on?

Two things always stay on: rules from `AGENTS.md` and `mitk/CLAUDE.md` that the scope touches
(comments say *why* and describe the current state, no patient data, no change to `mitk/`), and
a genuine bug or a factually wrong instruction you happen to see. Report those, marked as
outside the requested focus.

## 3. Report

Use `ReportFindings` when it is available, most severe first, one finding per concrete spot
with file and line; otherwise list them the same way in the reply. Each finding says what is
wrong, why it matters *for the stated focus or audience*, and what to write or do instead -
for prose, quote the replacement sentence.

Prefer few, high-signal findings over a sweep of nits. Nothing to report is a valid result;
say the scope is clean rather than inventing findings.

Close with the usual contract: **Answer** (verdict in a sentence or two), **Why**, **Not
checked** (files in scope you could not judge, and what would settle it).

## 4. Do not change anything

Review only - no edits, no commits, no git state. Only an explicit `--fix` in the instruction
lets you apply the findings afterwards, and then only the ones you reported.
