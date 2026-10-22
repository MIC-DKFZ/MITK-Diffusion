---
name: improvement-request
description: Write a short GitLab issue in plain language for a refactoring, cleanup, porting or other improvement of existing code, docs, tests or build in mitk-diffusion - one main topic under 200 words, smaller findings as one-line "other points" - plus title, markdown body (current state, why it hurts, proposed change, done-when checklist), commit message and a <nr>-branch-name branch; creates the issue in the mitk-diffusion GitLab project only when explicitly asked to. Use when the user asks for a refactoring issue, tech-debt, porting or cleanup ticket, improvement issue, or to turn review findings into a GitLab issue - i.e. nothing is broken for a user (that is bug-report) and no new user-visible behaviour is wanted (that is feature-request).
---

# GitLab improvement issue

Turn a review finding, a refactoring or porting idea or an improvement (or one already made)
into a short, copy-pasteable GitLab issue for `mic/personal/group6/neuro/mitk-diffusion` on
git.dkfz.de. Paths are relative to the workspace root (see `AGENTS.md`).

## Is it an improvement?

An improvement leaves what MITK Diffusion does alone and changes how well the code, docs,
tests or build hold up: duplicated logic, use of an outdated MITK, ITK, VTK or Qt interface,
dead code, a misleading name, a stale manual, a missing test, a slow step, a structure that
makes every change bigger than it should be.

- Something a user sees is **broken** - use the `bug-report` skill.
- Someone wants **new behaviour** they could ask for - use the `feature-request` skill.
- A review found both - write the improvement here and say in the note what belongs in a
  bug report or feature request instead.

If the source does not fit, say so in one line and switch rather than forcing it.

## Gather input

**With no argument, the source is `chat`** - write up the current conversation.

Otherwise the argument can be any of these:

- **A description in words** - write the issue from that.
- **Review findings** - the findings reported earlier in this session (by `/review`,
  `/code-review`, `/simplify` or in the chat) or a pasted list. They become the checklist
  of the issue; the issue text names the pattern they share (see "One topic").
- **A branch name** (e.g. `16-port-fibertracking-cmake`) - read the change with
  `git -C mitk-diffusion diff master...<branch>` and `git -C mitk-diffusion log master..<branch>`.
  The change is normally the improvement already made, so reconstruct the issue *as it would
  have been written before the work started*: what was awkward before it and what that cost.
- **`uncommitted`** (or "working tree", "my changes", "wip") - read `git -C mitk-diffusion
  status`, `diff` and `diff --cached`. Same idea: describe the state these edits leave behind,
  not the edits.
- **`chat`** (or "this conversation", "our chat", "the session", or nothing) - use the
  current conversation: what was found, what was discussed, and what was decided. Cover the
  whole thread, not just the last message.

Sources combine: `chat uncommitted` means take the reasoning from the conversation and
confirm the places against the working tree. Where they disagree, the conversation wins - it
is what the user actually wanted.

A diff shows the result, never the reason. Never paste changed code into the issue; translate
the diff back into the state it repairs and the cost of that state. If a diff gives you no
plausible cost, say so and ask.

Prefer real material: real paths, the real test that covers the spot, the real MITK commit
that changed an interface. If the issue number is unknown, use `<nr>` - do not invent one.

Ask a question only if you cannot tell *what is awkward today* or *what it costs*.
Everything else: make a sensible assumption and note it in one line under the output.

## Name the cost

Every improvement issue names a concrete cost of the current state: a build that fails against
current MITK, a class of bug it invites, a change that took longer than it should, two places
that have drifted apart, a check a developer has to remember by hand, a manual that misleads.
"Cleaner", "more modern" or "best practice" are not costs. If you cannot name one, say so and
ask before writing - an issue without a cost will not be picked up.

## One topic, plus a short list of the rest

One issue covers one topic - "the modules pass a CMake argument MITK no longer knows", not
"everything the review found". Pick the one that costs the most and write the whole body about
it. Findings under that topic become checklist lines in `## Done when`, one per spot, each with
its path.

Smaller or unrelated points are **not** dropped and **not** told as a second story: list them
at the end under `## Other points`, one line each, at most three, no explanation beyond the
place and what is off. Anything beyond that goes in the one-line note under the output as
"also found, not written up". If two topics are equally big, ask which one, or write the
bigger one and say so.

## Describe direction, not a plan

`## Proposed change` says where to go in two or three sentences: what should exist instead
and what disappears. It is not a step-by-step plan, not a design document, no code. The reader
decides how. What must **not** change is part of the issue: the behaviour that stays the
same and the check that proves it (a `ctest` test, a manual check in the Workbench or with a
command line tool).

## Style

Short and easy to read beats complete. Both are rules, not goals.

- **Simple language.** Everyday words, one idea per sentence, no jargon where a common word
  works. A term that cannot be avoided gets half a sentence of explanation. Someone who has
  never opened this part of the repo must understand the issue on first read.
- **Short.** The body stays under roughly 200 words: current state 2-3 sentences, cost 1-2
  sentences, proposed change 2-3 sentences, 3-5 done-when lines, at most three other points.
  Cut by leaving things out, not by packing more into each sentence.
- Only the main points. No background essay, no history of how it got this way, no
  debugging narrative, no list of everything the reader could also look at.
- Paths, module, plugin and class names are welcome as evidence - the reader has to find the
  spot - but the story is told in words, not in code. No diffs, no function bodies.
- Present tense, current state.
- If the improvement would have to happen in MITK itself, it belongs to MITK: say so in the
  note and keep it text only.

## Save the issue

Write the title and body to `gitlab_issues/<file>.md` **before** replying, using
the Write tool. The file is the branch name plus `.md` - `gitlab_issues/16-port-fibertracking-cmake.md`
- and `new-<short-branch-name>.md` while the issue number is still `<nr>`. `gitlab_issues/` lies
outside both repositories, so this is a scratch copy, not a repo change.

The file holds exactly the title as a top-level `# ` heading, a blank line, and then the
issue body (label line included). The file is the only place the body appears - it is
**never** printed in the reply, not even in part. Commit message and branch name stay in the
reply only. If the file already exists, overwrite it.

## Issue body

Markdown, written to the file above:

```markdown
## Current state

What exists today and where, in 2-3 sentences. Name the places involved.

## Why it hurts

The concrete cost of leaving it as it is, in 1-2 sentences.

## Proposed change

What exists instead and what goes away, in 2-3 sentences.

## Done when

- [ ] Checkable outcome, with the path where it applies
- [ ] ...
- [ ] Behaviour X is unchanged, shown by <test or check>

## Other points

- <place>: <what is off>, one line

/label ~Refactoring
```

`## Other points` is optional: leave the heading out entirely when there is nothing besides
the main topic. It never grows past three one-line bullets - it exists so a small finding is
not lost, not so it can be discussed here.

The `/label ~Refactoring` line is a GitLab quick action and is the last line of the body - it
tags the issue on creation. Use only labels the project has (`create_gitlab_issue.py labels`);
if it has no fitting one, drop the line and say so in the note.

Done-when lines are **checkable**: a reviewer can tick each one by reading the named place,
running the named test or trying the Workbench. At least one line states what stays the
same and how that is shown. Not "code is cleaner".

## Create the issue in GitLab - only on request

The default is text only: write the file, reply, done. Create the issue only when the user's
request explicitly says so ("create it", "post it to GitLab", "open the issue", or the argument
`create`); asking for an issue, ticket or issue text is not enough.

Issues are created only in the mitk-diffusion project.

1. Look for an existing issue first:
   `python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py search <two or three key words>`.
   If one already covers the same topic, do not create - name it in the note.
2. Create it from the saved file:
   `python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py create gitlab_issues/new-<short-branch-name>.md`.
   The script turns the `/label` lines into labels, refuses labels the project does not have,
   renames the file to `<nr>-<short-branch-name>.md` and prints the number and link. It reads
   `GITLAB_TOKEN` from `mitk-diffusion/.env`; never print or pass the token yourself.
3. Put the new number into the branch name in place of `<nr>`.

## Output format

Reply with one line naming the file just written (after a GitLab create: the renamed file and
the issue link), then exactly these three blocks, each fenced so it can be copied in one click,
and nothing else except an optional one-line note about assumptions or findings left out. The
issue body is not among them - the reader opens the file.

**Title** - one line, names the topic and where it lives, no issue number, no type prefix:

```text
Short descriptive title
```

**Commit message** - one line, imperative, at most 50 characters, capitalized, no period (as in
`mitk/CLAUDE.md`): `Refactor ...`, `Port ...`, `Remove ...`, `Move ...`. No conventional-commit
prefixes:

```text
Refactor ...
```

**Branch name** - kebab-case, 3-5 words, with the issue number placeholder:

```text
<nr>-short-branch-name
```

If the user gave an issue number, put it in instead of `<nr>`.

## Notes

- The body above is this skill's own, kept in the shape of the `feature-request` skill.
- To *do* the review that feeds this skill, use `review`, `code-review` or `simplify`; this
  skill only writes the issue.
- Add further labels, assignees or milestones only if the user asks.
