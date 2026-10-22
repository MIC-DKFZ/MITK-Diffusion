---
name: bug-report
description: Write a short GitLab bug report for mitk-diffusion - title, markdown body (reproduction, expected, actual, system, logs, possible solution), commit message and a <nr>-branch-name branch; creates the issue in the mitk-diffusion GitLab project only when explicitly asked to. Use when the user asks for a bug report, bug ticket, or GitLab issue about something broken.
---

# GitLab bug report

Turn a failure (or a fix already made) into a short, copy-pasteable GitLab bug report for
`mic/personal/group6/neuro/mitk-diffusion` on git.dkfz.de. Paths are relative to the workspace
root (see `AGENTS.md`).

## Gather input

**With no argument, the source is `chat`** - write up the failure from the current conversation.

Otherwise the argument can be any of these:

- **A description in words** - write the report from that.
- **A branch name** (e.g. `12-fix-tck-reader`) - read the change with
  `git -C mitk-diffusion diff master...<branch>` and `git -C mitk-diffusion log master..<branch>`.
  The change is normally the *fix*, so reconstruct the bug it repairs: what was broken before
  it, and how someone would have hit that.
- **`uncommitted`** (or "working tree", "my changes", "wip") - read `git -C mitk-diffusion
  status`, `diff` and `diff --cached`. Same idea: describe the defect these edits fix, not the
  edits.
- **`chat`** (or "this conversation", "our chat", "the session", or nothing) - use the
  current conversation: the symptom that was reported, the commands actually run, the real
  error output, and what was expected instead. If a fix was made in the session, still write
  the report as of *before* the fix.

Sources combine: `chat uncommitted` means take the symptom from the conversation and confirm
the details against the working tree. Where they disagree, the conversation wins - it is what
actually happened.

A diff shows a fix, never a symptom. Never paste changed code, file names or function names
into the report as if they were the bug - translate them into what went wrong for someone using
a Workbench view or a command line tool. For a build or test failure, the symptom is the failing
command and its error. If a diff gives you no plausible symptom, say so and ask.

Prefer real material over invented material: real commands, real log lines, the real versions
(`git -C mitk-diffusion describe --always`, `git -C mitk describe`), the build type (`debug/` or
`release/`), the real OS. If a heading needs something you do not have, write a short bracketed
placeholder such as `<paste output of ...>` instead of guessing. If the issue number is unknown,
use `<nr>` - do not invent one.

Ask a question only if you cannot tell *what went wrong* or *what should have happened*.

## Redaction

Never put patient data, image or tractogram content, patient or study identifiers, credentials
or tokens into the report. Replace them with `<redacted>` in log excerpts and commands.

## One problem per report

Pick the single problem that hurts most and write the report about that one. Smaller things
found along the way are mentioned in one sentence inside the section they belong to - never
as their own heading, list or second story. If two problems are equally big and unrelated,
ask which one, or write the bigger one and say in the note that the other is left out.

## Describe the problem, not the fix

The report says what breaks and what should happen instead. Keep causes, code paths and
repair ideas out of the reproduction, expected and actual sections - a reader must be able to
confirm the bug without knowing how it will be fixed. The only place a fix appears is the last
section, and it stays a single sentence about the idea, no technical detail.

## Style

- Short. The whole body should fit on one screen: reproduction as a numbered list of a few
  steps, expected and actual one or two sentences each.
- Plain language a junior developer understands without this repo in their head. Avoid
  jargon; where a term is unavoidable, explain it in half a sentence.
- No walls of text, no background essay, no debugging history. Present tense, current broken
  behaviour.
- Steps, commands and logs stay precise and literal - plain language applies to the prose,
  not to a command.
- Logs as text in a fenced code block, trimmed to the relevant lines plus the command that
  produced them - never a screenshot, never a full dump.

## Save the report

Write the title and body to `gitlab_issues/<file>.md` **before** replying, using
the Write tool. The file is the branch name plus `.md` - `gitlab_issues/12-fix-tck-reader.md` -
and `new-<short-branch-name>.md` while the issue number is still `<nr>`. `gitlab_issues/` lies
outside both repositories, so this is a scratch copy, not a repo change.

The file holds exactly the title as a top-level `# ` heading, a blank line, and then the
issue body (label line included). The file is the only place the body appears - it is
**never** printed in the reply, not even in part. Commit message and branch name stay in the
reply only. If the file already exists, overwrite it.

## Issue body

Markdown, written to the file above:

```markdown
/label ~Bug

## How to reproduce

1. ...
2. ...

## Expected behavior

## Actual behavior

## System information

## Logs

## Possible solution
```

The label line is a GitLab quick action at the top of the body, so the issue is tagged on
creation. Use only labels the project has (`create_gitlab_issue.py labels`); if it has no
fitting one, drop the line and say so in the note.

`## Possible solution` is the only place the report talks about a fix. Ideally one sentence,
two at most, naming the idea in plain words - no file names, no function names, no code, no
step-by-step plan. Leave the heading empty if you have no plausible idea; do not invent one.

## Create the issue in GitLab - only on request

The default is text only: write the file, reply, done. Create the issue only when the user's
request explicitly says so ("create it", "post it to GitLab", "open the issue", or the argument
`create`); asking for a report, ticket or issue text is not enough.

Issues are created only in the mitk-diffusion project. A bug in MITK itself belongs to MITK:
keep it text only and say so in the note.

1. Look for an existing issue first:
   `python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py search <two or three key words>`.
   If one already describes the same bug, do not create - name it in the note.
2. Create it from the saved file:
   `python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py create gitlab_issues/new-<short-branch-name>.md`.
   The script turns the `/label` lines into labels, refuses labels the project does not have,
   renames the file to `<nr>-<short-branch-name>.md` and prints the number and link. It reads
   `GITLAB_TOKEN` from `mitk-diffusion/.env`; never print or pass the token yourself.
3. Put the new number into the branch name in place of `<nr>`.

## Output format

Reply with one line naming the file just written (after a GitLab create: the renamed file and
the issue link), then exactly these three blocks, each fenced so it can be copied in one click,
and nothing else except an optional one-line note about assumptions or missing information. The
issue body is not among them - the reader opens the file.

**Title** - one line, names the symptom and where it happens, no issue number:

```text
Short descriptive title
```

**Commit message** - one line, imperative, at most 50 characters, capitalized, no period (as in
`mitk/CLAUDE.md`). No conventional-commit prefixes:

```text
Fix ...
```

**Branch name** - kebab-case, 3-5 words, with the issue number placeholder:

```text
<nr>-short-branch-name
```

If the user gave an issue number, put it in instead of `<nr>`.

## Notes

- For a feature or improvement rather than a defect, use the `feature-request` or
  `improvement-request` skill instead.
