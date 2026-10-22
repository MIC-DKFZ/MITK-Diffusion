---
name: feature-request
description: Write a short GitLab feature request for mitk-diffusion - title, markdown body (context, user story, acceptance criteria), commit message and a <nr>-branch-name branch; creates the issue in the mitk-diffusion GitLab project only when explicitly asked to. Use when the user asks for a feature request, GitLab issue, ticket, or issue text for a planned or finished change.
---

# GitLab feature request

Turn a rough idea (or a change already made) into a short, copy-pasteable GitLab feature request
for `mic/personal/group6/neuro/mitk-diffusion` on git.dkfz.de. Paths are relative to the
workspace root (see `AGENTS.md`).

## Gather input

**With no argument, the source is `chat`** - write up the current conversation.

Otherwise the argument can be any of three things:

- **An idea in words** - write the issue from that description.
- **The name of something already built** (a module, plugin view, command line tool, branch or
  directory - e.g. `FiberDissection`, `StreamlineTractography`). Find it, skim its user manual,
  UI strings and command line help, and reconstruct the feature request *as it would have been
  written before the work started*. Never describe the implementation back.
- **`chat`** (or "this conversation", "our chat", "the session", or nothing) - use the current
  conversation as the source: what the user asked for and what was actually built or
  decided in it. Cover the whole thread, not just the last message. If several unrelated
  things happened, ask which one to write up, or write the dominant one and say so.

Sources combine: `chat FiberDissection` means take the intent from the conversation and fill in
the details from the module. Where they disagree, the conversation wins - it is what the user
actually wanted; the built component only shows what came out of it.

Whatever the source, write from the **user's** point of view, not the implementation. If
the issue number is unknown, use `<nr>` as a placeholder - do not invent one.

Ask a question only if you cannot tell *who* wants the feature or *what problem it solves*.
Everything else: make a sensible assumption and note it in one line under the output.

## Style

- Simple language. A researcher who uses MITK Diffusion (the Workbench views or the command
  line tools) but has never read its code should understand it.
- Short. Context 2-4 sentences, user story one sentence, 3-6 acceptance criteria.
- Write from the perspective of that user, not a developer.
- Describe the wanted behaviour, not the solution. No class names, file names, CMake options or
  code unless the feature is genuinely about them.
- Present tense, no history ("the view does X" - not "we changed X to Y").
- Keep MITK and dMRI terms (view, data manager, tractogram, peak image) where a replacement
  would be unclear.

## Save the issue

Write the title and body to `gitlab_issues/<file>.md` **before** replying, using
the Write tool. The file is the branch name plus `.md` - `gitlab_issues/12-tract-envelope-export.md`
- and `new-<short-branch-name>.md` while the issue number is still `<nr>`. `gitlab_issues/` lies
outside both repositories, so this is a scratch copy, not a repo change.

The file holds exactly the title as a top-level `# ` heading, a blank line, and then the
issue body (label line included). The file is the only place the body appears - it is
**never** printed in the reply, not even in part. Commit message and branch name stay in the
reply only. If the file already exists, overwrite it.

## Issue body

Markdown, written to the file above:

```markdown
## Context

Why this is needed. What people do today and where it hurts.

## User story

As a <type of user>, I want to <goal> so that <benefit>.

## Acceptance criteria

- [ ] Observable thing a user can check
- [ ] ...

/label ~Feature
```

The `/label ~Feature` line is a GitLab quick action and is the last line of the body - it tags
the issue on creation. Use only labels the project has (`create_gitlab_issue.py labels`); if it
has no fitting one, drop the line and say so in the note.

Acceptance criteria are **functional**: what a person can see or do when the feature is
done. Not "class X has method Y". Each one must be checkable in the Workbench or with a command
line tool, not by reading code.

Tick a criterion (`- [x]`) when it is already fulfilled: the conversation built it, or the
component or branch in the source shows the behaviour exists in the current code. Leave it
open (`- [ ]`) when it is only planned, partly done, or you could not confirm it - a guess
stays unticked.

## Create the issue in GitLab - only on request

The default is text only: write the file, reply, done. Create the issue only when the user's
request explicitly says so ("create it", "post it to GitLab", "open the issue", or the argument
`create`); asking for a feature request, ticket or issue text is not enough.

Issues are created only in the mitk-diffusion project. A feature that needs MITK itself to
change belongs to MITK: keep it text only and say so in the note.

1. Look for an existing issue first:
   `python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py search <two or three key words>`.
   If one already asks for the same thing, do not create - name it in the note.
2. Create it from the saved file:
   `python3 mitk-diffusion/.claude/scripts/create_gitlab_issue.py create gitlab_issues/new-<short-branch-name>.md`.
   The script turns the `/label` lines into labels, refuses labels the project does not have,
   renames the file to `<nr>-<short-branch-name>.md` and prints the number and link. It reads
   `GITLAB_TOKEN` from `mitk-diffusion/.env`; never print or pass the token yourself.
3. Put the new number into the branch name in place of `<nr>`.

## Output format

Reply with one line naming the file just written (after a GitLab create: the renamed file and
the issue link), then exactly these three blocks, each fenced so it can be copied in one click,
and nothing else except an optional one-line note about assumptions. The issue body is not among
them - the reader opens the file.

**Title** - one line, plain, no issue number, no type prefix:

```text
Short descriptive title
```

**Commit message** - one line, imperative, at most 50 characters, capitalized, no period (as in
`mitk/CLAUDE.md`). No conventional-commit prefixes:

```text
Add ...
```

**Branch name** - kebab-case, 3-5 words, with the issue number placeholder:

```text
<nr>-short-branch-name
```

If the user gave an issue number, put it in instead of `<nr>`.

## Notes

- Add further labels, assignees or milestones only if the user asks.
