---
name: investigate
description: Treat everything after /investigate as an investigation-only task - research and answer, make no code changes. Read-only mode stays on for the rest of the session until the user explicitly asks for an edit. Use when the user asks a question about the codebase, a behaviour, a cause or a design, and wants an answer rather than a change.
---

# Investigation only

Everything after `/investigate` is the question - with no argument, the one just asked in the
conversation. Answer it, change nothing.

## The rule

Until the user explicitly asks for an edit, change nothing the user owns: no `Edit`/`Write`,
no file-modifying shell (`sed -i`, redirects into repo paths, `rm`, `mv`), no git state changes
in `mitk-diffusion/` or `mitk/`, no CMake cache changes, no GitLab writes.

Allowed: reading, searching, read-only `git` and GitLab/web lookups, scratchpad notes, and
building or running `ctest` in an existing build tree (`debug/`, `release/`) when that is the
only way to answer (say so in the reply).

If answering would require a change, describe it and ask instead of making it. If a fix is
the obvious next step, sketch it in one short paragraph and ask whether to implement it. Do
not start.

## Depth

- Read the actual code - never answer from memory or a file name. Name the file and line the
  answer rests on.
- Find the cause, not the symptom; say so if it sits in MITK or a third-party library.
- When MITK changed an interface, find the MITK commit that did it
  (`git -C mitk log -S'<symbol>'`) and how MITK migrated its own code.
- Check whether an issue or `<issue-number>-…` branch already exists.

## Leaving the mode

Only an explicit "implement it" / "fix that" ends the mode. A follow-up question, agreement
with your analysis, or a fix looking small is not permission.
