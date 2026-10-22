# AGENTS.md

This workspace builds **MITK Diffusion**, a set of MITK modules, BlueBerry plugins and command
line tools for diffusion-weighted MRI, as an **extension of MITK**. This file is the operating
guide; it lives in the mitk-diffusion repository, and both this repo's `CLAUDE.md` and the
workspace root's import it, so put guidance here. Keep it short: state only what an agent cannot
read from the code, a README or MITK's own `CLAUDE.md`.

Paths below are relative to the **workspace root**, the folder that holds `mitk-diffusion/`,
`mitk/`, `debug/` and `release/` side by side; sessions start there.

## Project context

MITK Diffusion covers dMRI preprocessing, tensor/ODF/SH reconstruction, IVIM and kurtosis
fitting, streamline and global (Gibbs) tractography, fiber processing, and Fiberfox dMRI
simulation. It has no application of its own: its build configurations choose which diffusion
plugins the MITK Workbench is built with (all, or a smaller set such as Fiberfox or IVIM) and
which command line tools are built. Current MITK has no setting for the product name, so the
application is the "MITK Workbench". The code follows MITK's `develop` branch; it no longer
builds against the v2026.06 release. "Porting state" lists what is still open.

## Two repositories

| Path | Repository | Role |
| --- | --- | --- |
| `mitk-diffusion/` | https://git.dkfz.de/mic/personal/group6/neuro/mitk-diffusion, branch `master` | The project. All changes go here; GitLab issues are filed here. |
| `mitk/` | https://git.dkfz.de/mic/mitk.git, local branch `my-working` | The base. **Do not change it.** |

The workspace root is not a git repository. Besides the build trees it holds `prompts/`
(handoff prompts for new sessions), `gitlab_issues/` (issue drafts and handoff plans from the
skills), `mitk-data/` (a clone of https://git.dkfz.de/mic/mitk-data for committing test data;
the build trees read their own `MITK-Data` checkouts), a one-line `CLAUDE.md` importing this
file and `.claude`, a symlink to `mitk-diffusion/.claude`: Claude Code
reads `.claude/settings.json` only from the folder a session starts in, and loads skills from a
subfolder only once a file there is touched. Edit both under `mitk-diffusion/`. Check
`git -C mitk-diffusion status --short --branch` before and after working.

Git worktrees go under `mitk-diffusion/.worktrees/` (git-ignored). A build tree is bound to one
`MITK_EXTENSION_DIRS` path, so `debug/` and `release/` build only the main checkout.

When MITK changed an interface, adapt mitk-diffusion to it. MITK usually migrated its own modules
in the same commit, so `git -C mitk log -S'<old symbol>'` shows the intended replacement. If a
change in `mitk/` looks unavoidable, stop and report what it would be and why.

MITK's `mitk/CLAUDE.md` and `mitk/.claude/rules/` (`modules.md`, `plugins.md`, `cmake.md`)
describe the module and plugin macros, `files.cmake`, tests and C++ conventions; they apply to
mitk-diffusion code too. Read the matching rule before touching CMake, a module or a plugin.

## Map

| Path | What it is |
| --- | --- |
| `mitk-diffusion/Modules/` | Modules, built in the order of `ModuleList.cmake`. `DiffusionImage`: a DWI is a plain `mitk::Image` plus gradient, b-value and measurement-frame properties (`Properties/mitkDiffusionPropertyHelper.h`), with NRRD/NIfTI/DICOM IO. `DiffusionModelling`: tensor, ODF, SH and peak image types, their IO and mappers, reconstruction filters. `FiberBundle`: `mitk::FiberBundle` (a `vtkPolyData` of streamlines) and fib/vtk, tck, trk and DICOM tractography IO. `FiberProcessing`, `FiberTracking` (streamline + Gibbs), `MriSimulation` (Fiberfox), `FiberDissection` (OpenCV), `DiffusionPreprocessing`, `DiffusionRegistration` (MatchPoint), `IVIM`. |
| `mitk-diffusion/Modules/DiffusionCmdApps/` | Command line tools, one subfolder per group, each group switched by a `BUILD_Diffusion*CmdApps` option; an app is a name in its group's list in that `CMakeLists.txt`. |
| `mitk-diffusion/Plugins/` | `org.mitk.gui.qt.diffusionimaging.*` views; `PluginList.cmake` registers them, all off by default. User manuals in `<plugin>/documentation/UserManual/`. |
| `mitk-diffusion/CMake/BuildConfigurations/` | `DiffusionDebug` (everything, tests on), `DiffusionRelease` (= DiffusionDebug + release settings), `Tractography`, `Fiberfox`, `IVIM`, `atTRACTive`, `DiffusionCmdApps` (no GUI). Each starts with `set_always.txt`, which switches every diffusion module, cmd-app group and plugin off; `set_release.txt` forces `Release` and `BUILD_TESTING OFF`. |
| `mitk-diffusion/CMakeExternals/`, `CMake/PackageDepends/` | OpenCV, which MITK no longer provides, added to the superbuild by the extension. |
| `mitk-diffusion/.claude/` | Skills (`bug-report`, `feature-request`, `improvement-request`, `handoff`, `investigate`, `review`), `settings.json`, and `scripts/create_gitlab_issue.py`, which searches and creates issues in the GitLab project with `GITLAB_TOKEN` from `mitk-diffusion/.env` (git-ignored). |
| `mitk/` | MITK itself, see above. |
| `debug/` | Superbuild tree, Debug, `MITK_BUILD_CONFIGURATION=DiffusionDebug`; externals in `debug/ep/`, the MITK + diffusion build in `debug/MITK-build/`, test data checkout in `debug/MITK-Data/`. |
| `release/` | The same as a Release build, `DiffusionRelease`. |

## How the build works

`MITK_EXTENSION_DIRS` points the MITK superbuild at `mitk-diffusion/`. MITK then takes from the
extension: `CMake/` (module path), `CMakeExternals/ExternalProjectList.cmake` (superbuild
dependencies), `CMake/BuildConfigurations/*.cmake` (offered as `MITK_BUILD_CONFIGURATION`),
`CMake/PackageDepends/`, `Modules/ModuleList.cmake` and `Plugins/PluginList.cmake` (see the
`MITK_EXTENSION_DIR` loops in `mitk/CMakeLists.txt`). The superbuild builds the externals, then
configures and builds `MITK-build` itself (its `MITK-Configure` step).

The selected build configuration is re-included on **every** configure and `FORCE`s its cache
values, so a variable it sets cannot be changed in the cache; change the configuration file.

**Wiring for a new module or plugin.** Module: its folder, an `option(BUILD_MODULE_<Name>)` guard
in its `CMakeLists.txt`, an entry in `ModuleList.cmake` after its dependencies, an `OFF` line in
`set_always.txt` and an `ON` line in each configuration that should contain it. Plugin: the same
via `PluginList.cmake` and `MITK_BUILD_<plugin>`. A forgotten line silently leaves it unbuilt.

Facts the layout does not show:

- Modules, plugins and command line apps compile with `-Werror` unless they pass
  `WARNINGS_NO_ERRORS` (`mitk/CMake/mitkFunctionCreateModule.cmake`). Fix warnings rather than
  adding that option.
- MITK's superbuild fetches only the Boost libraries listed in `MITK_USE_Boost_HEADER_LIBRARIES`
  and `MITK_USE_Boost_COMPILED_LIBRARIES`; a feature adds its own with
  `mitkFunctionRequireBoostLibraries` (example: `MITK_USE_Python3` in `mitk/CMakeLists.txt`).
  mitk-diffusion requests its libraries in `CMakeExternals/ExternalProjectList.cmake`; a new
  Boost include needs its library added there before the next superbuild.
- Debug builds and `MITK_USE_Python3=ON` do not go together in MITK, so `DiffusionDebug` forces
  Python off; `DiffusionRelease` includes it and switches Python on again. Its build type is
  forced to `Debug`, `set_release.txt` forces `Release`.
- Configuring `FiberTracking` downloads a lookup-table tarball from mitk.org into its binary
  `resource/` folder, so it needs network. MITK embeds module resources only from `resource/`.
- Module tests use `MITK_CREATE_MODULE_TESTS` / `mitkAddCustomModuleTest` in each `Testing/`
  folder and read input from `${MITK_DATA_DIR}/DiffusionImaging/`.

## Porting state

In the `DiffusionDebug` and `DiffusionRelease` builds, MITK and all of mitk-diffusion (modules,
test drivers, command line apps, plugins, user manuals) compile and link without warnings, and
all 14 diffusion tests pass. The Workbench has not been started. Remove an item once it is fixed.

- Three regenerated test references are committed in `mitk-data/` (`c962623`) but not pushed;
  `debug/MITK-Data` and `release/MITK-Data` hold them as uncommitted copies until upstream has
  them:
  `DiffusionImaging/Fiberfox/uniform.fib` and `gaussian.fib` (MITK clamps ellipse control points
  to the plane since `git -C mitk show bbc5c77b27`, which moves a point of `Fiducial_2.pf`) and
  `DiffusionImaging/gibbsTractogram.fib` (from 2014; Gibbs tracking changed since, e.g. its random
  generator in 2019).
- The two `DiffusionRegistration` tests are compiled but deliberately not registered; they were
  commented out before the 2024 test refactoring.
- MITK bugs stop the inner configure of the smaller configurations; mitk-diffusion itself
  configures in all of them. With `BUILD_TESTING` off (`set_release.txt`: `Tractography`,
  `Fiberfox`, `IVIM`, `atTRACTive`; MITK's own `Minimal` too),
  `mitk/Modules/QtWidgets/test/CMakeLists.txt:108` configures a test driver that does not exist.
  With Qt off (`DiffusionCmdApps`), the same happens there and in
  `mitk/Modules/RESTAPI/test/CMakeLists.txt:5-13` and `mitk/Modules/Agent/test/CMakeLists.txt:4`.
- MITK's `org.mitk.gui.qt.segmentation` does not compile with Python off: `MakeError` in
  `mitk/Plugins/org.mitk.gui.qt.segmentation/src/internal/QmitkSegmentationRestRoutes.cpp:465`
  is only used inside `#if MITK_HAS_PYTHON` (`-Werror=unused-function`). A bug in MITK.
  `DiffusionDebug` switches the plugin off, `DiffusionRelease` (Python on) keeps it;
  `Tractography`, `Fiberfox` and `IVIM` switch it on without switching Python on.
- `FiberfoxParameters` (`Modules/MriSimulation/mitkFiberfoxParameters.h`) deep-copies its signal
  models on copy construction but shares them on assignment, and never deletes them.

## Commands

Run from the workspace root. Both superbuild trees use Unix Makefiles; for Release, use
`release` in place of `debug`.

| Command | Effect |
| --- | --- |
| `cmake --build debug -j$(nproc)` | Superbuild: externals, then MITK + diffusion. Needed once, and after changing a `MITK_USE_*` option or an external. Takes long; every run repeats each external's update, configure, build and install steps and recompiles parts of some. |
| `cmake --build debug/MITK-build -j$(nproc)` | Rebuild MITK + diffusion: the normal edit-compile loop. |
| `cmake --build debug/MITK-build --target MitkFiberTracking -j$(nproc)` | One module (`Mitk<Module>`); a plugin target is its folder name with `_` for `.`, e.g. `org_mitk_gui_qt_diffusionimaging_tractography`. |
| `cmake debug/MITK-build` | Re-configure the inner build after editing CMake files in `mitk-diffusion/`. |
| `ctest --test-dir debug/MITK-build -R Fiber` | Run tests matching a pattern; `-N` lists without running. Never call a test driver directly. |
| `debug/MITK-build/bin/MitkWorkbench` | Start the Workbench with the diffusion plugins; command line apps are in the same `bin/`. |

## Working rules

### Changes and bugs

- Fix the cause, not the symptom, and report the cause. When the cause is in MITK, say so and
  state what the change in mitk-diffusion does and does not cover.
- For every bug, first look for an existing issue and a branch named `<issue-number>-…` in the
  mitk-diffusion GitLab project and build on it. GitLab is read-only: never create or edit
  anything there; print issue text for the user to post. One exception: when the user
  explicitly asks an issue skill (`bug-report`, `feature-request`, `improvement-request`) to
  create the issue, it creates it in the mitk-diffusion project.
- Port to the current MITK API instead of wrapping old calls in compatibility shims. Keep each
  change small: no collateral refactoring or formatting churn in files touched only for the port.
- Keep this file current when structure, commands or the porting state change.
- Long tasks go in major steps; ask for confirmation after each. When a session grows long, hand
  over at a checkpoint: write `prompts/NN-<topic>.md` (next number) with the state, decisions,
  anchored facts and remaining steps, and give the user the start sentence.

### Git: never commit or push unasked

- Commit only when asked, push only when asked. Otherwise leave the work in the working tree
  and say what changed. A commit message or issue text is text: print it, do not run it.
- Commit message: follow `mitk/CLAUDE.md` (imperative subject up to 50 characters, blank line,
  body with what and why). No `Co-Authored-By` or "generated with" lines. Stage only the task's
  files, never `git add -A`. Never force-push; no tags, merges or rebases unasked.

### Code and comments

- New code follows the C++ conventions in `mitk/CLAUDE.md`; in existing files match the
  surrounding code. No tabs, newline at end of file. New files get the header from
  `mitk/Templates/copyright_header.txt`.
- Comment non-obvious changes with *why*, in terms of this project, describing the current
  state; git is the changelog.

### Tests and docs

- After a change, build the affected targets and run their tests with `ctest`. A failing test
  means a bug in the code, unless its assertion is wrong; never loosen a test to get it green.
- If a build or test step is skipped (time, missing data, no display), say exactly what was not
  run.
- Update a plugin's `documentation/UserManual/` when its view changes behaviour.

### Data

- Never add patient data, DWI/NIfTI/NRRD/DICOM files or tractograms to the repository; test
  input belongs in MITK-Data.
