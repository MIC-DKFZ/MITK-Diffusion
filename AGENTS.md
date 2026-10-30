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
| `mitk/` | https://git.dkfz.de/mic/mitk.git | The base. **Do not change it.** |

The workspace root is not a git repository. Besides the build trees it holds `prompts/`
(handoff prompts for new sessions), `gitlab_issues/` (issue drafts and handoff plans from the
skills), `mitk-diffusion-test-data/` (a clone of the test data repository, see "Test data";
new references are committed there), `mitk-data/` (a clone of MITK-Data, which mitk-diffusion
no longer reads), a one-line `CLAUDE.md` importing this file and `.claude`, a symlink to
`mitk-diffusion/.claude`: Claude Code reads `.claude/settings.json` only from the folder a session
starts in, and loads skills from a subfolder only once a file there is touched. Edit both under
`mitk-diffusion/`. Check `git -C mitk-diffusion status --short --branch` before and after working.

Git worktrees go under `mitk-diffusion/.worktrees/` (git-ignored). A build tree is bound to one
`MITK_EXTENSION_DIRS` path, so `debug/` and `release/` build only the main checkout. A worktree
gets a Debug inner build of its own in `wt-build/<worktree name>`, configured by
`mitk-diffusion/.claude/scripts/configure-worktree-build.sh` (see "Commands"); it reuses the
externals in `debug/ep` and MITK's test data in `debug/MITK-Data`, and reads the diffusion test
data from `mitk-diffusion-test-data/`. `mitkConfig.h` embeds the extension path, so the first
build of a tree compiles nearly all of MITK: keep the tree as long as its worktree exists, and
rerun the script after a superbuild in `debug/`.

When MITK changed an interface, adapt mitk-diffusion to it. MITK usually migrated its own modules
in the same commit, so `git -C mitk log -S'<old symbol>'` shows the intended replacement. If a
change in `mitk/` looks unavoidable, stop and report what it would be and why.

MITK's `mitk/CLAUDE.md` and `mitk/.claude/rules/` (`modules.md`, `plugins.md`, `cmake.md`)
describe the module and plugin macros, `files.cmake`, tests and C++ conventions; they apply to
mitk-diffusion code too. Read the matching rule before touching CMake, a module or a plugin.

## Map

| Path | What it is |
| --- | --- |
| `mitk-diffusion/Modules/` | Modules, built in the order of `ModuleList.cmake`. `DiffusionImage`: a DWI is a plain `mitk::Image` plus gradient, b-value and measurement-frame properties (`Properties/mitkDiffusionPropertyHelper.h`), with NRRD/NIfTI/DICOM IO; `mitkMaskGeometry.h` holds the voxel-grid rule (`HasVoxelGridOf`, 1e-4 mm) of the IVIM, kurtosis and peak filters' masks. `DiffusionModelling`: tensor, ODF, SH and peak image types, their IO and mappers, reconstruction filters; in `Algorithms/Reconstruction/` the functions the reconstruction views and apps call (`mitkDiffusionReconstruction.h`: tensors, Q-ball, tensors to DWI; `mitkDiffusionPeaks.h`; `mitkDiffusionIndices.h`) and the ODF direction tables (`mitkPointShell.h`). `FiberBundle`: `mitk::FiberBundle` (a `vtkPolyData` of streamlines) and fib/vtk, tck, trk and DICOM tractography IO. `FiberProcessing`, `FiberTracking` (streamline + Gibbs), `MriSimulation` (Fiberfox), `FiberDissection` (OpenCV), `DiffusionPreprocessing`, `DiffusionRegistration` (MatchPoint), `IVIM`. |
| `mitk-diffusion/Modules/DiffusionCmdApps/` | Command line tools, one subfolder per group, each group switched by a `BUILD_Diffusion*CmdApps` option. A group's `CMakeLists.txt` is one call of `mitk_diffusion_create_cmdapps` (defined in `Modules/DiffusionCmdApps/CMakeLists.txt`) with the option, the modules of the group and its apps, `<App>` or `<App>:<module>,...` for modules only that app links. The apps link with `--no-as-needed`: a diffusion module registers its readers and writers when it is loaded, and an app that only loads and saves files uses none of its symbols. |
| `mitk-diffusion/Plugins/` | `org.mitk.gui.qt.diffusionimaging.*` views; `PluginList.cmake` registers them, all off by default. User manuals in `<plugin>/documentation/UserManual/`. atTRACTive (`...attractive`, on `FiberDissection`) is a plugin of its own, so the other fiber views (`...fiberprocessing`) build without OpenCV. |
| `mitk-diffusion/CMake/BuildConfigurations/` | `DiffusionDebug` (everything, tests on), `DiffusionRelease` (= DiffusionDebug + release settings), `Tractography`, `Fiberfox`, `IVIM`, `atTRACTive`, `DiffusionCmdApps` (no GUI); `Tractography`, `Fiberfox` and `IVIM` include the fiber processing views, only `DiffusionDebug`/`DiffusionRelease` and `atTRACTive` build OpenCV and atTRACTive. Each starts with `set_always.txt`, which resets every variable a configuration changes (external packages; every diffusion module, cmd-app group and plugin off), and then switches on what it contains; `set_release.txt` forces `Release` and `BUILD_TESTING OFF`. |
| `mitk-diffusion/CMakeExternals/`, `CMake/PackageDepends/` | OpenCV, which MITK no longer provides, added to the superbuild by the extension; `MITK-Diffusion-Data`, the test data (see "Test data"). |
| `mitk-diffusion/.claude/` | Skills (`bug-report`, `feature-request`, `improvement-request`, `handoff`, `investigate`, `review`), `settings.json`, `scripts/create_gitlab_issue.py`, which searches and creates issues in the GitLab project with `GITLAB_TOKEN` from `mitk-diffusion/.env` (git-ignored), and `scripts/configure-worktree-build.sh` for a worktree's build tree. |
| `mitk/` | MITK itself, see above. |
| `debug/` | Superbuild tree, Debug, `MITK_BUILD_CONFIGURATION=DiffusionDebug`; externals in `debug/ep/`, the MITK + diffusion build in `debug/MITK-build/`, MITK's test data checkout in `debug/MITK-Data/`. |
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
`set_always.txt` first resets every variable any configuration changes, so switching a tree to
another configuration gives the same result as a fresh tree. MITK evaluates some packages before
it applies the configuration, so the first configure after a switch can fail; the second is clean.

**Wiring for a new module, cmd-app group or plugin.** Module: its folder, an
`option(BUILD_MODULE_<Name>)` guard in its `CMakeLists.txt`, an entry in `ModuleList.cmake` after
its dependencies, an `OFF` line in `set_always.txt` and an `ON` line in each configuration that
should contain it. Cmd-app group: its subfolder whose `CMakeLists.txt` calls
`mitk_diffusion_create_cmdapps(OPTION BUILD_Diffusion<Group>CmdApps ...)`, an `add_subdirectory`
in `Modules/DiffusionCmdApps/CMakeLists.txt`, and the same `OFF` and `ON` lines; a new app is
one more name in its group's `APPS` list. Plugin: the same via `PluginList.cmake` and
`MITK_BUILD_<plugin>`. A forgotten line silently leaves it unbuilt. A plugin whose modules are
off fails the configure, but a module or cmd app with a missing dependency is only listed as
disabled in the configure's feature summary.

Facts the layout does not show:

- Modules, plugins and command line apps compile with `-Werror` unless they pass
  `WARNINGS_NO_ERRORS` (`mitk/CMake/mitkFunctionCreateModule.cmake`). Fix warnings rather than
  adding that option. GCC 13 at -O3 (`release/`) reports warnings the Debug trees never show, e.g.
  `-Warray-bounds` on an inlined `std::vector` copy; to check a change without a Release build,
  compile its translation unit with the flags of the target's `flags.make` in
  `release/MITK-build`, output to a scratch folder.
- MITK's superbuild fetches only the Boost libraries listed in `MITK_USE_Boost_HEADER_LIBRARIES`
  and `MITK_USE_Boost_COMPILED_LIBRARIES`; a feature adds its own with
  `mitkFunctionRequireBoostLibraries` (example: `MITK_USE_Python3` in `mitk/CMakeLists.txt`).
  mitk-diffusion requests its libraries in `CMakeExternals/ExternalProjectList.cmake`; a new
  Boost include needs its library added there before the next superbuild.
- Debug builds and `MITK_USE_Python3=ON` do not go together in MITK, so `set_always.txt` switches
  Python off; `DiffusionRelease` (which includes `DiffusionDebug`), `Tractography`, `Fiberfox`
  and `IVIM` switch it on for MITK's segmentation plugin. `DiffusionDebug` forces the build type
  `Debug`, `set_release.txt` forces `Release`.
- Configuring `FiberTracking` downloads a lookup-table tarball from mitk.org into its binary
  `resource/` folder, so it needs network. MITK embeds module resources only from `resource/`.
- Test data: module tests use `MITK_CREATE_MODULE_TESTS` / `mitkAddCustomModuleTest` in each
  `Testing/` folder and read their input and references from a checkout of
  https://git.dkfz.de/mic/personal/group6/neuro/mitk-diffusion-test-data (branch `main`; paths as
  below `DiffusionImaging/` in MITK-Data, which no test reads), named by `MITK_DIFFUSION_DATA_DIR`.
  The superbuild clones it at the commit pinned in `CMakeExternals/MITK-Diffusion-Data.cmake`
  unless the variable is given, and passes it to the inner build. A `Testing/CMakeLists.txt`
  includes `mitkDiffusionTestData` after `MITK_CREATE_MODULE_TESTS()`: it fails the configure
  when the variable is unset and gives C++ `mitk::GetDiffusionTestDataFilePath("<path>")`
  (`mitkDiffusionTestData.h`); CMake test arguments use `${MITK_DIFFUSION_DATA_DIR}/<path>`. New
  or changed references are committed in `mitk-diffusion-test-data/` (pushed only when the user
  asks; then raise the pinned commit). A changed reference gets a new file name
  (`<name>_v2.<ext>`, `_v3`, ...) and the old file stays, so other branches keep finding theirs.
- Gradient convention: the modified gradients (`DiffusionPropertyHelper::GetGradientContainer`),
  which all filters, mappers, trackers and Fiberfox use, are world (LPS) vectors,
  `modified = MF^T * original` with the stored measurement frame `MF` (frame axes as rows) and
  without the image rotation. NIfTI bvals/bvecs follow FSL, as MRtrix does (voxel frame, x negated
  when the image matrix has a positive determinant); NRRD uses the standard `measurement frame`
  field. NRRD files of earlier MITK versions (only MITK's `measurement frame:=` key) are read as
  before, with the image rotation applied. DICOM import is unchanged;
  `mitkDiffusionDicomSeriesTest` reads synthetic series of each vendor parser (`DicomDwi/`:
  Siemens mosaic and non-mosaic CSA2, GE, Philips; identity orientation; values in that
  folder's README). Untested: the gradient frame of oblique Siemens/Philips
  scans, GE b-values stored as b + 1e9, Siemens CSA1 and enhanced multi-frame series.
  `ctest -R 'mitkDiffusion(Property|Image|Dicom)'` runs the `DiffusionImage` tests of readers,
  writers and property helper.
- SH and peak files: MITK's SH basis (`mitk::sh`, `mrtrix = true`) is MRtrix 0.2's (functions
  with m != 0 are MRtrix 3's divided by sqrt(2)) at world (LPS) directions. `ShImageWriter` writes
  MRtrix 3's orthonormal basis and marks the file `MITK_SH_convention=MRtrix3` (NRRD key, NIfTI
  descrip); an image in MITK's FSL definition (`ShImage::FSL`) is written unchanged and marked
  `FSL`. `ShImageReader` converts files so marked or whose descrip starts with
  `MRtrix version: 3`, and reads unmarked files (the 4D SH files of earlier MITK versions) in the
  legacy basis. For unmarked files its option "SH convention" forces one; "FSL/Dipy" converts
  dipy's default basis (`descoteaux07`, legacy) from the FSL bvec frame, which is not MITK's FSL
  definition (taken unchanged). MITK's image reader stays the default for 4D SH NIfTI files
  (`MimeTypes_Sh4dNifti_ImageReaderIsDefault`); the SH reader ("SH Image") is chosen in the open
  dialog or with a `PreferenceListReaderOptionsFunctor`. That functor takes the first reader on
  its list in ascending confidence, not in list order, so a list naming "SH Image" and a
  `Supported` reader of the same file gets the SH reader. Older coefficient files with the
  coefficients as vector pixels of a 3D image
  (`qBallImage_SHCoeffs.nrrd`, `StreamlineTractography/qball_image_sh.nii.gz`) get the SH reader
  offered too (confidence `PartiallySupported`, NRRD DWIs excluded), which reads them, with
  "auto", in MITK's FSL definition, the basis of those MITK versions. NIfTI peak images are
  written in RAS with the descrip `MITK_peak_convention=RAS` and read back to LPS, as are MRtrix
  peak files (descrip `MRtrix ...`, reader option on by default); NRRD peaks are LPS.
- Reconstruction: `mitk::pointshell::Directions(n)` holds the ODF sampling directions as data;
  their order is part of the ODF file format (252 values per voxel) and of Fiberfox references.
  The single-shell Q-ball methods take volumes whose b-value rounds to 0 as b0 and need one
  non-zero shell; the multi-shell method needs b, 2b, 3b. The multi-shell Q-ball and radial
  multi-shell filters take their shells from the b-value map they are given, without one from
  `BValueMapProperty::CreateBValueMap`, the DWI's own map. `ReconstructQball` computes the sampled
  ODFs only when asked (`QballParameters::sampledOdf`). Peak extraction
  (`itk::OdfMaximaExtractionFilter`) keeps a refined maximum only if it is farther than the
  angular threshold from every larger one kept, so each axis appears once, and gives no peak
  where the sampled ODF varies by at most 1e-5 of its maximum (isotropic).
  `ctest -R 'ImageReconstruction|DiffusionReconstruction'` runs their tests.
- `ctest -R 'mitkDiffusion(Modelling|Reconstruction)|PeakShImageReader|ImageReconstruction'` runs
  the `DiffusionModelling` tests. `mitkDiffusionModellingMapperTest` renders offscreen (EGL, no
  display needed) and is skipped with exit code 77 where no OpenGL context can be created.
- The netlib L-BFGS-B code behind `vnl_lbfgsb` keeps its working variables in static storage.
  Threaded code holds `mitk::GetLbfgsbMutex()` (`Modules/DiffusionModelling/Algorithms/mitkLbfgsbMutex.h`)
  while `minimize` runs; the cost function may release it during its evaluation.
- `ctest -R mitkFiberBundle` runs the `FiberBundle` tests; `mitkFiberBundleMapperTest` renders
  offscreen (EGL) and is skipped with exit code 77 where no OpenGL context can be created.
  `FiberBundle` operations that keep, split or move fibers carry each fiber's weight along and
  give the same result at any number of OpenMP threads; fiber references in the test data
  depend on exact point positions, so a change there is checked bit for bit against the old output.
  A bundle keeps its fibers in the order of the polydata (file), one-point fibers as one-point
  lines, and drops fibers without points with a warning; `SetFiberPolyData` gives the new fibers
  weight 1. The mask cut (`RemoveFibersOutside`) cuts segments at the mask border;
  `ResamplePolylineToNumPoints` (equal arc length) serves `ResampleToNumPoints` and the
  clustering filter. The `FiberBundle` and `FiberProcessing` tests that change the OpenMP thread
  count use `mitk::OmpThreadGuard` (`CMake/mitkOmpThreadGuard.h`, in the include path of each
  test driver whose `Testing/CMakeLists.txt` includes `mitkDiffusionTestData`).
- `ctest -R 'mitkFiber(Transformation|Extraction|Processing|WeightsFile|DerivedImages|Clustering|ModificationFilter|ExtractionFilter|RoiImage)Test'`
  runs the `FiberProcessing` tests. Its filters read const bundles and give the same result at
  any number of OpenMP threads; bundles they derive keep the fiber weights, the RAS flag and a
  set TrackVis header (`mitkDerivedFiberBundle.h`), and the four derived-image filters share one
  voxel grid (`mitkFiberImageGrid.h`). The clustering filter prepares each fiber's metric data
  once (`ClusteringMetric::PrepareFiber`); the metric values equal the former vnl formulation
  bit for bit (`Metrics_SameValuesAsVnlFormulation`), which the clustering outputs depend on.
- `ctest -R Tractography` runs the `FiberTracking` tests. Streamline tracking seeds its random
  engines per seed point, so a tractogram with a fixed seed (`--fix_seed`, "Fix Random Seed") is
  the same at any number of OpenMP threads; the handlers take tensors, ODF directions and peaks as
  world (LPS) vectors and never apply the image rotation. Gibbs tracking samples the ODF along the
  particle in the image-axis frame and looks it up at the world direction. Tractography
  references in the test data (`StreamlineTractography/ReferenceFibs/`, `gibbsTractogram_v3.fib`)
  depend on the exact random engines and seeds, so a change there is checked bit for bit against
  the old output.
- `ctest -R 'FiberGeneration|Fiberfox|FiberFit'` runs the `MriSimulation` tests. Gradients in
  `.ffp` files (and their `.ffp.bvals`/`.ffp.bvecs`) are FSL bvecs of the simulated geometry;
  `FiberfoxParameters::ApplyDirectionMatrix` turns them into world vectors and
  `SetGradientsFromImage` stores a DWI's gradients in that form. Noise is complex Gaussian in
  k-space, or a noise model on the magnitude image: Rician or noncentral chi (`.ffp`
  `artifacts.noisetype` `noncentralchi` with `artifacts.noisecoils`; the `chisquare` type of earlier
  files is read as noncentral chi with the file's `numberofcoils`). With relaxation, an inversion
  time must be shorter than TR (`CheckInversionTime`); without, TI is not used. With a constant seed
  (`SetUseConstantRandSeed`, `MitkFiberfox --fix_seed`) a simulation gives the same image at any
  number of OpenMP threads (`mitkFiberfoxThreadDeterminismTest`: output, coil, phase and k-space
  images, spike log), and so does a fiber fit (`mitkFiberFitTest`; `--max_iter` counts L-BFGS-B
  iterations). The simulation log (status text, log file, progress) is `mitk::FiberfoxLog`, which
  reports progress through `itk::ProgressEvent`. Fiberfox references in the
  test data depend on the random engines and seeds, so a change there is checked bit for bit
  against the old output.
- `ctest -R mitkFiberDissection` runs the `FiberDissection` tests (label bundle moves of the
  labeling brush, sphere ROI selection, the classifier's candidates), all on synthetic bundles
  without a display. Label bundles are changed only through `StreamlineInteractorBrush::MoveFiber`,
  so fiber count and weights stay consistent; labeled fibers are found in the classifier's pool
  by their exact geometry. The module depends only on `FiberBundle` and OpenCV (`opencv_ml`);
  its one user is the atTRACTive plugin.

## Porting state

In the `DiffusionDebug` build, MITK and all of mitk-diffusion (modules, test drivers, command
line apps, plugins, user manuals) compile and link without warnings, and all 58 diffusion tests
pass. `DiffusionRelease` is built less often; its last full build stopped at the GCC 13 warning
named under "Facts", which is fixed but not rebuilt there. The Workbench has not been started.
Remove an item once it is fixed.

- MITK bugs stop the inner configure of the smaller configurations; mitk-diffusion itself
  configures in all of them. With `BUILD_TESTING` off (`set_release.txt`: `Tractography`,
  `Fiberfox`, `IVIM`, `atTRACTive`; MITK's own `Minimal` too),
  `mitk/Modules/QtWidgets/test/CMakeLists.txt:108` configures a test driver that does not exist.
  With Qt off (`DiffusionCmdApps`), the same happens there and in
  `mitk/Modules/RESTAPI/test/CMakeLists.txt:5-13` and `mitk/Modules/Agent/test/CMakeLists.txt:4`.
  A throwaway tree of `configure-worktree-build.sh` reuses `debug/ep`, which has no Python, so
  `Tractography`, `Fiberfox` and `IVIM` stop there earlier, at `FindPython3` (NumPy missing).
- MITK's `org.mitk.gui.qt.segmentation` does not compile with Python off: `MakeError` in
  `mitk/Plugins/org.mitk.gui.qt.segmentation/src/internal/QmitkSegmentationRestRoutes.cpp:465`
  is only used inside `#if MITK_HAS_PYTHON` (`-Werror=unused-function`). A bug in MITK.
  `DiffusionDebug` leaves the plugin off; `DiffusionRelease`, `Tractography`, `Fiberfox` and
  `IVIM` switch Python on and build it.

## Commands

Run from the workspace root. Both superbuild trees use Unix Makefiles; for Release, use
`release` in place of `debug`.

| Command | Effect |
| --- | --- |
| `cmake --build debug -j$(nproc)` | Superbuild: externals, then MITK + diffusion. Needed once, and after changing a `MITK_USE_*` option or an external. Takes long; every run repeats each external's update, configure, build and install steps and recompiles parts of some. |
| `cmake --build debug/MITK-build -j$(nproc)` | Rebuild MITK + diffusion: the normal edit-compile loop. |
| `cmake --build debug/MITK-build --target MitkFiberTracking -j$(nproc)` | One module (`Mitk<Module>`); a plugin target is its folder name with `_` for `.`, e.g. `org_mitk_gui_qt_diffusionimaging_tractography`. |
| `cmake debug/MITK-build` | Re-configure the inner build after editing CMake files in `mitk-diffusion/`. |
| `cmake -DMITK_DIFFUSION_DATA_DIR=$PWD/mitk-diffusion-test-data debug` | Let the superbuild pass the workspace clone of the test data instead of cloning its own; the inner build gets it at the next superbuild, or at once with the same `-D` on `debug/MITK-build`. |
| `ctest --test-dir debug/MITK-build -R Fiber` | Run tests matching a pattern; `-N` lists without running. Never call a test driver directly. |
| `debug/MITK-build/bin/MitkWorkbench` | Start the Workbench with the diffusion plugins; command line apps are in the same `bin/`. |
| `mitk-diffusion/.claude/scripts/configure-worktree-build.sh mitk-diffusion/.worktrees/<name> wt-build/<name>` | Configure a worktree's own build tree (about 20 s), and re-configure it after editing CMake files. Extra arguments override the initial cache, e.g. `-DMITK_BUILD_CONFIGURATION=Fiberfox` in a throwaway tree to check another configuration. Build it like `debug/MITK-build`. |
| `mkdir -p wt-build/<name>/tmp && TMPDIR=wt-build/<name>/tmp ctest --test-dir wt-build/<name> -R Fiber` | Tests in a worktree's tree. Failing tests write to the system temp folder; a `TMPDIR` per tree keeps parallel sessions apart. |

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
- Every handoff (a `prompts/` file or a `handoff` plan) is started in a new chat in the
  workspace root with `carry out this plan: <absolute path>`. A plan for work in a worktree
  therefore names the worktree and pins every edit, git command (`git -C <worktree>`), build
  and test to it by absolute path, and the main checkout stays untouched.

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
  input belongs in mitk-diffusion-test-data.
