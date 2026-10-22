#!/usr/bin/env bash
# Configures a Debug MITK + diffusion build tree for a mitk-diffusion worktree.
#
#   configure-worktree-build.sh <worktree> <build-dir> [extra cmake arguments]
#
# debug/ and release/ are bound to the main checkout through MITK_EXTENSION_DIRS, so a worktree
# needs a tree of its own. This one is the inner build only: it repeats the superbuild's
# MITK-Configure step with the worktree as extension and reuses the externals in debug/ep and
# the test data in debug/MITK-Data. Rerun it after the superbuild in debug/ was re-run.
# Extra arguments come last and override the initial cache, e.g.
# -DMITK_BUILD_CONFIGURATION=Fiberfox to check another build configuration.
set -euo pipefail

worktree=$(realpath "$1")
build=$(realpath -m "$2")
shift 2
[[ -f "$worktree/Modules/ModuleList.cmake" ]] || { echo "not a mitk-diffusion checkout: $worktree" >&2; exit 1; }

# The workspace root holds the main checkout next to mitk/ and debug/. The worktree's common git
# directory is <main checkout>/.git, wherever the worktree and this script live.
common_dir=$(git -C "$worktree" rev-parse --path-format=absolute --git-common-dir)
ws=$(dirname "$(dirname "$common_dir")")
[[ -f "$ws/debug/ep/tmp/MITK-Configure-cfgcmd.txt" ]] || { echo "no superbuild in $ws/debug" >&2; exit 1; }

mkdir -p "$build"
cache="$build/initial-cache.cmake"
sed "s|^set(MITK_EXTENSION_DIRS .*|set(MITK_EXTENSION_DIRS \"$worktree\" CACHE STRING \"Initial cache\" FORCE)|" \
  "$ws/debug/ep/tmp/MITK-Configure-cache-Debug.cmake" > "$cache"
grep -q "^set(MITK_EXTENSION_DIRS \"$worktree\"" "$cache"

# The -D arguments the superbuild passes to MITK-Configure (paths into debug/ep).
mapfile -t defs < <(sed -e "s/^cmd='//" -e "s/'$//" "$ws/debug/ep/tmp/MITK-Configure-cfgcmd.txt" \
  | tr ';' '\n' | grep '^-D')

cmake "${defs[@]}" -G "Unix Makefiles" -C "$cache" -S "$ws/mitk" -B "$build" "$@"
