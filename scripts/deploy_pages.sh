#!/usr/bin/env bash
# Build the v2 site and publish it to the gh-pages branch.
#
# Usage:
#   scripts/deploy_pages.sh            # full build + deploy
#   scripts/deploy_pages.sh --limit 3  # quick sanity-deploy (3 reports only)
#
# Pre-requisites:
#   - `origin` remote pointing at the GitHub repo
#   - Push access to that repo
#   - The batch analysis CSV pipeline must have been run
#     (so The_batch_RRCS_analyzer/batch_analysis_full/ exists)
#   - Optional: conservation cache pre-populated via
#     scripts/fetch_conservation.py (otherwise snake-plot conservation view
#     falls back to variant-only positions)
#   - The browser smoke tests must be runnable: see
#     docs/superpowers/RELEASE_CHECKLIST.md for the one-time venv setup.
#     They run automatically below, against the freshly built site, before
#     anything is pushed. A skipped run is treated as a failure, not a pass,
#     and there is no flag to bypass them.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT="$ROOT/GPCompaReports_v2/output"
WORKTREE="$ROOT/.gh-pages-worktree"
BRANCH="gh-pages"

cd "$ROOT"

# 0. Sanity check: remote exists
if ! git remote | grep -qx "origin"; then
  echo "ERROR: no 'origin' remote. Set one with:" >&2
  echo "  git remote add origin git@github.com:tggr-lab/gpcompreports.git" >&2
  exit 1
fi

# 1. Build the site
echo "==> Building v2 site..."
python3 "$ROOT/GPCompaReports_v2/generate_site.py" "$@"

if [ ! -f "$OUTPUT/index.html" ]; then
  echo "ERROR: build did not produce $OUTPUT/index.html" >&2
  exit 1
fi

# 1b. Required pre-deployment check: browser smoke tests against the build
#     that is about to be published.
#
#     This runs after the build so it tests the real artefact, and before the
#     push so a failure costs nothing. run_browser_tests.sh exits non-zero if
#     the venv is missing, if Playwright is not installed, if there is no
#     build, or if pytest collects nothing, so a skip cannot masquerade as a
#     pass here.
#
#     There is deliberately no bypass flag. This is a publication companion
#     site, so no emergency justifies publishing unverified JavaScript. If an
#     exceptional case ever needs one, edit this script: that leaves a record
#     in git, which an environment variable would not.
echo "==> Required pre-deployment check: browser smoke tests..."
if ! SMOKE_BUILD_DIR=output bash "$ROOT/scripts/run_browser_tests.sh" -q; then
  echo "" >&2
  echo "ERROR: browser smoke tests did not pass against the build." >&2
  echo "       Nothing has been pushed. The gh-pages branch is untouched." >&2
  echo "" >&2
  echo "       If they could not run at all, that is still a failure, not a" >&2
  echo "       pass: see docs/superpowers/RELEASE_CHECKLIST.md for the" >&2
  echo "       one-time venv setup." >&2
  exit 1
fi

# 2. Clean any previous worktree
if [ -d "$WORKTREE" ]; then
  git worktree remove "$WORKTREE" --force || rm -rf "$WORKTREE"
fi

# 3. Ensure the gh-pages orphan branch exists locally
if ! git show-ref --verify --quiet "refs/heads/$BRANCH"; then
  echo "==> Creating orphan gh-pages branch..."
  git worktree add --detach "$WORKTREE"
  cd "$WORKTREE"
  git checkout --orphan "$BRANCH"
  git rm -rf . >/dev/null 2>&1 || true
  touch .nojekyll
  git add .nojekyll
  git commit -m "Initialize gh-pages"
  cd "$ROOT"
else
  git worktree add "$WORKTREE" "$BRANCH"
fi

# 4. Sync built output into the worktree. --delete drops files that no longer
#    exist in the build so the branch stays clean.
echo "==> Syncing output to gh-pages..."
rsync -a --delete \
  --exclude='.git' --exclude='.nojekyll' \
  --exclude='/data/' \
  --exclude='/poster/' --exclude='/dossier/' \
  "$OUTPUT/" "$WORKTREE/"
# poster/ and dossier/ are not part of the website. They are outputs of
# scripts/poster_figures.py and a one-off PAR dossier export that happen to
# live under output/, which is long-lived. Without these excludes rsync
# publishes them at guessable URLs even though no page links to them.
# data/ holds build-time JSON caches (conservation_*.json, alphamissense_*.json)
# that the HTML doesn't reference at runtime. Drop any stale data/ left from a
# previous deploy so gh-pages stays lean.
rm -rf "$WORKTREE/data"
touch "$WORKTREE/.nojekyll"

# 5. Commit + push
cd "$WORKTREE"
git add -A
if git diff --cached --quiet; then
  echo "==> No changes to deploy"
else
  TAG="$(cd "$ROOT" && git describe --tags --always)"
  STAMP="$(date -Iseconds)"
  git commit -m "Deploy: $TAG ($STAMP)"
  git push origin "$BRANCH"
fi

cd "$ROOT"
git worktree remove "$WORKTREE" --force
echo ""
echo "==> Done."
echo "    If not already enabled: GitHub repo Settings -> Pages -> Source = gh-pages / root"
echo "    Live URL: https://tggr-lab.github.io/gpcompreports/"
