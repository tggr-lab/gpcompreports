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

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTPUT="$ROOT/GPCompReports_v2/output"
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
python3 "$ROOT/GPCompReports_v2/generate_site.py" "$@"

if [ ! -f "$OUTPUT/index.html" ]; then
  echo "ERROR: build did not produce $OUTPUT/index.html" >&2
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
rsync -a --delete --exclude='.git' --exclude='.nojekyll' "$OUTPUT/" "$WORKTREE/"
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
