#!/usr/bin/env bash
# Put the previously deployed site back, immediately.
#
#   bash scripts/rollback_pages.sh            # roll back to live-before-v3
#   bash scripts/rollback_pages.sh <commit>   # roll back to a specific deploy
#
# This force-updates the gh-pages branch on origin to a previous deploy commit.
# It does not rebuild anything and does not touch your working tree, so it is
# fast and cannot be affected by the state of your checkout.
#
# Nothing is lost: every deploy is a commit on gh-pages, so the version you are
# rolling back FROM stays in the branch history and can be restored the same
# way. `git log origin/gh-pages --oneline` lists them.
set -euo pipefail

TARGET="${1:-live-before-v3}"
BRANCH="gh-pages"

command -v git >/dev/null || { echo "git not found" >&2; exit 1; }
git remote | grep -qx origin || { echo "no 'origin' remote" >&2; exit 1; }

git fetch origin --tags --quiet

if ! COMMIT=$(git rev-parse --verify --quiet "${TARGET}^{commit}"); then
  echo "ERROR: '$TARGET' is not a commit or tag this repository knows." >&2
  echo "Recent deploys:" >&2
  git log origin/$BRANCH --oneline -10 >&2
  exit 1
fi

CURRENT=$(git rev-parse origin/$BRANCH)
echo "gh-pages is at : $CURRENT  $(git log -1 --format=%s "$CURRENT")"
echo "rolling back to: $COMMIT  $(git log -1 --format=%s "$COMMIT")"
if [ "$CURRENT" = "$COMMIT" ]; then
  echo "Already there. Nothing to do."
  exit 0
fi

printf "Type ROLLBACK to confirm: "
read -r reply
[ "$reply" = "ROLLBACK" ] || { echo "Aborted."; exit 1; }

git push origin "$COMMIT:$BRANCH" --force
echo ""
echo "Done. GitHub Pages usually republishes within a minute."
echo "The version you rolled back from is still in the branch history:"
echo "  git log origin/$BRANCH --oneline"
