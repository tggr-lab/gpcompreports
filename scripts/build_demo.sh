#!/usr/bin/env bash
# Build the five-receptor V3 demo.
#
# Conservation and AlphaMissense caches live under the OUTPUT dir, so a fresh
# output directory silently loses the full snake-plot colour views. Seed them
# from the main build before generating.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/GPCompaReports_v2/output_v3_demo"
SRC_DATA="$ROOT/GPCompaReports_v2/output/data"
RECEPTORS="par2_human,par1_human,adrb2_human,5ht2a_human,cxcr4_human"

mkdir -p "$OUT/data"

# A silent seeding failure is the whole hazard: without these caches the snake
# plot's Conservation and AlphaMissense views fall from ~400 coloured positions
# to ~50, with no error and nothing on the page saying so. So count what we
# actually copied and stop if it is nothing.
echo "==> Seeding conservation / AlphaMissense caches..."
SEEDED=0
if [ -d "$SRC_DATA" ]; then
  for gid in ${RECEPTORS//,/ }; do
    for kind in conservation alphamissense; do
      if cp -f "$SRC_DATA/${kind}_${gid}.json" "$OUT/data/" 2>/dev/null; then
        SEEDED=$((SEEDED + 1))
      else
        echo "    missing: ${kind}_${gid}.json"
      fi
    done
  done
else
  echo "    no cache directory at $SRC_DATA"
fi

if [ "$SEEDED" -eq 0 ]; then
  echo "ERROR: seeded 0 cache files. The demo would build with degraded snake" >&2
  echo "       plots (Conservation and AlphaMissense views mostly empty) and" >&2
  echo "       would not represent the real site." >&2
  echo "       Run a full build first so $SRC_DATA exists, or re-run" >&2
  echo "       scripts/fetch_conservation.py and scripts/fetch_alphamissense.py." >&2
  exit 1
fi
echo "    seeded $SEEDED cache files"

echo "==> Building demo..."
python3 "$ROOT/GPCompaReports_v2/generate_site.py" --output "$OUT" --only "$RECEPTORS"

echo ""
echo "==> Done. Open:"
echo "    $OUT/index.html"

echo ""
echo "======================================================================"
echo " Demo report pages built (5):"
for gid in ${RECEPTORS//,/ }; do
  echo "    $OUT/reports/${gid}.html"
done
echo "======================================================================"
echo "The landing and browse pages list all 283 receptors, but only these five have pages in this demo, so links to the others will not resolve."
