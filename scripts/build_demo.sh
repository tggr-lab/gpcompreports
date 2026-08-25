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
if [ -d "$SRC_DATA" ]; then
  echo "==> Seeding conservation / AlphaMissense caches..."
  for gid in ${RECEPTORS//,/ }; do
    cp -f "$SRC_DATA/conservation_${gid}.json"  "$OUT/data/" 2>/dev/null || true
    cp -f "$SRC_DATA/alphamissense_${gid}.json" "$OUT/data/" 2>/dev/null || true
  done
fi

echo "==> Building demo..."
python3 "$ROOT/GPCompaReports_v2/generate_site.py" --output "$OUT" --only "$RECEPTORS"

echo ""
echo "==> Done. Open:"
echo "    $OUT/index.html"
