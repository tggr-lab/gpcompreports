#!/usr/bin/env bash
# Run the report-page browser smoke tests.
#
# These need Playwright, which cannot go in the system Python: it is marked
# externally managed (PEP 668), and we do not modify the system install. So
# they live in a project venv and run as their own command. The ordinary
# suite is unaffected and still runs with:
#
#     python3 -m pytest GPCompaReports_v2/tests
#
# where these tests skip with a message rather than failing.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV="$ROOT/.venv-test"
BUILD="${SMOKE_BUILD_DIR:-output_v3_demo}"

if [ ! -x "$VENV/bin/python" ]; then
  echo "No test venv at $VENV. Create it once with:" >&2
  echo "    python3 -m venv \"$VENV\"" >&2
  echo "    \"$VENV/bin/python\" -m pip install -r GPCompaReports_v2/requirements-test.txt" >&2
  echo "    \"$VENV/bin/python\" -m playwright install chromium" >&2
  echo "  (pins playwright==1.62.0; the browser download is a separate step)" >&2
  exit 1
fi

if ! "$VENV/bin/python" -c "import playwright" 2>/dev/null; then
  echo "The venv exists but Playwright is missing. Install it with:" >&2
  echo "    \"$VENV/bin/python\" -m pip install -r GPCompaReports_v2/requirements-test.txt" >&2
  echo "    \"$VENV/bin/python\" -m playwright install chromium" >&2
  echo "  (pins playwright==1.62.0; the browser download is a separate step)" >&2
  exit 1
fi

if [ ! -d "$ROOT/GPCompaReports_v2/$BUILD/reports" ] && [ ! -d "$BUILD/reports" ]; then
  echo "No built reports to test at GPCompaReports_v2/$BUILD/reports." >&2
  echo "Build the demo first: bash scripts/build_demo.sh" >&2
  echo "Or point at another build: SMOKE_BUILD_DIR=output $0" >&2
  exit 1
fi

echo "==> Browser smoke tests against $BUILD"
cd "$ROOT/GPCompaReports_v2"
SMOKE_BUILD_DIR="$BUILD" exec "$VENV/bin/python" -m pytest \
  tests/test_smoke_browser.py "$@"
