#!/usr/bin/env python3
"""Generate the frozen-content manifest used by GPCompaReports_v2/tests/test_freeze.py.

This is a deliberate, human-run command. It is never invoked by
generate_site.py, build_demo.sh, or any test, and it never runs as a build
step. Its output is a git-tracked fixture
(GPCompaReports_v2/tests/fixtures/freeze_manifest.json) that pins the
PI-approved content of every report page, so that a routine full site
rebuild into the gitignored output/ directory can never silently become the
new baseline.

Usage:
    python3 scripts/build_freeze_manifest.py <source_dir> [--force]

<source_dir> must be a pristine, pre-V3 build (normally
GPCompaReports_v2/output/, the directory this project treats as the source
of truth for the frozen content). There is no default: the human running
this must name the directory, so that pointing it at the wrong build is a
deliberate choice, not an accident of a default flag.

The script refuses to run against anything that looks like a V3 build (see
V3_MARKERS below), and refuses to overwrite an existing manifest unless
given --force. Its output is a proposal, not a trusted artifact: review the
printed summary before `git add`-ing and committing the result.
"""
import argparse
import hashlib
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_PATH = (
    REPO_ROOT / 'GPCompaReports_v2' / 'tests' / 'fixtures' / 'freeze_manifest.json')

SCHEMA_VERSION = 1
BRANCH_BASE_COMMIT = '3a10a53'
GH_PAGES_DEPLOY = 'b892874'
BUILD_DATE = '2026-05-07'

# Any of these appearing anywhere under source_dir means it is a V3 build, or
# at least carries V3 markup, and must never be used to generate the freeze
# manifest. This is the guard against accidentally approving the very
# changes the freeze test exists to detect.
V3_MARKERS = (
    'v3-secnav',
    'v3-nav.js',
    'v3-deeplink',
    'v3.css',
    'sec-snake',
    'GPCompaRe database',
)

# "Approved content" = the inner HTML of every <section class="report-section">
# element, in document order, not including the opening tag. The V3 branch
# adds id= and data-section-title= to that opening tag deliberately (sanctioned
# wayfinding); everything inside the tag is the PI-approved content that must
# not move. This regex and SECTION_SEP must be kept identical to the copy in
# GPCompaReports_v2/tests/test_freeze.py, or the manifest and the test that
# reads it will silently disagree about what "approved content" means.
SECTION_RE = re.compile(
    r'<section class="report-section"[^>]*>(.*?)</section>', re.S)
SECTION_SEP = b'\x00SECTION\x00'

# Stylesheets whose pre-V3 fingerprint is recorded for provenance only. No
# test asserts on these; styling is explicitly out of scope for the freeze
# (see GPCompaReports_v2/tests/fixtures/README.md). Paths are relative to
# static/css/ in both source_dir and the current source tree.
STYLESHEETS = (
    'landing.css',
    'site.css',
    'primer/primitives.css',
    'primer/primitives-dark.css',
)


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def scan_for_v3_markers(source_dir):
    """Return {marker: [relative file paths]} for markers found under source_dir."""
    hits = {}
    for path in sorted(source_dir.rglob('*')):
        if not path.is_file():
            continue
        if path.suffix.lower() not in ('.html', '.css', '.js'):
            continue
        try:
            text = path.read_text(encoding='utf-8', errors='ignore')
        except OSError:
            continue
        for marker in V3_MARKERS:
            if marker in text:
                hits.setdefault(marker, []).append(str(path.relative_to(source_dir)))
    return hits


def approved_content_hash(html):
    sections = SECTION_RE.findall(html)
    blob = SECTION_SEP.join(section.encode('utf-8') for section in sections)
    return hashlib.sha256(blob).hexdigest(), len(sections)


def git_head():
    try:
        result = subprocess.run(
            ['git', 'rev-parse', 'HEAD'], cwd=str(REPO_ROOT),
            capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except Exception:
        return None


def build_manifest(source_dir):
    reports_dir = source_dir / 'reports'
    report_files = sorted(reports_dir.glob('*.html'))
    if not report_files:
        print('ERROR: no *.html files found under %s' % reports_dir, file=sys.stderr)
        sys.exit(1)

    receptors = {}
    for path in report_files:
        gpcr_id = path.stem
        html = path.read_text(encoding='utf-8')
        digest, count = approved_content_hash(html)
        if count == 0:
            print('ERROR: %s has zero <section class="report-section"> elements, '
                  'refusing to record a hash for it' % gpcr_id, file=sys.stderr)
            sys.exit(1)
        receptors[gpcr_id] = {'sha256': digest, 'section_count': count}

    src_css_dir = REPO_ROOT / 'GPCompaReports_v2' / 'static' / 'css'
    pre_css_dir = source_dir / 'static' / 'css'
    stylesheets = {}
    for rel in STYLESHEETS:
        pre_path = pre_css_dir / rel
        cur_path = src_css_dir / rel
        pre_hash = sha256_file(pre_path) if pre_path.exists() else None
        cur_hash = sha256_file(cur_path) if cur_path.exists() else None
        stylesheets[rel] = {
            'pre_v3_sha256': pre_hash,
            'current_source_sha256': cur_hash,
            'current_source_path': 'GPCompaReports_v2/static/css/%s' % rel,
            'changed': pre_hash != cur_hash,
        }
    # v3.css is new: it did not exist before the branch, so there is no
    # pre-V3 fingerprint to compare it against.
    v3css_cur = src_css_dir / 'v3.css'
    if v3css_cur.exists():
        stylesheets['v3.css'] = {
            'pre_v3_sha256': None,
            'current_source_sha256': sha256_file(v3css_cur),
            'current_source_path': 'GPCompaReports_v2/static/css/v3.css',
            'changed': True,
            'note': 'did not exist in the pre-V3 build',
        }

    return {
        'schema_version': SCHEMA_VERSION,
        # The definition of "approved content" is recorded here so the test can
        # assert its own copy matches. Without this, the generator and the test
        # each hold their own regex and separator, and a change to one would
        # silently redefine what is being frozen: the exact silent-failure mode
        # this manifest replaced.
        'content_definition': {
            'section_regex': SECTION_RE.pattern,
            'separator_hex': SECTION_SEP.hex(),
            'hash_algorithm': 'sha256',
            'scope': 'inner HTML of each <section class="report-section">, '
                     'in document order, excluding the opening tag',
        },
        'provenance': {
            'source_dir': str(source_dir),
            'branch_base_commit': BRANCH_BASE_COMMIT,
            'gh_pages_deploy': GH_PAGES_DEPLOY,
            'build_date': BUILD_DATE,
            'generated_by_commit': git_head(),
            'generated_at': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        },
        'receptors': receptors,
        'stylesheet_provenance': {
            '_note': ('Informational only. No test asserts on these hashes; '
                      'styling is explicitly out of scope for the freeze, see '
                      'GPCompaReports_v2/tests/fixtures/README.md.'),
            'files': stylesheets,
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description='Generate the report-content freeze manifest. A deliberate, '
                     'human-run command; never invoked by a build or a test.')
    parser.add_argument(
        'source_dir', type=str,
        help='Pristine pre-V3 build directory to hash, e.g. GPCompaReports_v2/output '
             '(no default: you must name it)')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite an existing manifest file')
    args = parser.parse_args()

    source_dir = Path(args.source_dir).resolve()
    if not source_dir.is_dir():
        print('ERROR: not a directory: %s' % source_dir, file=sys.stderr)
        sys.exit(1)

    print('==> Scanning %s for V3 markers...' % source_dir)
    hits = scan_for_v3_markers(source_dir)
    if hits:
        print('ERROR: refusing to generate a manifest from what looks like a '
              'V3 build.', file=sys.stderr)
        print('       Found markers:', file=sys.stderr)
        for marker, files in sorted(hits.items()):
            sample = ', '.join(files[:3])
            more = '' if len(files) <= 3 else ' (+%d more)' % (len(files) - 3)
            print('         %-20s in %s%s' % (marker, sample, more), file=sys.stderr)
        print('       The manifest must come from a pristine pre-V3 build, e.g. '
              'GPCompaReports_v2/output/. Do not point this at output_v3_demo/, '
              'output_v3_alt/, or any build you just made.', file=sys.stderr)
        sys.exit(1)
    print('    no V3 markers found')

    if MANIFEST_PATH.exists() and not args.force:
        print('ERROR: manifest already exists at %s' % MANIFEST_PATH, file=sys.stderr)
        print('       Pass --force to overwrite it.', file=sys.stderr)
        sys.exit(1)

    manifest = build_manifest(source_dir)

    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + '\n', encoding='utf-8')

    receptors = manifest['receptors']
    total_sections = sum(entry['section_count'] for entry in receptors.values())
    print('')
    print('==> Wrote %s' % MANIFEST_PATH)
    print('    receptors:      %d' % len(receptors))
    print('    total sections: %d' % total_sections)
    print('    first fingerprints:')
    for gpcr_id in sorted(receptors)[:5]:
        entry = receptors[gpcr_id]
        print('      %-16s sections=%d sha256=%s...'
              % (gpcr_id, entry['section_count'], entry['sha256'][:16]))
    print('')
    print('REMINDER: this file is a proposed baseline, not yet a trusted one.')
    print('          Review it by hand (diff against the previous manifest if')
    print('          there was one, sanity-check the receptor count and a few')
    print('          hashes) before `git add`-ing and committing it.')


if __name__ == '__main__':
    main()
