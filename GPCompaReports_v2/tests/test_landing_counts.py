"""The landing page's scientific counts must be the frozen dataset's, rendered
server side, never animated from a placeholder.

Every intermediate frame of a count-up animation is a wrong number on a page whose
whole purpose is to be quoted, so the animation is gone and this pins it.
"""
import os
import pathlib
import re
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
_env = os.environ.get('FREEZE_BUILD_DIR')
if _env:
    BUILD_DIR = pathlib.Path(_env).expanduser()
    if not BUILD_DIR.is_absolute():
        BUILD_DIR = BASE / _env
else:
    BUILD_DIR = BASE / 'output_v3_demo'

INDEX = BUILD_DIR / 'index.html'
LANDING_JS = BASE / 'static' / 'js' / 'landing.js'

# Verified against the build on 2026-08-25. If a rebuild changes any of these,
# the dataset changed and that is a release event, not a test to update.
FROZEN = {
    '283': 'receptors',
    # '60' alone occurs many times in the page, so pin it to its own row instead.
    '213,456': 'contact-pair records',
    '23,025': 'threshold-passing changes',
    '566': 'models',
}


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_frozen_counts_appear_verbatim():
    html = INDEX.read_text(encoding='utf-8')
    for value, what in FROZEN.items():
        assert value in html, 'missing %s count: %s' % (what, value)
    # A bare '60' appears all over the page, so a substring check cannot pin the
    # family count. Anchor it to its own table row.
    assert re.search(r'Receptor families</th>\s*<td>60</td>', html), \
        'the receptor-family count is not 60'


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_the_promotional_stat_band_is_gone():
    html = INDEX.read_text(encoding='utf-8')
    assert 'stat-band' not in html, 'the stat band survived'
    assert 'stat-num' not in html, 'stat tiles survived'
    assert 'data-target' not in html, 'an animated counter survived'


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_the_limitation_statement_is_present():
    html = INDEX.read_text(encoding='utf-8')
    assert 'model-derived structural candidates' in html
    assert 'not direct experimental evidence of function' in html


def test_no_count_up_animation_remains():
    js = LANDING_JS.read_text(encoding='utf-8')
    assert 'data-target' not in js, 'the count-up animation is still wired up'
    assert 'animateCount' not in js, 'the count-up animation function survives'
