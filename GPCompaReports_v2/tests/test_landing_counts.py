"""The landing page's scientific counts must be the frozen dataset's, rendered
server side, never animated from a placeholder.

Every intermediate frame of a count-up animation is a wrong number on a page whose
whole purpose is to be quoted, so the animation is gone and this pins it.
"""
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
INDEX = BASE / 'output_v3_demo' / 'index.html'
LANDING_JS = BASE / 'static' / 'js' / 'landing.js'

# Verified against the build on 2026-08-25. If a rebuild changes any of these,
# the dataset changed and that is a release event, not a test to update.
FROZEN = {
    '283': 'receptors',
    '60': 'receptor families',
    '213,456': 'contact-pair records',
    '23,025': 'threshold-passing changes',
    '566': 'models',
}


@pytest.mark.skipif(not INDEX.exists(), reason='needs a demo build on disk')
def test_frozen_counts_appear_verbatim():
    html = INDEX.read_text(encoding='utf-8')
    for value, what in FROZEN.items():
        assert value in html, 'missing %s count: %s' % (what, value)


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
