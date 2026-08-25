"""The PI approved the v2 report pages. With every V3 toggle off, the text
inside the report sections must be what it was. This test compares a freshly
built demo report against the last full build in GPCompaReports_v2/output/.
"""
import re
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
OLD = BASE / 'output' / 'reports' / 'par2_human.html'
NEW = BASE / 'output_v3_demo' / 'reports' / 'par2_human.html'

SECTION_RE = re.compile(
    r'<section class="report-section"[^>]*>(.*?)</section>', re.S)
TAG_RE = re.compile(r'<[^>]+>')
WS_RE = re.compile(r'\s+')


def _section_text(html):
    out = []
    for body in SECTION_RE.findall(html):
        text = WS_RE.sub(' ', TAG_RE.sub(' ', body)).strip()
        out.append(text)
    return out


@pytest.mark.skipif(not NEW.exists(),
                    reason='no demo build on disk, run scripts/build_demo.sh')
def test_report_sections_are_unchanged_with_features_off():
    # Skipping when there is no demo build is fine: there is nothing to check.
    # Skipping when the demo EXISTS but the baseline does not is not fine. That
    # is the case where a green run would mean "we never looked" while reading
    # as "the freeze holds", on the one property this whole branch preserves.
    # output/ is gitignored, so this is the likely state on a fresh clone.
    assert OLD.exists(), (
        'demo build present but no baseline at %s. Run a full build first, '
        'this test cannot vouch for the freeze without one.' % OLD)
    old = _section_text(OLD.read_text(encoding='utf-8'))
    new = _section_text(NEW.read_text(encoding='utf-8'))
    assert old, 'baseline produced no sections, the regex or the fixture is wrong'
    assert len(new) == len(old), 'section count changed'
    for i, (a, b) in enumerate(zip(old, new)):
        assert a == b, 'section %d text changed' % i


V3_TOKEN_RE = re.compile(r'v3-[a-z-]+')


@pytest.mark.skipif(not NEW.exists(), reason='needs a demo build on disk')
def test_no_v3_markup_inside_approved_sections():
    """Structural half of the freeze rule.

    The text comparison above catches added words. This catches added markup
    that carries no text yet, which is how a V3 feature would most likely
    leak: an empty mount point or a hidden span inside approved content.
    """
    html = NEW.read_text(encoding='utf-8')
    for i, body in enumerate(SECTION_RE.findall(html)):
        found = V3_TOKEN_RE.findall(body)
        assert not found, 'section %d contains V3 markup: %s' % (i, sorted(set(found)))
