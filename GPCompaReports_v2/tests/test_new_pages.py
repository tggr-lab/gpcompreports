"""The Contact and Downloads pages must not ship invented credentials or dead links."""
import re
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / 'output_v3_demo'
CONTACT = OUT / 'contact.html'
DOWNLOADS = OUT / 'downloads.html'

pytestmark = pytest.mark.skipif(not OUT.exists(), reason='needs a demo build on disk')


def test_contact_form_has_every_required_field():
    html = CONTACT.read_text(encoding='utf-8')
    for field in ('name', 'email', 'institution', 'enquiry_type', 'message'):
        assert 'name="%s"' % field in html, 'missing field: %s' % field


def test_contact_form_accepts_no_file_uploads():
    html = CONTACT.read_text(encoding='utf-8')
    assert 'type="file"' not in html


def test_contact_page_invents_no_credentials():
    html = CONTACT.read_text(encoding='utf-8')
    assert 'formspree.io/f/' not in html, 'a Formspree endpoint was invented'
    assert not re.search(r'6L[\w-]{20,}', html), 'a reCAPTCHA site key was invented'


def test_downloads_page_offers_no_dead_links():
    html = DOWNLOADS.read_text(encoding='utf-8')
    lowered = html.lower()
    for ext in ('.zip', '.tar.gz', '.tgz', '.7z'):
        assert ext not in lowered, 'a download link was invented: %s' % ext
    assert 'DOI pending release' in html


def test_downloads_page_invents_no_doi():
    """Presence of the pending line is not enough: assert no real DOI sneaks in beside it."""
    html = DOWNLOADS.read_text(encoding='utf-8')
    assert 'zenodo.org' not in html.lower(), 'a Zenodo record was invented'
    assert not re.search(r'10\.\d{4,9}/\S+', html), 'a DOI-shaped string was invented'
    assert 'doi.org/' not in html.lower(), 'a resolvable DOI link was invented'


def test_both_pages_are_reachable_from_the_nav():
    """Reachability means the target exists, not that its name appears in the markup.

    The nav markup carried these strings before either page was built, so a
    substring check passes even when the page is missing. Assert the files.
    """
    index = (OUT / 'index.html').read_text(encoding='utf-8')
    for name, path in (('contact.html', CONTACT), ('downloads.html', DOWNLOADS)):
        assert name in index, '%s is not linked from the landing page' % name
        assert path.exists(), '%s is linked but was never built' % name
