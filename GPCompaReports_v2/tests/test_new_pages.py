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
    assert '.zip' not in html and '.tar.gz' not in html, 'a download link was invented'
    assert 'DOI pending release' in html
    assert 'zenodo.org/record' not in html, 'a Zenodo DOI was invented'


def test_both_pages_are_reachable_from_the_nav():
    index = (OUT / 'index.html').read_text(encoding='utf-8')
    assert 'contact.html' in index
    assert 'downloads.html' in index
