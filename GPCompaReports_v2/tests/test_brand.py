from pathlib import Path

from GPCompaReports_v2.website.brand import BRAND_NAME, BRAND_SHORT

TEMPLATES = Path(__file__).resolve().parent.parent / 'templates'
BASE = Path(__file__).resolve().parent.parent.parent


def test_brand_constants():
    assert BRAND_NAME == 'GPCompaRe database'
    assert BRAND_SHORT == 'GPCompaRe'


def test_no_old_brand_text_left_in_templates():
    for path in TEMPLATES.rglob('*.html'):
        assert 'GPCompaReports' not in path.read_text(), path


def test_urls_and_hosts_are_untouched():
    base_html = (TEMPLATES / 'base.html').read_text()
    assert 'gpcompreports.goatcounter.com' in base_html
    deploy = (BASE / 'scripts' / 'deploy_pages.sh').read_text()
    assert 'tggr-lab.github.io/gpcompreports/' in deploy
