"""Generate the Downloads page.

Neither archive exists yet. This module invents no file, no link, and no DOI:
it only names what each archive will contain once it is packaged. See
docs/superpowers/EXTERNAL_SETUP.md for what remains before either can ship.

The two placeholders are not equivalent. The database archive placeholder is
temporary: this page's companion database must ship a real archive at the
final public release, so this placeholder must not survive to publication.
The software placeholder may persist past publication, since there is no
releasable user-facing analysis program yet.
"""

RELEASE_NAME = 'GPCompaRe database release 2026.08'
SOFTWARE_NAME = 'GPCompaRe software v1.0.0'
DOI_STATUS = 'DOI pending release'


def generate_downloads_page(env, store, output_dir):
    template = env.get_template('downloads.html')
    html = template.render(
        static_prefix='',
        active_page='downloads',
        nav_home_url='index.html',
        nav_browse_url='browse/index.html',
        nav_stats_url='statistics.html',
        nav_downloads_url='downloads.html',
        nav_contact_url='contact.html',
        page_title='Downloads · GPCompaRe',
        total_gpcrs=len(store.gpcr_ids),
        release_name=RELEASE_NAME,
        software_name=SOFTWARE_NAME,
        doi_status=DOI_STATUS,
    )
    (output_dir / 'downloads.html').write_text(html, encoding='utf-8')
    print("  Generated: downloads.html")
