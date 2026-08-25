"""Generate the Contact page.

The form posts to Formspree. The endpoint is a placeholder until the lab creates
the account, see docs/superpowers/EXTERNAL_SETUP.md. Nothing here invents a
credential.
"""

FORMSPREE_ENDPOINT_PLACEHOLDER = 'FORMSPREE_ENDPOINT_NOT_YET_CONFIGURED'
CONTACT_EMAIL = 'tggrlab@gmail.com'


def generate_contact_page(env, store, output_dir):
    template = env.get_template('contact.html')
    html = template.render(
        static_prefix='',
        active_page='contact',
        nav_home_url='index.html',
        nav_browse_url='browse/index.html',
        nav_stats_url='statistics.html',
        nav_downloads_url='downloads.html',
        nav_contact_url='contact.html',
        page_title='Contact · GPCompaRe',
        total_gpcrs=len(store.gpcr_ids),
        formspree_endpoint=FORMSPREE_ENDPOINT_PLACEHOLDER,
        contact_email=CONTACT_EMAIL,
    )
    (output_dir / 'contact.html').write_text(html, encoding='utf-8')
    print("  Generated: contact.html")
