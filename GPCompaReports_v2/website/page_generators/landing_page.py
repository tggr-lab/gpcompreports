"""Generate the landing page for GPCompaReports v2."""

import json
from pathlib import Path

from jinja2 import Environment


def generate_landing_page(env: Environment, store, output_dir: Path):
    """Render the landing to output/index.html."""
    from ..page_generators import gpcr_report_helpers as helpers

    info_df = store.get_all_info_df()

    total_gpcrs = len(store.gpcr_ids)
    n_families = info_df['receptor_family'].nunique()

    n_contact_records = 0
    n_threshold_changes = 0
    for gid in store.gpcr_ids:
        delta_df = store.delta_data.get(gid)
        if delta_df is None or delta_df.empty:
            continue
        n_contact_records += len(delta_df)
        thr = helpers._calc_significance_threshold(delta_df)
        n_threshold_changes += int((delta_df['abs_delta'] >= thr).sum())

    n_models = 2 * len(store.gpcr_ids)
    n_excluded = len(store.metadata) - len(store.gpcr_ids)

    search_records = (
        info_df[['gpcr_id', 'uniprot_name', 'gene_name', 'receptor_family']]
        .fillna('')
        .to_dict('records')
    )

    template = env.get_template('landing.html')
    html = template.render(
        static_prefix='',
        active_page='home',
        nav_home_url='index.html',
        nav_browse_url='browse/index.html',
        nav_stats_url='statistics.html',
        nav_downloads_url='downloads.html',
        nav_contact_url='contact.html',
        page_title='GPCompaRe database: active-inactive contact changes across human Class A GPCRs',
        extra_css=['static/css/landing.css'],
        total_gpcrs=total_gpcrs,
        n_families=n_families,
        n_contact_records=n_contact_records,
        n_threshold_changes=n_threshold_changes,
        n_models=n_models,
        n_excluded=n_excluded,
        search_json=json.dumps(search_records, separators=(',', ':')),
    )

    out_path = output_dir / 'index.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: {out_path.name}")
