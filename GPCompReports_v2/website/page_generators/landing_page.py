"""Generate the landing page for GPCompReports v2."""

import json
from pathlib import Path

from jinja2 import Environment

from ..format_helpers import DELTA_CLIP, SPARKLINE_BINS, bin_signed_delta


def generate_landing_page(env: Environment, store, output_dir: Path):
    """Render the landing to output/index.html."""
    info_df = store.get_all_info_df()

    total_gpcrs = len(store.gpcr_ids)
    n_families = info_df['receptor_family'].nunique()
    n_ligand_types = info_df['ligand_type'].nunique()
    total_contacts = int(info_df['total_contacts'].sum())
    total_contacts_k = f"{total_contacts // 1000}"

    top5 = info_df.nlargest(5, 'sum_abs_delta').to_dict('records')
    for g in top5:
        delta_df = store.delta_data.get(g['gpcr_id'])
        if delta_df is None or delta_df.empty:
            g['delta_bins'] = [0] * SPARKLINE_BINS
        else:
            g['delta_bins'] = bin_signed_delta(delta_df['delta_rrcs'].values)

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
        page_title='GPCompReports: active-inactive contact changes across human Class A GPCRs',
        extra_css=['static/css/landing.css'],
        total_gpcrs=total_gpcrs,
        n_families=n_families,
        n_ligand_types=n_ligand_types,
        total_contacts_k=total_contacts_k,
        top5=top5,
        search_json=json.dumps(search_records, separators=(',', ':')),
        sparkline_bins=SPARKLINE_BINS,
        delta_clip=DELTA_CLIP,
    )

    out_path = output_dir / 'index.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: {out_path.name}")
