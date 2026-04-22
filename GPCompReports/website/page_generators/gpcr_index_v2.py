"""Generate the v2 GPCR browser page.

Sibling of gpcr_index.py. Writes output/browse/index_v2.html alongside the
existing output/browse/index.html. Both versions reuse the same
output/data/gpcr_index.json (also produced by gpcr_index.py).
"""

import json

from jinja2 import Environment


def generate_gpcr_index_v2(env: Environment, store, output_dir):
    info_df = store.get_all_info_df()
    info_df = info_df.sort_values('sum_abs_delta', ascending=False).reset_index(drop=True)
    info_df['rank'] = info_df.index + 1

    records = info_df[[
        'gpcr_id', 'uniprot_name', 'gene_name', 'receptor_family', 'ligand_type',
        'total_contacts', 'significant_changes', 'sum_abs_delta', 'mean_abs_delta',
        'max_abs_delta', 'variants_found', 'rank'
    ]].to_dict('records')

    browse_dir = output_dir / 'browse'
    browse_dir.mkdir(parents=True, exist_ok=True)

    total_gpcrs = len(store.gpcr_ids)
    gpcr_json = json.dumps(records, separators=(',', ':'))

    template = env.get_template('gpcr_index_v2.html')
    html = template.render(
        static_prefix='../',
        active_page='browse',
        nav_home_url='../index_v2.html',
        nav_browse_url='index_v2.html',
        nav_stats_url='../statistics.html',
        page_title='Browse GPCRs · GPCompReports',
        total_gpcrs=total_gpcrs,
        gpcr_json=gpcr_json,
    )

    out_path = browse_dir / 'index_v2.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: browse/index_v2.html")
