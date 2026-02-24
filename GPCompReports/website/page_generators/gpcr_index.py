"""Generate the GPCR browser index page and JSON data."""

import json
from pathlib import Path

from jinja2 import Environment


def generate_gpcr_index(env: Environment, store, output_dir):
    """Generate browse/index.html and data/gpcr_index.json."""
    # Generate JSON data for client-side filtering
    info_df = store.get_all_info_df()
    # Add rank
    info_df = info_df.sort_values('sum_abs_delta', ascending=False).reset_index(drop=True)
    info_df['rank'] = info_df.index + 1

    records = info_df[[
        'gpcr_id', 'uniprot_name', 'gene_name', 'receptor_family', 'ligand_type',
        'total_contacts', 'significant_changes', 'sum_abs_delta', 'mean_abs_delta',
        'max_abs_delta', 'variants_found', 'rank'
    ]].to_dict('records')

    # Write JSON data
    data_dir = output_dir / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)
    json_path = data_dir / 'gpcr_index.json'
    json_path.write_text(json.dumps(records, indent=None), encoding='utf-8')
    print(f"  Generated: data/gpcr_index.json ({len(records)} GPCRs)")

    # Generate browse HTML
    browse_dir = output_dir / 'browse'
    browse_dir.mkdir(parents=True, exist_ok=True)

    # Serialize JSON for inline embedding (avoids file:// fetch issues)
    gpcr_json = json.dumps(records, indent=None)

    template = env.get_template('gpcr_index.html')
    html = template.render(
        static_prefix='../',
        active_page='browse',
        total_gpcrs=len(store.gpcr_ids),
        gpcr_json=gpcr_json,
    )

    out_path = browse_dir / 'index.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: browse/index.html")
