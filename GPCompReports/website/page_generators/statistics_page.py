"""Generate the statistics page with all cross-GPCR analysis charts."""

import json

from jinja2 import Environment


def generate_statistics_page(env: Environment, store, analysis_results, output_dir):
    """Generate statistics.html with all 12 interactive charts."""
    cross = analysis_results.get('cross_gpcr', {})
    tm = analysis_results.get('tm_domain', {})
    cfr = analysis_results.get('cfr', {})
    variant = analysis_results.get('variant', {})

    # Collect all chart JSONs
    charts = {}
    chart_map = {
        'chart-ranking': cross.get('fig_ranking_bar'),
        'chart-histogram': cross.get('fig_delta_histogram'),
        'chart-box-ligand': cross.get('fig_box_ligand'),
        'chart-box-family': cross.get('fig_box_family'),
        'chart-pca': cross.get('fig_pca'),
        'chart-tm-bar': tm.get('fig_tm_bar'),
        'chart-tm-heatmap': tm.get('fig_tm_heatmap'),
        'chart-cv-scatter': tm.get('fig_conserved_variable'),
        'chart-cfr-dotplot': cfr.get('fig_cfr_dotplot'),
        'chart-variant-violin': variant.get('fig_variant_violin'),
        'chart-pathogenicity': variant.get('fig_pathogenicity_bar'),
        'chart-conservation': variant.get('fig_conservation_scatter'),
    }

    for chart_id, fig in chart_map.items():
        if fig is not None:
            charts[chart_id] = fig.to_json()

    # CFR table data
    cfr_table = cfr.get('cfr_table')
    cfr_table_data = []
    if cfr_table is not None and not cfr_table.empty:
        cfr_table_data = cfr_table.head(30).to_dict('records')

    # CFR network data
    cfr_network = cfr.get('cfr_network')
    cfr_network_data = []
    if cfr_network is not None and not cfr_network.empty:
        cfr_network_data = cfr_network.head(20).to_dict('records')

    # Variant stats
    freq_stats = variant.get('freq_comparison', {}).get('stats', {})
    path_result = variant.get('pathogenicity', {})
    path_stats = path_result.get('stats', {})

    # High impact variants
    hi_variants = variant.get('high_impact_variants')
    hi_variant_data = []
    if hi_variants is not None and not hi_variants.empty:
        hi_variant_data = hi_variants.head(30).to_dict('records')

    # PCA variance
    pca_variance = cross.get('pca_variance', [])

    template = env.get_template('statistics.html')
    html = template.render(
        static_prefix='',
        active_page='statistics',
        total_gpcrs=len(store.gpcr_ids),
        charts=charts,
        cfr_table=cfr_table_data,
        cfr_network=cfr_network_data,
        freq_stats=freq_stats if freq_stats else None,
        path_stats=path_stats if path_stats else None,
        cfr_pathogenic_pct=path_result.get('cfr_pathogenic_pct', 0),
        non_cfr_pathogenic_pct=path_result.get('non_cfr_pathogenic_pct', 0),
        high_impact_variants=hi_variant_data,
        pca_variance=pca_variance,
    )

    out_path = output_dir / 'statistics.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: statistics.html")
