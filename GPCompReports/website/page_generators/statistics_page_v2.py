"""Generate the v2 statistics page.

Sibling of statistics_page.py. Writes output/statistics_v2.html. Reuses the
same Plotly figures produced by analysis/* modules; theme overrides are
applied client-side via Plotly.relayout() driven by the data-chart-layout-*
attributes on each chart div.
"""

import json

from jinja2 import Environment

from ..plotly_theming import theme_overrides


def generate_statistics_page_v2(env: Environment, store, analysis_results, output_dir):
    cross = analysis_results.get('cross_gpcr', {})
    tm = analysis_results.get('tm_domain', {})
    cfr = analysis_results.get('cfr', {})
    variant = analysis_results.get('variant', {})

    chart_map = {
        'chart-ranking': cross.get('fig_ranking_bar'),
        'chart-histogram': cross.get('fig_delta_histogram'),
        'chart-box-ligand': cross.get('fig_box_ligand'),
        'chart-box-family': cross.get('fig_box_family'),
        'chart-tm-bar': tm.get('fig_tm_bar'),
        'chart-tm-heatmap': tm.get('fig_tm_heatmap'),
        'chart-cv-scatter': tm.get('fig_conserved_variable'),
        'chart-cfr-dotplot': cfr.get('fig_cfr_dotplot'),
        'chart-pathogenicity': variant.get('fig_pathogenicity_bar'),
        'chart-conservation': variant.get('fig_conservation_scatter'),
    }

    charts = {}
    for chart_id, fig in chart_map.items():
        if fig is not None:
            charts[chart_id] = fig.to_json()

    cfr_table = cfr.get('cfr_table')
    cfr_table_data = []
    if cfr_table is not None and not cfr_table.empty:
        cfr_table_data = cfr_table.head(30).to_dict('records')

    cfr_network = cfr.get('cfr_network')
    cfr_network_data = []
    if cfr_network is not None and not cfr_network.empty:
        cfr_network_data = cfr_network.head(20).to_dict('records')

    path_result = variant.get('pathogenicity', {})
    path_stats = path_result.get('stats', {})

    hi_variants = variant.get('high_impact_variants')
    hi_variant_data = []
    if hi_variants is not None and not hi_variants.empty:
        hi_variant_data = hi_variants.head(30).to_dict('records')

    light, dark = theme_overrides()

    template = env.get_template('statistics_v2.html')
    html = template.render(
        static_prefix='',
        active_page='statistics',
        nav_home_url='index_v2.html',
        nav_browse_url='browse/index_v2.html',
        nav_stats_url='statistics_v2.html',
        page_title='Database Statistics · GPCompReports',
        total_gpcrs=len(store.gpcr_ids),
        charts=charts,
        cfr_table=cfr_table_data,
        cfr_network=cfr_network_data,
        path_stats=path_stats if path_stats else None,
        cfr_pathogenic_pct=path_result.get('cfr_pathogenic_pct', 0),
        non_cfr_pathogenic_pct=path_result.get('non_cfr_pathogenic_pct', 0),
        high_impact_variants=hi_variant_data,
        layout_light_json=json.dumps(light, separators=(',', ':')),
        layout_dark_json=json.dumps(dark, separators=(',', ':')),
    )

    out_path = output_dir / 'statistics_v2.html'
    out_path.write_text(html, encoding='utf-8')
    print(f"  Generated: statistics_v2.html")
