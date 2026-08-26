"""Cross-receptor recurrence statistics, as submitted with the manuscript.

The site computes ΔRRCS per receptor from its own batch data, and that is what
every one of the 283 report pages shows. It cannot compute the cross-receptor
recurrence statistics the same way the manuscript did, because mapping a
residue to a shared generic position needs a GPCRdb number, and the site's
per-receptor annotation files only number residues that took part in a
high-magnitude contact. The manuscript analysis used the full GPCRdb residue
maps for all 283 receptors instead, which changes the answers:

    recurrent positions (CFR in >= 3 receptors)   site 356   manuscript 368
    3.50x50                                       site 261   manuscript 274
    CFR-CFR contact pair universe                 site 179   manuscript 2,234

So the statistics page reads the submitted tables rather than recomputing.
Provenance, including the source workbook's sha256, is in
data/manuscript_2026-08-24/PROVENANCE.md.

This module feeds the statistics page ONLY. Report pages keep their own
report-generation data and stay frozen against
tests/fixtures/freeze_manifest.json.
"""
import json
import pathlib

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from .cfr_analysis import CFR_TOP_N

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / 'data' / 'manuscript_2026-08-24'

TOP50_CSV = DATA_DIR / 'S3_top50_recurrent_cfr_positions.csv'
ALL_CSV = DATA_DIR / 'S4_all_recurrent_positions.csv'
PAIRS_CSV = DATA_DIR / 'S5_top20_recurrent_cfr_pairs.csv'
FIG3_JSON = DATA_DIR / 'figure3_enrichment.json'

#: Class colours, unchanged from the previous enrichment chart.
AM_COLORS = {'Pathogenic': '#E8820C', 'Ambiguous': '#F0C27A', 'Benign': '#008080'}


def _require(path):
    if not path.is_file():
        raise FileNotFoundError(
            '%s is missing. The statistics page needs the manuscript tables; '
            'they are tracked in the repository so a build never depends on '
            'the source workbook. Re-extract with '
            'scripts/extract_manuscript_tables.py.' % path)
    return path


def load_manuscript_stats():
    """Everything the statistics page needs, straight from the submitted tables."""
    top50 = pd.read_csv(_require(TOP50_CSV))
    all_positions = pd.read_csv(_require(ALL_CSV))
    pairs = pd.read_csv(_require(PAIRS_CSV))
    figure3 = json.loads(_require(FIG3_JSON).read_text(encoding='utf-8'))

    if len(top50) != CFR_TOP_N:
        raise ValueError(
            'the submitted top-position table has %d rows but the site reports '
            'a top %d. These must agree.' % (len(top50), CFR_TOP_N))

    # Preserve the submitted order exactly. Re-sorting here would silently
    # diverge from the published ranking, which is the whole point of reading
    # these files rather than recomputing.
    top50 = top50.sort_values('cfr_score_rank').reset_index(drop=True)

    return {
        'top50': top50,
        'all_positions': all_positions,
        'pairs': pairs,
        'figure3': figure3,
        'n_recurrent_positions': len(all_positions),
        'n_pairs_universe': figure3.get('n_pairs_universe'),
        'fig_recurrence_dotplot': make_recurrence_dotplot(top50),
        'fig_enrichment_bar': make_enrichment_bar(figure3),
    }


def make_recurrence_dotplot(top50):
    """The ranked recurrent positions, in the submitted order."""
    if top50.empty:
        return go.Figure()
    fig = px.scatter(
        top50, x='mean_abs_delta', y='n_receptors',
        size='cfr_score', color='segment',
        hover_data=['generic_number', 'cfr_score_rank', 'recurrence_rank',
                    'pct_receptors'],
        text='generic_number',
        title='Top %d recurrent Core Functional Residue positions' % CFR_TOP_N,
        labels={
            'mean_abs_delta': 'Mean |ΔRRCS|',
            'n_receptors': 'Receptors with an above-threshold change',
            'segment': 'Segment',
            'cfr_score': 'CFR score',
        },
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_traces(textposition='top center', textfont_size=8)
    fig.update_layout(height=600)
    return fig


def make_enrichment_bar(figure3):
    """Figure 3a: AlphaMissense class composition, recurrent against the rest.

    Ambiguous predictions are drawn but excluded from the odds ratio, exactly
    as the published figure does.
    """
    groups = figure3['panel_a']['groups']
    if not groups:
        return go.Figure()
    labels = ['%s<br>n = %s' % (g['label'], format(g['n'], ',')) for g in groups]
    fig = go.Figure()
    for cls in ('Pathogenic', 'Ambiguous', 'Benign'):
        key = cls.lower()
        fig.add_trace(go.Bar(
            name=cls, x=labels, y=[g[key] for g in groups],
            marker_color=AM_COLORS[cls],
            text=['%.1f%%' % g[key] for g in groups],
            textposition='inside',
        ))
    a = figure3['panel_a']
    fig.update_layout(
        barmode='stack',
        title=('AlphaMissense class at the top %d recurrent positions and at '
               'other generic-numbered positions (odds ratio %.2f)'
               % (CFR_TOP_N, a['odds_ratio'])),
        yaxis_title='Share of annotated missense substitutions (%)',
        height=450,
    )
    return fig
