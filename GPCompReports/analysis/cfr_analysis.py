"""Core Functional Residue (CFR) identification from cross-GPCR significant changes."""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def run_cfr_analysis(store):
    """Identify Core Functional Residues and build CFR contact network."""
    results = {}

    cfr_table = identify_cfrs(store)
    results['cfr_table'] = cfr_table
    results['fig_cfr_dotplot'] = make_cfr_dotplot(cfr_table)

    cfr_network = build_cfr_network(store, cfr_table)
    results['cfr_network'] = cfr_network

    return results


def identify_cfrs(store):
    """For each generic number, count how many GPCRs show it as significant.

    Rank by frequency + mean abs_delta -> top 50 CFR positions.
    """
    gn_stats = {}  # generic_number -> {count, deltas, segments, gpcr_ids}
    total_gpcrs = 0

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        sig_df = store.significant_data.get(gid, pd.DataFrame())
        if sig_df.empty or not annot_map:
            continue

        total_gpcrs += 1

        # Build position -> generic_number mapping
        pos_gn = {}
        pos_seg = {}
        for pos, info in annot_map.items():
            gn = info.get('display_number', '') or info.get('generic_number', '')
            if gn:
                pos_gn[pos] = gn
                pos_seg[pos] = info.get('protein_segment', '')

        # Count significant residues at each generic number
        seen_gn = set()
        for _, row in sig_df.iterrows():
            for res_col in ['res1', 'res2']:
                pos = int(row[res_col])
                gn = pos_gn.get(pos, '')
                if gn and gn not in seen_gn:
                    seen_gn.add(gn)
                    if gn not in gn_stats:
                        gn_stats[gn] = {
                            'count': 0, 'deltas': [], 'segment': pos_seg.get(pos, ''),
                            'gpcr_ids': []
                        }
                    gn_stats[gn]['count'] += 1
                    gn_stats[gn]['gpcr_ids'].append(gid)

            # Also collect delta values for this generic number
            for res_col in ['res1', 'res2']:
                pos = int(row[res_col])
                gn = pos_gn.get(pos, '')
                if gn:
                    gn_stats[gn]['deltas'].append(abs(row['delta_rrcs']))

    # Build CFR table
    rows = []
    for gn, stats in gn_stats.items():
        if stats['count'] >= 3:  # At least 3 GPCRs
            rows.append({
                'generic_number': gn,
                'frequency': stats['count'],
                'frequency_pct': stats['count'] / total_gpcrs * 100 if total_gpcrs > 0 else 0,
                'mean_abs_delta': np.mean(stats['deltas']) if stats['deltas'] else 0,
                'max_abs_delta': np.max(stats['deltas']) if stats['deltas'] else 0,
                'segment': stats['segment'],
                'n_gpcrs': stats['count'],
            })

    cfr_df = pd.DataFrame(rows)
    if cfr_df.empty:
        return cfr_df

    # Composite score: normalize frequency and mean_abs_delta, then average
    max_freq = cfr_df['frequency'].max()
    max_delta = cfr_df['mean_abs_delta'].max()
    cfr_df['norm_frequency'] = cfr_df['frequency'] / max_freq if max_freq > 0 else 0
    cfr_df['norm_delta'] = cfr_df['mean_abs_delta'] / max_delta if max_delta > 0 else 0
    cfr_df['cfr_score'] = (cfr_df['norm_frequency'] + cfr_df['norm_delta']) / 2
    cfr_df = cfr_df.sort_values('cfr_score', ascending=False).reset_index(drop=True)
    cfr_df['rank'] = cfr_df.index + 1

    return cfr_df


def make_cfr_dotplot(cfr_table):
    """Dot plot: generic number vs frequency, sized by mean delta."""
    if cfr_table.empty:
        return go.Figure()

    top50 = cfr_table.head(50)

    fig = px.scatter(
        top50, x='mean_abs_delta', y='frequency',
        size='cfr_score', color='segment',
        hover_data=['generic_number', 'rank'],
        text='generic_number',
        title='Top 50 Core Functional Residue Positions',
        labels={
            'mean_abs_delta': 'Mean |ΔRRCS|',
            'frequency': 'Number of GPCRs with Significant Change',
            'segment': 'Segment',
            'cfr_score': 'CFR Score',
        },
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_traces(textposition='top center', textfont_size=8)
    fig.update_layout(height=600)
    return fig


def build_cfr_network(store, cfr_table):
    """Find contact pairs where both residues are CFR positions."""
    if cfr_table.empty:
        return pd.DataFrame()

    cfr_gns = set(cfr_table.head(50)['generic_number'].tolist())
    pair_counts = {}  # (gn1, gn2) -> count of GPCRs

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        sig_df = store.significant_data.get(gid, pd.DataFrame())
        if sig_df.empty or not annot_map:
            continue

        pos_gn = {}
        for pos, info in annot_map.items():
            gn = info.get('display_number', '') or info.get('generic_number', '')
            if gn:
                pos_gn[pos] = gn

        for _, row in sig_df.iterrows():
            gn1 = pos_gn.get(int(row['res1']), '')
            gn2 = pos_gn.get(int(row['res2']), '')
            if gn1 in cfr_gns and gn2 in cfr_gns and gn1 != gn2:
                pair = tuple(sorted([gn1, gn2]))
                pair_counts[pair] = pair_counts.get(pair, 0) + 1

    rows = []
    for (gn1, gn2), count in sorted(pair_counts.items(), key=lambda x: -x[1]):
        rows.append({'cfr_1': gn1, 'cfr_2': gn2, 'n_gpcrs': count})

    return pd.DataFrame(rows)
