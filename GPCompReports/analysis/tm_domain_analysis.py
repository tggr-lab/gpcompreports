"""TM domain pattern analysis: per-domain changes, domain-domain interactions, conserved regions."""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats


TM_SEGMENTS = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']
ALL_SEGMENTS = TM_SEGMENTS + ['H8', 'ICL1', 'ICL2', 'ICL3', 'ECL1', 'ECL2', 'ECL3']


def run_tm_domain_analysis(store):
    """Run TM domain analyses and return results + figures."""
    results = {}

    domain_stats = compute_domain_stats(store)
    results['domain_stats'] = domain_stats
    results['fig_tm_bar'] = make_tm_bar(domain_stats)

    interaction_matrix = compute_domain_interactions(store)
    results['interaction_matrix'] = interaction_matrix
    results['fig_tm_heatmap'] = make_interaction_heatmap(interaction_matrix)

    cv_data = compute_generic_number_variation(store)
    results['cv_data'] = cv_data
    results['fig_conserved_variable'] = make_conserved_variable_scatter(cv_data)

    return results


def compute_domain_stats(store):
    """Per-TM-domain mean |delta| aggregated across all GPCRs."""
    domain_deltas = {seg: [] for seg in ALL_SEGMENTS}

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        delta_df = store.delta_data.get(gid, pd.DataFrame())
        if delta_df.empty or not annot_map:
            continue

        for _, row in delta_df.iterrows():
            # Assign interaction to the segment of each residue
            for res_col in ['res1', 'res2']:
                pos = int(row[res_col])
                seg = annot_map.get(pos, {}).get('protein_segment', '')
                if seg in domain_deltas:
                    domain_deltas[seg].append(abs(row['delta_rrcs']))

    stats_list = []
    for seg in ALL_SEGMENTS:
        vals = domain_deltas[seg]
        if vals:
            stats_list.append({
                'segment': seg,
                'mean_abs_delta': np.mean(vals),
                'median_abs_delta': np.median(vals),
                'total_abs_delta': np.sum(vals),
                'count': len(vals),
                'std_abs_delta': np.std(vals),
            })
        else:
            stats_list.append({
                'segment': seg, 'mean_abs_delta': 0, 'median_abs_delta': 0,
                'total_abs_delta': 0, 'count': 0, 'std_abs_delta': 0,
            })

    return pd.DataFrame(stats_list)


def make_tm_bar(domain_stats):
    """Bar chart of per-TM-domain mean |delta| (TM1-7 only)."""
    df = domain_stats[domain_stats['segment'].isin(TM_SEGMENTS)].copy()
    df = df.sort_values('mean_abs_delta', ascending=True)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=df['mean_abs_delta'], y=df['segment'],
        orientation='h',
        marker_color='#008080',
        text=[f'{v:.2f}' for v in df['mean_abs_delta']],
        textposition='outside',
    ))
    fig.update_layout(
        title='Mean Conformational Change by TM Helix',
        xaxis_title='Mean |ΔRRCS|',
        yaxis_title='Segment',
        height=400,
        margin=dict(l=80, r=80),
    )
    return fig


def compute_domain_interactions(store):
    """Compute domain-domain interaction change matrix (TM1-7 only)."""
    seg_list = TM_SEGMENTS
    seg_idx = {s: i for i, s in enumerate(seg_list)}
    n = len(seg_list)
    matrix = np.zeros((n, n))
    counts = np.zeros((n, n))

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        delta_df = store.delta_data.get(gid, pd.DataFrame())
        if delta_df.empty or not annot_map:
            continue

        for _, row in delta_df.iterrows():
            seg1 = annot_map.get(int(row['res1']), {}).get('protein_segment', '')
            seg2 = annot_map.get(int(row['res2']), {}).get('protein_segment', '')
            if seg1 in seg_idx and seg2 in seg_idx:
                i, j = seg_idx[seg1], seg_idx[seg2]
                val = abs(row['delta_rrcs'])
                matrix[i, j] += val
                matrix[j, i] += val
                counts[i, j] += 1
                counts[j, i] += 1

    # Convert to mean
    with np.errstate(divide='ignore', invalid='ignore'):
        mean_matrix = np.where(counts > 0, matrix / counts, 0)

    return pd.DataFrame(mean_matrix, index=seg_list, columns=seg_list)


def make_interaction_heatmap(interaction_matrix):
    """Heatmap of domain-domain interaction changes."""
    fig = go.Figure(data=go.Heatmap(
        z=interaction_matrix.values,
        x=interaction_matrix.columns.tolist(),
        y=interaction_matrix.index.tolist(),
        colorscale=[[0, '#F5F7FA'], [0.5, '#008080'], [1, '#E8820C']],
        text=np.round(interaction_matrix.values, 2),
        texttemplate='%{text}',
        textfont=dict(size=10),
        hovertemplate='%{y} - %{x}<br>Mean |ΔRRCS|: %{z:.2f}<extra></extra>',
    ))
    fig.update_layout(
        title='TM Helix Interaction Change Heatmap',
        height=500, width=550,
        xaxis_title='Segment', yaxis_title='Segment',
        yaxis_autorange='reversed',
    )
    return fig


def compute_generic_number_variation(store):
    """Coefficient of variation of delta per generic number across GPCRs."""
    gn_values = {}  # generic_number -> list of abs_delta values

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        delta_df = store.delta_data.get(gid, pd.DataFrame())
        if delta_df.empty or not annot_map:
            continue

        # Build position -> generic_number mapping
        pos_gn = {}
        for pos, info in annot_map.items():
            gn = info.get('display_number', '') or info.get('generic_number', '')
            if gn:
                pos_gn[pos] = gn

        for _, row in delta_df.iterrows():
            for res_col in ['res1', 'res2']:
                pos = int(row[res_col])
                gn = pos_gn.get(pos, '')
                if gn:
                    if gn not in gn_values:
                        gn_values[gn] = []
                    gn_values[gn].append(abs(row['delta_rrcs']))

    # Compute stats per generic number, filter to TM1-7 only
    tm_segment_ids = {'1', '2', '3', '4', '5', '6', '7'}
    rows = []
    for gn, vals in gn_values.items():
        if len(vals) >= 5:  # Need enough data points
            # Extract segment number from generic number (e.g., "3.50x50" -> "3")
            seg_id = gn.split('.')[0].replace('x', '') if '.' in gn else ''
            if seg_id not in tm_segment_ids:
                continue

            mean_val = np.mean(vals)
            std_val = np.std(vals)
            cv = std_val / mean_val if mean_val > 0 else 0
            rows.append({
                'generic_number': gn,
                'mean_abs_delta': mean_val,
                'std_abs_delta': std_val,
                'cv': cv,
                'n_observations': len(vals),
                'segment': f'TM{seg_id}',
            })

    return pd.DataFrame(rows)


def make_conserved_variable_scatter(cv_data):
    """Scatter plot: mean_delta vs CV for each generic number position (TM1-7 only)."""
    if cv_data.empty:
        return go.Figure()

    fig = px.scatter(
        cv_data, x='mean_abs_delta', y='cv',
        color='segment',
        hover_data=['generic_number', 'n_observations'],
        title='Conserved vs Variable Positions (Mean |Δ| vs Coefficient of Variation)',
        labels={
            'mean_abs_delta': 'Mean |ΔRRCS|',
            'cv': 'Coefficient of Variation',
            'segment': 'TM Segment',
        },
        color_discrete_sequence=px.colors.qualitative.Set2,
    )

    # Add Spearman correlation + trendline
    r, p = stats.spearmanr(cv_data['mean_abs_delta'], cv_data['cv'])
    fig.add_annotation(
        x=0.02, y=0.98, xref='paper', yref='paper',
        text=f"Spearman r = {r:.3f}, p = {p:.2e}",
        showarrow=False, font=dict(size=11),
        xanchor='left', yanchor='top',
    )

    # Add OLS trendline
    m, b = np.polyfit(cv_data['mean_abs_delta'], cv_data['cv'], 1)
    x_range = np.linspace(cv_data['mean_abs_delta'].min(), cv_data['mean_abs_delta'].max(), 50)
    fig.add_trace(go.Scatter(
        x=x_range, y=m * x_range + b,
        mode='lines', line=dict(color='gray', dash='dash', width=1.5),
        showlegend=False, hoverinfo='skip',
    ))

    fig.update_layout(height=550)
    # Annotate quadrants
    fig.add_annotation(x=0.95, y=0.05, xref='paper', yref='paper',
                       text='High Δ, Low CV: Conserved functional', showarrow=False,
                       font=dict(size=10, color='#008080'))
    fig.add_annotation(x=0.95, y=0.95, xref='paper', yref='paper',
                       text='High Δ, High CV: Variable functional', showarrow=False,
                       font=dict(size=10, color='#E8820C'))
    return fig
