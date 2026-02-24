"""Cross-GPCR comparison: rankings, distributions, clustering, family comparisons."""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_cross_gpcr_analysis(store):
    """Run all cross-GPCR analyses and return dict of results + plotly figures."""
    results = {}
    info_df = store.get_all_info_df()

    results['rankings'] = compute_rankings(info_df)
    results['fig_ranking_bar'] = make_ranking_bar(info_df)
    results['fig_delta_histogram'] = make_delta_histogram(store)
    results['fig_box_ligand'] = make_box_by_ligand(info_df)
    results['fig_box_family'] = make_box_by_family(info_df)

    pca_result = run_pca_analysis(store)
    results['pca_data'] = pca_result['data']
    results['fig_pca'] = pca_result['fig']
    results['pca_variance'] = pca_result['variance_explained']

    return results


def compute_rankings(info_df):
    """Rank GPCRs by sum_abs_delta, mean_abs_delta, max_abs_delta."""
    rankings = {}
    for metric in ['sum_abs_delta', 'mean_abs_delta', 'max_abs_delta']:
        ranked = info_df.sort_values(metric, ascending=False).reset_index(drop=True)
        ranked['rank'] = ranked.index + 1
        rankings[metric] = ranked[['gpcr_id', 'uniprot_name', 'gene_name',
                                    'receptor_family', 'ligand_type', metric, 'rank']].to_dict('records')
    # Add rank columns to info_df for later use
    for metric in ['sum_abs_delta', 'mean_abs_delta', 'max_abs_delta']:
        rank_col = f'rank_{metric}'
        info_df[rank_col] = info_df[metric].rank(ascending=False, method='min').astype(int)
    rankings['info_df_with_ranks'] = info_df
    return rankings


def make_ranking_bar(info_df):
    """Top 30 GPCRs bar chart by sum |delta|, colored by ligand type."""
    top30 = info_df.nlargest(30, 'sum_abs_delta').sort_values('sum_abs_delta', ascending=True)

    fig = px.bar(
        top30, x='sum_abs_delta', y='uniprot_name',
        color='ligand_type', orientation='h',
        title='Top 30 GPCRs by Total Conformational Change (Σ|ΔRRCS|)',
        labels={'sum_abs_delta': 'Sum |ΔRRCS|', 'uniprot_name': 'GPCR',
                'ligand_type': 'Ligand Type'},
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(
        height=700, yaxis_tickfont_size=10,
        margin=dict(l=100, r=20, t=50, b=50),
        legend=dict(orientation='h', y=-0.15),
    )
    return fig


def make_delta_histogram(store):
    """Global distribution histogram of all delta_rrcs values."""
    all_deltas = []
    for gid in store.gpcr_ids:
        df = store.delta_data.get(gid, pd.DataFrame())
        if not df.empty:
            all_deltas.extend(df['delta_rrcs'].tolist())

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=all_deltas, nbinsx=100,
        marker_color='#008080', opacity=0.8,
        name='ΔRRCS'
    ))
    fig.update_layout(
        title='Global Distribution of ΔRRCS Values Across All GPCRs',
        xaxis_title='ΔRRCS (Active − Inactive)',
        yaxis_title='Count',
        height=450,
    )
    return fig


def make_box_by_ligand(info_df):
    """Box plots of sum |delta| by ligand type."""
    fig = px.box(
        info_df, x='ligand_type', y='sum_abs_delta',
        color='ligand_type',
        title='Conformational Change Distribution by Ligand Type',
        labels={'sum_abs_delta': 'Sum |ΔRRCS|', 'ligand_type': 'Ligand Type'},
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(
        height=500, showlegend=False,
        xaxis_tickangle=-45,
        margin=dict(b=120),
    )
    return fig


def make_box_by_family(info_df):
    """Box plots by receptor family (top 15 by count)."""
    top_families = info_df['receptor_family'].value_counts().head(15).index.tolist()
    subset = info_df[info_df['receptor_family'].isin(top_families)]

    fig = px.box(
        subset, x='receptor_family', y='sum_abs_delta',
        color='receptor_family',
        title='Conformational Change by Receptor Family (Top 15)',
        labels={'sum_abs_delta': 'Sum |ΔRRCS|', 'receptor_family': 'Receptor Family'},
        color_discrete_sequence=px.colors.qualitative.Pastel,
    )
    fig.update_layout(
        height=550, showlegend=False,
        xaxis_tickangle=-45,
        margin=dict(b=150),
    )
    return fig


def run_pca_analysis(store):
    """PCA on per-GPCR TM-domain feature vectors, colored by family."""
    # Build feature matrix: for each GPCR, get per-TM-domain mean abs_delta
    tm_labels = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7',
                 'H8', 'ICL1', 'ICL2', 'ICL3', 'ECL1', 'ECL2', 'ECL3']
    rows = []
    gpcr_ids = []
    families = []
    info = store.gpcr_info

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        delta_df = store.delta_data.get(gid, pd.DataFrame())
        if delta_df.empty or not annot_map:
            continue

        # Map residues to segments
        tm_deltas = {tm: [] for tm in tm_labels}
        for _, row in delta_df.iterrows():
            for res_col in ['res1', 'res2']:
                pos = int(row[res_col])
                seg = annot_map.get(pos, {}).get('protein_segment', '')
                if seg in tm_deltas:
                    tm_deltas[seg].append(abs(row['delta_rrcs']))

        feature = [np.mean(tm_deltas[tm]) if tm_deltas[tm] else 0 for tm in tm_labels]
        rows.append(feature)
        gpcr_ids.append(gid)
        families.append(info.get(gid, {}).get('receptor_family', 'Unknown'))

    if len(rows) < 5:
        return {'data': None, 'fig': go.Figure(), 'variance_explained': []}

    X = np.array(rows)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=min(3, X_scaled.shape[1]))
    X_pca = pca.fit_transform(X_scaled)

    pca_df = pd.DataFrame({
        'PC1': X_pca[:, 0],
        'PC2': X_pca[:, 1],
        'gpcr_id': gpcr_ids,
        'receptor_family': families,
    })

    # Limit legend to top 10 families for readability
    top10 = pca_df['receptor_family'].value_counts().head(10).index.tolist()
    pca_df['family_display'] = pca_df['receptor_family'].apply(
        lambda x: x if x in top10 else 'Other'
    )

    var_explained = pca.explained_variance_ratio_

    fig = px.scatter(
        pca_df, x='PC1', y='PC2',
        color='family_display',
        hover_data=['gpcr_id'],
        title='PCA of GPCR Conformational Change Profiles',
        labels={
            'PC1': f'PC1 ({var_explained[0]:.1%} variance)',
            'PC2': f'PC2 ({var_explained[1]:.1%} variance)',
            'family_display': 'Receptor Family',
        },
        color_discrete_sequence=px.colors.qualitative.Set3,
    )
    fig.update_layout(height=600)

    return {
        'data': pca_df,
        'fig': fig,
        'variance_explained': var_explained.tolist(),
    }
