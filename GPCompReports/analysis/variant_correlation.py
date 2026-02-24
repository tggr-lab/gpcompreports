"""Variant correlation analysis: CFR vs gnomAD variant frequency and pathogenicity."""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from scipy import stats


def run_variant_analysis(store, cfr_table):
    """Run variant correlation analyses against CFR positions."""
    results = {}

    cfr_positions = _get_cfr_position_map(store, cfr_table)
    results['cfr_positions'] = cfr_positions

    freq_result = compare_allele_frequencies(store, cfr_positions)
    results['freq_comparison'] = freq_result
    results['fig_variant_violin'] = make_variant_violin(freq_result)

    path_result = pathogenicity_enrichment(store, cfr_positions)
    results['pathogenicity'] = path_result
    results['fig_pathogenicity_bar'] = make_pathogenicity_bar(path_result)

    conservation_result = conservation_vs_delta(store)
    results['conservation'] = conservation_result
    results['fig_conservation_scatter'] = make_conservation_scatter(conservation_result)

    high_impact = find_high_impact_variants(store, cfr_positions)
    results['high_impact_variants'] = high_impact

    return results


def _get_cfr_position_map(store, cfr_table):
    """Build mapping: for each GPCR, which positions are CFR positions.

    Returns dict: gpcr_id -> set of positions that map to top-50 CFR generic numbers.
    """
    if cfr_table.empty:
        return {}

    cfr_gns = set(cfr_table.head(50)['generic_number'].tolist())
    result = {}

    for gid in store.gpcr_ids:
        annot_map = store.get_annotation_map(gid)
        cfr_pos = set()
        for pos, info in annot_map.items():
            gn = info.get('display_number', '') or info.get('generic_number', '')
            if gn in cfr_gns:
                cfr_pos.add(pos)
        if cfr_pos:
            result[gid] = cfr_pos

    return result


def compare_allele_frequencies(store, cfr_positions):
    """Compare mean allele frequency at CFR positions vs non-CFR (Mann-Whitney U)."""
    cfr_afs = []
    non_cfr_afs = []

    for gid in store.gpcr_ids:
        var_df = store.variant_data.get(gid, pd.DataFrame())
        if var_df.empty:
            continue

        cfr_pos = cfr_positions.get(gid, set())
        for _, row in var_df.iterrows():
            af = row.get('af', None)
            if af is None or pd.isna(af):
                continue
            try:
                af = float(af)
            except (ValueError, TypeError):
                continue

            protein_pos = int(row['protein_position']) if pd.notna(row.get('protein_position')) else None
            if protein_pos is None:
                continue

            if protein_pos in cfr_pos:
                cfr_afs.append(af)
            else:
                non_cfr_afs.append(af)

    # Mann-Whitney U test
    stat_result = {}
    if len(cfr_afs) >= 5 and len(non_cfr_afs) >= 5:
        u_stat, p_value = stats.mannwhitneyu(cfr_afs, non_cfr_afs, alternative='two-sided')
        stat_result = {
            'test': 'Mann-Whitney U',
            'u_statistic': float(u_stat),
            'p_value': float(p_value),
            'cfr_median_af': float(np.median(cfr_afs)),
            'non_cfr_median_af': float(np.median(non_cfr_afs)),
            'cfr_n': len(cfr_afs),
            'non_cfr_n': len(non_cfr_afs),
        }

    return {
        'cfr_afs': cfr_afs,
        'non_cfr_afs': non_cfr_afs,
        'stats': stat_result,
    }


def make_variant_violin(freq_result):
    """Violin plot: allele frequency at CFR vs non-CFR positions."""
    cfr_afs = freq_result.get('cfr_afs', [])
    non_cfr_afs = freq_result.get('non_cfr_afs', [])

    if not cfr_afs and not non_cfr_afs:
        return go.Figure()

    # Sample non-CFR to prevent rendering issues (keep all CFR)
    if len(non_cfr_afs) > 5000:
        rng = np.random.RandomState(42)
        non_cfr_sample = rng.choice(non_cfr_afs, 5000, replace=False).tolist()
    else:
        non_cfr_sample = non_cfr_afs

    df = pd.DataFrame([
        {'Allele Frequency': af, 'Position Type': 'CFR'} for af in cfr_afs
    ] + [
        {'Allele Frequency': af, 'Position Type': 'Non-CFR'} for af in non_cfr_sample
    ])

    # Log transform for better visualization
    df['Log10 Allele Frequency'] = np.log10(df['Allele Frequency'].clip(lower=1e-10))

    fig = px.violin(
        df, x='Position Type', y='Log10 Allele Frequency',
        color='Position Type',
        box=True, points='outliers',
        title='Variant Allele Frequency: CFR vs Non-CFR Positions',
        color_discrete_map={'CFR': '#008080', 'Non-CFR': '#E8820C'},
    )

    stat = freq_result.get('stats', {})
    if stat:
        fig.add_annotation(
            x=0.5, y=1.08, xref='paper', yref='paper',
            text=f"Mann-Whitney U p = {stat['p_value']:.2e} | "
                 f"CFR median AF: {stat['cfr_median_af']:.2e} | "
                 f"Non-CFR median AF: {stat['non_cfr_median_af']:.2e}",
            showarrow=False, font=dict(size=11),
        )

    fig.update_layout(height=500, showlegend=False)
    return fig


def pathogenicity_enrichment(store, cfr_positions):
    """AlphaMissense pathogenicity enrichment at CFRs (chi-squared test)."""
    cfr_counts = {'pathogenic': 0, 'benign': 0, 'ambiguous': 0}
    non_cfr_counts = {'pathogenic': 0, 'benign': 0, 'ambiguous': 0}

    for gid in store.gpcr_ids:
        var_df = store.variant_data.get(gid, pd.DataFrame())
        if var_df.empty or 'am_class' not in var_df.columns:
            continue

        cfr_pos = cfr_positions.get(gid, set())
        for _, row in var_df.iterrows():
            am_class = str(row.get('am_class', '')).lower().strip()
            if am_class not in ('pathogenic', 'benign', 'ambiguous'):
                continue

            protein_pos = int(row['protein_position']) if pd.notna(row.get('protein_position')) else None
            if protein_pos is None:
                continue

            if protein_pos in cfr_pos:
                cfr_counts[am_class] += 1
            else:
                non_cfr_counts[am_class] += 1

    # Chi-squared test
    stat_result = {}
    obs = np.array([
        [cfr_counts['pathogenic'], cfr_counts['benign']],
        [non_cfr_counts['pathogenic'], non_cfr_counts['benign']],
    ])
    if obs.min() >= 0 and obs.sum() > 0:
        try:
            chi2, p_value, dof, expected = stats.chi2_contingency(obs)
            stat_result = {
                'test': 'Chi-squared',
                'chi2': float(chi2),
                'p_value': float(p_value),
                'dof': int(dof),
            }
        except ValueError:
            pass

    cfr_total = sum(cfr_counts.values())
    non_cfr_total = sum(non_cfr_counts.values())

    return {
        'cfr_counts': cfr_counts,
        'non_cfr_counts': non_cfr_counts,
        'cfr_total': cfr_total,
        'non_cfr_total': non_cfr_total,
        'cfr_pathogenic_pct': cfr_counts['pathogenic'] / cfr_total * 100 if cfr_total > 0 else 0,
        'non_cfr_pathogenic_pct': non_cfr_counts['pathogenic'] / non_cfr_total * 100 if non_cfr_total > 0 else 0,
        'stats': stat_result,
    }


def make_pathogenicity_bar(path_result):
    """Stacked bar chart: pathogenicity class distribution at CFR vs non-CFR."""
    cfr = path_result['cfr_counts']
    non_cfr = path_result['non_cfr_counts']
    cfr_total = path_result['cfr_total']
    non_cfr_total = path_result['non_cfr_total']

    if cfr_total == 0 and non_cfr_total == 0:
        return go.Figure()

    categories = ['Pathogenic', 'Ambiguous', 'Benign']
    cfr_pcts = [cfr.get(c.lower(), 0) / cfr_total * 100 if cfr_total > 0 else 0 for c in categories]
    non_cfr_pcts = [non_cfr.get(c.lower(), 0) / non_cfr_total * 100 if non_cfr_total > 0 else 0 for c in categories]

    colors = {'Pathogenic': '#E8820C', 'Ambiguous': '#F0C27A', 'Benign': '#008080'}

    fig = go.Figure()
    for i, cat in enumerate(categories):
        fig.add_trace(go.Bar(
            name=cat, x=['CFR Positions', 'Non-CFR Positions'],
            y=[cfr_pcts[i], non_cfr_pcts[i]],
            marker_color=colors[cat],
            text=[f'{cfr_pcts[i]:.1f}%', f'{non_cfr_pcts[i]:.1f}%'],
            textposition='inside',
        ))

    stat = path_result.get('stats', {})
    title = 'AlphaMissense Pathogenicity: CFR vs Non-CFR Positions'
    if stat:
        title += f" (χ² p = {stat['p_value']:.2e})"

    fig.update_layout(
        barmode='stack', title=title,
        yaxis_title='Percentage of Variants',
        height=450,
    )
    return fig


def conservation_vs_delta(store):
    """Scatter plot data: conservation score vs delta_rrcs across all variants."""
    rows = []

    for gid in store.gpcr_ids:
        var_df = store.variant_data.get(gid, pd.DataFrame())
        delta_df = store.delta_data.get(gid, pd.DataFrame())
        if var_df.empty or delta_df.empty:
            continue
        if 'conservation' not in var_df.columns:
            continue

        # Build position -> max abs_delta map
        pos_delta = {}
        for _, row in delta_df.iterrows():
            for res_col in ['res1', 'res2']:
                pos = int(row[res_col])
                val = abs(row['delta_rrcs'])
                if pos not in pos_delta or val > pos_delta[pos]:
                    pos_delta[pos] = val

        for _, row in var_df.iterrows():
            conservation = row.get('conservation')
            if pd.isna(conservation):
                continue
            try:
                conservation = float(conservation)
            except (ValueError, TypeError):
                continue

            protein_pos = int(row['protein_position']) if pd.notna(row.get('protein_position')) else None
            if protein_pos is None:
                continue

            delta = pos_delta.get(protein_pos, 0)
            rows.append({
                'conservation': conservation,
                'max_abs_delta': delta,
                'gpcr_id': gid,
            })

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def make_conservation_scatter(conservation_data):
    """Scatter: conservation vs delta_rrcs."""
    if conservation_data.empty:
        return go.Figure()

    # Sample if too many points
    df = conservation_data
    if len(df) > 5000:
        df = df.sample(5000, random_state=42)

    fig = px.scatter(
        df, x='conservation', y='max_abs_delta',
        opacity=0.3,
        title='Conservation Score vs Conformational Change',
        labels={
            'conservation': 'Conservation Score',
            'max_abs_delta': 'Max |ΔRRCS| at Position',
        },
        color_discrete_sequence=['#008080'],
    )
    fig.update_traces(marker_size=4)
    fig.update_layout(height=500)
    return fig


def find_high_impact_variants(store, cfr_positions):
    """Table of high-impact variants at CFR positions."""
    rows = []

    for gid in store.gpcr_ids:
        var_df = store.variant_data.get(gid, pd.DataFrame())
        if var_df.empty:
            continue

        cfr_pos = cfr_positions.get(gid, set())
        if not cfr_pos:
            continue

        for _, row in var_df.iterrows():
            protein_pos = int(row['protein_position']) if pd.notna(row.get('protein_position')) else None
            if protein_pos is None or protein_pos not in cfr_pos:
                continue

            am_class = str(row.get('am_class', '')).lower().strip()
            if am_class != 'pathogenic':
                continue

            rows.append({
                'gpcr_id': gid,
                'uniprot_name': store.name_map.get(gid, gid),
                'position': protein_pos,
                'ref_aa': row.get('ref_aa', ''),
                'alt_aa': row.get('alt_aa', ''),
                'allele_frequency': row.get('af', ''),
                'am_score': row.get('am_score', ''),
                'am_class': row.get('am_class', ''),
                'conservation': row.get('conservation', ''),
            })

    df = pd.DataFrame(rows) if rows else pd.DataFrame()
    if not df.empty:
        df = df.sort_values('am_score', ascending=False).head(100)
    return df
