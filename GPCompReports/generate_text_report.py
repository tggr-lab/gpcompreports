#!/usr/bin/env python3
"""Generate a comprehensive plain-text statistics report for GPCompReports.

Produces a standalone .txt file with maximum data depth: top 100 tables,
full 283-GPCR rankings, ASCII tables, and every available statistic.
No plots — just text and tables.

Usage:
    python GPCompReports/generate_text_report.py
    python GPCompReports/generate_text_report.py --output /path/to/report.txt
"""

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

# Add parent to path for package imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from GPCompReports.analysis.data_loader import GPCRDataStore
from GPCompReports.analysis.cross_gpcr_analysis import compute_rankings, run_pca_analysis
from GPCompReports.analysis.tm_domain_analysis import (
    compute_domain_stats, compute_domain_interactions, compute_generic_number_variation,
)
from GPCompReports.analysis.cfr_analysis import identify_cfrs, build_cfr_network
from GPCompReports.analysis.variant_correlation import run_variant_analysis


# ═══════════════════════════════════════════════════════════════════════════════
# ASCII TABLE FORMATTER
# ═══════════════════════════════════════════════════════════════════════════════

def format_table(headers, rows, alignments=None, max_col_width=40):
    """Format data as an ASCII table with auto-sized columns.

    Args:
        headers: list of column header strings
        rows: list of lists/tuples of cell values
        alignments: list of '<' (left) or '>' (right) per column; defaults to '<'
        max_col_width: maximum column width before truncation
    """
    if not rows:
        return "(no data)\n"

    n_cols = len(headers)
    if alignments is None:
        alignments = ['<'] * n_cols

    # Convert all cells to strings and compute column widths
    str_rows = []
    for row in rows:
        str_row = []
        for i, cell in enumerate(row):
            s = _fmt_cell(cell)
            if len(s) > max_col_width:
                s = s[:max_col_width - 1] + '~'
            str_row.append(s)
        str_rows.append(str_row)

    col_widths = [len(h) for h in headers]
    for row in str_rows:
        for i, cell in enumerate(row):
            if i < n_cols:
                col_widths[i] = max(col_widths[i], len(cell))

    # Build format strings
    sep = '+' + '+'.join('-' * (w + 2) for w in col_widths) + '+'
    header_cells = []
    for i, h in enumerate(headers):
        header_cells.append(f' {h:<{col_widths[i]}} ')
    header_line = '|' + '|'.join(header_cells) + '|'

    lines = [sep, header_line, sep]
    for row in str_rows:
        cells = []
        for i in range(n_cols):
            val = row[i] if i < len(row) else ''
            if alignments[i] == '>':
                cells.append(f' {val:>{col_widths[i]}} ')
            else:
                cells.append(f' {val:<{col_widths[i]}} ')
        lines.append('|' + '|'.join(cells) + '|')
    lines.append(sep)
    return '\n'.join(lines) + '\n'


def _fmt_cell(val):
    """Format a cell value for display."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ''
    if isinstance(val, float):
        if abs(val) >= 1000:
            return f'{val:,.1f}'
        if abs(val) >= 1:
            return f'{val:.3f}'
        if abs(val) >= 0.001:
            return f'{val:.4f}'
        return f'{val:.2e}'
    return str(val)


# ═══════════════════════════════════════════════════════════════════════════════
# REPORT SECTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def section_header(num, title):
    """Return a formatted section header."""
    rule = '=' * 80
    return f"\n\n{rule}\n  {num}. {title}\n{rule}\n\n"


def write_header_and_toc(f, store):
    """Section 0: Title, date, and table of contents."""
    f.write('=' * 80 + '\n')
    f.write('  GPCompReports — Comprehensive Statistics Report\n')
    f.write('=' * 80 + '\n\n')
    f.write(f'  Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n')
    f.write(f'  GPCRs analyzed: {len(store.gpcr_ids)}\n')
    f.write(f'  Data source: {store.batch_dir}\n')
    f.write(f'  Metadata: {store.metadata_csv}\n\n')
    f.write('-' * 80 + '\n')
    f.write('  TABLE OF CONTENTS\n')
    f.write('-' * 80 + '\n\n')
    toc = [
        ' 1. Dataset Overview',
        ' 2. Top 100 GPCRs by Conformational Change',
        ' 3. Complete GPCR Rankings (All GPCRs)',
        ' 4. Transmembrane Domain Statistics',
        ' 5. Domain-Domain Interaction Matrix',
        ' 6. Position-Level Variation (Top 100 Generic Numbers)',
        ' 7. Core Functional Residues (Top 100)',
        ' 8. CFR Contact Network (Top 50 Pairs)',
        ' 9. High-Impact Pathogenic Variants',
        '10. Statistical Tests (Mann-Whitney U, Chi-Squared)',
        '11. PCA Variance Explained',
        '12. Analysis by Ligand Type & Receptor Family',
        '13. Interpretation Guide',
    ]
    for item in toc:
        f.write(f'  {item}\n')
    f.write('\n')


def write_dataset_overview(f, store, info_df):
    """Section 1: Key summary statistics."""
    f.write(section_header(1, 'DATASET OVERVIEW'))

    n_gpcrs = len(store.gpcr_ids)
    n_with_variants = sum(1 for gid in store.gpcr_ids if len(store.variant_data.get(gid, pd.DataFrame())) > 0)
    total_contacts = int(info_df['total_contacts'].sum())
    total_sig = int(info_df['significant_changes'].sum())
    total_variants = int(info_df['variants_found'].sum())

    kv = [
        ('Total GPCRs analyzed', str(n_gpcrs)),
        ('GPCRs with variant data', str(n_with_variants)),
        ('Total contacts across all GPCRs', f'{total_contacts:,}'),
        ('Total significant changes', f'{total_sig:,}'),
        ('Total variants catalogued', f'{total_variants:,}'),
        ('', ''),
        ('Sum |delta| — mean', f'{info_df["sum_abs_delta"].mean():.2f}'),
        ('Sum |delta| — median', f'{info_df["sum_abs_delta"].median():.2f}'),
        ('Sum |delta| — std', f'{info_df["sum_abs_delta"].std():.2f}'),
        ('Sum |delta| — min', f'{info_df["sum_abs_delta"].min():.2f}'),
        ('Sum |delta| — max', f'{info_df["sum_abs_delta"].max():.2f}'),
        ('', ''),
        ('Mean |delta| — mean', f'{info_df["mean_abs_delta"].mean():.4f}'),
        ('Mean |delta| — median', f'{info_df["mean_abs_delta"].median():.4f}'),
        ('Mean |delta| — std', f'{info_df["mean_abs_delta"].std():.4f}'),
        ('Mean |delta| — min', f'{info_df["mean_abs_delta"].min():.4f}'),
        ('Mean |delta| — max', f'{info_df["mean_abs_delta"].max():.4f}'),
        ('', ''),
        ('Max |delta| — mean', f'{info_df["max_abs_delta"].mean():.4f}'),
        ('Max |delta| — median', f'{info_df["max_abs_delta"].median():.4f}'),
        ('Max |delta| — std', f'{info_df["max_abs_delta"].std():.4f}'),
        ('Max |delta| — min', f'{info_df["max_abs_delta"].min():.4f}'),
        ('Max |delta| — max', f'{info_df["max_abs_delta"].max():.4f}'),
    ]
    for label, val in kv:
        if label:
            f.write(f'  {label:.<45s} {val}\n')
        else:
            f.write('\n')


def write_top100_gpcrs(f, info_df):
    """Section 2: Top 100 GPCRs by sum |delta|."""
    f.write(section_header(2, 'TOP 100 GPCRs BY CONFORMATIONAL CHANGE'))

    top100 = info_df.nlargest(100, 'sum_abs_delta').reset_index(drop=True)
    headers = ['Rank', 'GPCR', 'Gene', 'Family', 'Ligand Type',
               'Contacts', 'Sig.Chg', 'Sum|dR|', 'Mean|dR|']
    aligns = ['>', '<', '<', '<', '<', '>', '>', '>', '>']
    rows = []
    for i, (_, r) in enumerate(top100.iterrows(), 1):
        rows.append([
            i,
            r['uniprot_name'],
            r['gene_name'],
            str(r['receptor_family'])[:20],
            str(r['ligand_type'])[:15],
            int(r['total_contacts']),
            int(r['significant_changes']),
            f'{r["sum_abs_delta"]:.2f}',
            f'{r["mean_abs_delta"]:.4f}',
        ])
    f.write(format_table(headers, rows, aligns))


def write_complete_rankings(f, info_df):
    """Section 3: All GPCRs ranked by sum |delta|."""
    f.write(section_header(3, 'COMPLETE GPCR RANKINGS — ALL GPCRs'))

    ranked = info_df.sort_values('sum_abs_delta', ascending=False).reset_index(drop=True)
    headers = ['Rank', 'GPCR', 'Gene', 'Family', 'Sum|dR|', 'Mean|dR|',
               'Max|dR|', 'Contacts', 'Sig.Chg', 'Variants']
    aligns = ['>', '<', '<', '<', '>', '>', '>', '>', '>', '>']
    rows = []
    for i, (_, r) in enumerate(ranked.iterrows(), 1):
        rows.append([
            i,
            r['uniprot_name'],
            r['gene_name'],
            str(r['receptor_family'])[:20],
            f'{r["sum_abs_delta"]:.2f}',
            f'{r["mean_abs_delta"]:.4f}',
            f'{r["max_abs_delta"]:.4f}',
            int(r['total_contacts']),
            int(r['significant_changes']),
            int(r['variants_found']),
        ])
    f.write(format_table(headers, rows, aligns))
    f.write(f'\n  Total GPCRs: {len(ranked)}\n')


def write_tm_domain_stats(f, domain_stats):
    """Section 4: Per-segment statistics."""
    f.write(section_header(4, 'TRANSMEMBRANE DOMAIN STATISTICS'))

    ds = domain_stats.sort_values('mean_abs_delta', ascending=False)
    headers = ['Segment', 'Mean|dR|', 'Median|dR|', 'Std|dR|', 'Total|dR|', 'Contact Count']
    aligns = ['<', '>', '>', '>', '>', '>']
    rows = []
    for _, r in ds.iterrows():
        rows.append([
            r['segment'],
            f'{r["mean_abs_delta"]:.4f}',
            f'{r["median_abs_delta"]:.4f}',
            f'{r["std_abs_delta"]:.4f}',
            f'{r["total_abs_delta"]:.1f}',
            int(r['count']),
        ])
    f.write(format_table(headers, rows, aligns))


def write_interaction_matrix(f, interaction_matrix):
    """Section 5: Domain-domain interaction matrix + top 20 pairs."""
    f.write(section_header(5, 'DOMAIN-DOMAIN INTERACTION MATRIX'))

    segments = interaction_matrix.columns.tolist()
    # Abbreviated labels for the matrix grid
    f.write('  Mean |delta| for interactions between segments:\n\n')

    # Column header
    label_width = 5
    f.write(' ' * (label_width + 1))
    for s in segments:
        f.write(f'{s:>6s}')
    f.write('\n')
    f.write(' ' * (label_width + 1) + '-' * (6 * len(segments)) + '\n')

    for seg in segments:
        f.write(f'{seg:>{label_width}s}|')
        for col in segments:
            val = interaction_matrix.loc[seg, col]
            if val == 0:
                f.write('     .')
            else:
                f.write(f'{val:6.2f}')
        f.write('\n')

    # Top 20 interacting pairs
    f.write('\n\n  Top 20 Domain-Domain Interaction Pairs (by mean |delta|):\n\n')
    pairs = []
    for i, seg1 in enumerate(segments):
        for j, seg2 in enumerate(segments):
            if j > i:  # Upper triangle only
                val = interaction_matrix.loc[seg1, seg2]
                if val > 0:
                    pairs.append((seg1, seg2, val))
    pairs.sort(key=lambda x: -x[2])

    headers = ['Rank', 'Segment 1', 'Segment 2', 'Mean |dR|']
    aligns = ['>', '<', '<', '>']
    rows = []
    for k, (s1, s2, v) in enumerate(pairs[:20], 1):
        rows.append([k, s1, s2, f'{v:.4f}'])
    f.write(format_table(headers, rows, aligns))


def write_position_variation(f, cv_data):
    """Section 6: Top 100 generic numbers by mean |delta| + sub-tables."""
    f.write(section_header(6, 'POSITION-LEVEL VARIATION — TOP 100 GENERIC NUMBERS'))

    if cv_data.empty:
        f.write('  (no data available)\n')
        return

    top100 = cv_data.nlargest(100, 'mean_abs_delta').reset_index(drop=True)
    headers = ['Rank', 'Generic Number', 'Segment', 'Mean|dR|', 'Std|dR|', 'CV', 'N Obs']
    aligns = ['>', '<', '<', '>', '>', '>', '>']
    rows = []
    for i, (_, r) in enumerate(top100.iterrows(), 1):
        rows.append([
            i,
            r['generic_number'],
            r['segment'],
            f'{r["mean_abs_delta"]:.4f}',
            f'{r["std_abs_delta"]:.4f}',
            f'{r["cv"]:.3f}',
            int(r['n_observations']),
        ])
    f.write(format_table(headers, rows, aligns))

    # Sub-table: Most Consistently Functional (high delta, low CV)
    f.write('\n  Most Consistently Functional Positions (high |delta|, low CV):\n')
    f.write('  (Sorted by mean |delta| descending, filtered to CV < median CV)\n\n')
    median_cv = cv_data['cv'].median()
    consistent = cv_data[cv_data['cv'] < median_cv].nlargest(20, 'mean_abs_delta')
    headers2 = ['Rank', 'Generic Number', 'Segment', 'Mean|dR|', 'CV', 'N Obs']
    aligns2 = ['>', '<', '<', '>', '>', '>']
    rows2 = []
    for i, (_, r) in enumerate(consistent.iterrows(), 1):
        rows2.append([
            i, r['generic_number'], r['segment'],
            f'{r["mean_abs_delta"]:.4f}', f'{r["cv"]:.3f}', int(r['n_observations']),
        ])
    f.write(format_table(headers2, rows2, aligns2))

    # Sub-table: Most Variable (high CV, above-median delta)
    f.write('\n  Most Variable Positions (highest CV among above-median |delta|):\n\n')
    median_delta = cv_data['mean_abs_delta'].median()
    variable = cv_data[cv_data['mean_abs_delta'] > median_delta].nlargest(20, 'cv')
    rows3 = []
    for i, (_, r) in enumerate(variable.iterrows(), 1):
        rows3.append([
            i, r['generic_number'], r['segment'],
            f'{r["mean_abs_delta"]:.4f}', f'{r["cv"]:.3f}', int(r['n_observations']),
        ])
    f.write(format_table(headers2, rows3, aligns2))


def write_cfr_table(f, cfr_table):
    """Section 7: Top 100 Core Functional Residues."""
    f.write(section_header(7, 'CORE FUNCTIONAL RESIDUES — TOP 100'))

    if cfr_table.empty:
        f.write('  (no data available)\n')
        return

    top100 = cfr_table.head(100)
    headers = ['Rank', 'Generic Number', 'Segment', 'Freq (N)', 'Freq (%)',
               'Mean|dR|', 'Max|dR|', 'CFR Score']
    aligns = ['>', '<', '<', '>', '>', '>', '>', '>']
    rows = []
    for _, r in top100.iterrows():
        rows.append([
            int(r['rank']),
            r['generic_number'],
            r['segment'],
            int(r['frequency']),
            f'{r["frequency_pct"]:.1f}%',
            f'{r["mean_abs_delta"]:.4f}',
            f'{r["max_abs_delta"]:.4f}',
            f'{r["cfr_score"]:.4f}',
        ])
    f.write(format_table(headers, rows, aligns))
    f.write(f'\n  Total CFR positions identified (freq >= 3): {len(cfr_table)}\n')


def write_cfr_network(f, cfr_network, n_gpcrs_total):
    """Section 8: Top 50 CFR contact pairs."""
    f.write(section_header(8, 'CFR CONTACT NETWORK — TOP 50 PAIRS'))

    if cfr_network.empty:
        f.write('  (no data available)\n')
        return

    top50 = cfr_network.head(50)
    headers = ['Rank', 'CFR Position 1', 'CFR Position 2', 'N GPCRs', '% of GPCRs']
    aligns = ['>', '<', '<', '>', '>']
    rows = []
    for i, (_, r) in enumerate(top50.iterrows(), 1):
        pct = r['n_gpcrs'] / n_gpcrs_total * 100 if n_gpcrs_total > 0 else 0
        rows.append([
            i,
            r['cfr_1'],
            r['cfr_2'],
            int(r['n_gpcrs']),
            f'{pct:.1f}%',
        ])
    f.write(format_table(headers, rows, aligns))
    f.write(f'\n  Total CFR-CFR contact pairs found: {len(cfr_network)}\n')


def write_high_impact_variants(f, high_impact):
    """Section 9: High-impact pathogenic variants at CFR positions."""
    f.write(section_header(9, 'HIGH-IMPACT PATHOGENIC VARIANTS AT CFR POSITIONS'))

    if high_impact.empty:
        f.write('  (no high-impact variants found)\n')
        return

    headers = ['Rank', 'GPCR', 'Position', 'Ref', 'Alt', 'Allele Freq',
               'AM Score', 'AM Class', 'Conservation']
    aligns = ['>', '<', '>', '<', '<', '>', '>', '<', '>']
    rows = []
    for i, (_, r) in enumerate(high_impact.iterrows(), 1):
        af = r.get('allele_frequency', '')
        af_str = f'{float(af):.2e}' if af != '' and not pd.isna(af) else ''
        am = r.get('am_score', '')
        am_str = f'{float(am):.4f}' if am != '' and not pd.isna(am) else ''
        cons = r.get('conservation', '')
        cons_str = f'{float(cons):.4f}' if cons != '' and not pd.isna(cons) else ''
        rows.append([
            i,
            r['uniprot_name'],
            int(r['position']),
            r.get('ref_aa', ''),
            r.get('alt_aa', ''),
            af_str,
            am_str,
            r.get('am_class', ''),
            cons_str,
        ])
    f.write(format_table(headers, rows, aligns))
    f.write(f'\n  Total high-impact variants shown: {len(high_impact)}\n')


def write_statistical_tests(f, variant_results):
    """Section 10: Mann-Whitney U and Chi-Squared test results."""
    f.write(section_header(10, 'STATISTICAL TESTS'))

    # Mann-Whitney U
    freq_stats = variant_results.get('freq_comparison', {}).get('stats', {})
    f.write('  --- Mann-Whitney U Test: Allele Frequency at CFR vs Non-CFR ---\n\n')
    if freq_stats:
        kv = [
            ('Test', freq_stats.get('test', 'Mann-Whitney U')),
            ('U statistic', f'{freq_stats["u_statistic"]:.1f}'),
            ('p-value', f'{freq_stats["p_value"]:.2e}'),
            ('CFR median allele frequency', f'{freq_stats["cfr_median_af"]:.2e}'),
            ('Non-CFR median allele frequency', f'{freq_stats["non_cfr_median_af"]:.2e}'),
            ('CFR sample size (N)', f'{freq_stats["cfr_n"]:,}'),
            ('Non-CFR sample size (N)', f'{freq_stats["non_cfr_n"]:,}'),
        ]
        for label, val in kv:
            f.write(f'  {label:.<50s} {val}\n')
        p = freq_stats['p_value']
        if p < 0.001:
            interp = 'Highly significant (p < 0.001)'
        elif p < 0.05:
            interp = 'Significant (p < 0.05)'
        else:
            interp = 'Not significant (p >= 0.05)'
        cfr_med = freq_stats['cfr_median_af']
        non_cfr_med = freq_stats['non_cfr_median_af']
        if cfr_med < non_cfr_med:
            direction = 'CFR positions have LOWER allele frequency (more constrained)'
        else:
            direction = 'CFR positions have HIGHER allele frequency'
        f.write(f'\n  Interpretation: {interp}\n')
        f.write(f'  Direction: {direction}\n')
    else:
        f.write('  (insufficient data for Mann-Whitney U test)\n')

    # Chi-Squared
    path_result = variant_results.get('pathogenicity', {})
    path_stats = path_result.get('stats', {})
    f.write('\n\n  --- Chi-Squared Test: Pathogenicity Enrichment at CFRs ---\n\n')
    if path_stats:
        kv = [
            ('Test', path_stats.get('test', 'Chi-squared')),
            ('Chi-squared statistic', f'{path_stats["chi2"]:.2f}'),
            ('p-value', f'{path_stats["p_value"]:.2e}'),
            ('Degrees of freedom', str(path_stats['dof'])),
            ('CFR pathogenic %', f'{path_result.get("cfr_pathogenic_pct", 0):.1f}%'),
            ('Non-CFR pathogenic %', f'{path_result.get("non_cfr_pathogenic_pct", 0):.1f}%'),
        ]
        for label, val in kv:
            f.write(f'  {label:.<50s} {val}\n')
        p = path_stats['p_value']
        if p < 0.001:
            interp = 'Highly significant (p < 0.001)'
        elif p < 0.05:
            interp = 'Significant (p < 0.05)'
        else:
            interp = 'Not significant (p >= 0.05)'
        f.write(f'\n  Interpretation: {interp}\n')
    else:
        f.write('  (insufficient data for Chi-squared test)\n')

    # Contingency table
    cfr_c = path_result.get('cfr_counts', {})
    non_cfr_c = path_result.get('non_cfr_counts', {})
    if cfr_c or non_cfr_c:
        f.write('\n\n  Contingency Table:\n\n')
        headers = ['', 'Pathogenic', 'Ambiguous', 'Benign', 'Total']
        aligns = ['<', '>', '>', '>', '>']
        cfr_total = sum(cfr_c.get(k, 0) for k in ['pathogenic', 'ambiguous', 'benign'])
        non_cfr_total = sum(non_cfr_c.get(k, 0) for k in ['pathogenic', 'ambiguous', 'benign'])
        rows = [
            ['CFR', cfr_c.get('pathogenic', 0), cfr_c.get('ambiguous', 0),
             cfr_c.get('benign', 0), cfr_total],
            ['Non-CFR', non_cfr_c.get('pathogenic', 0), non_cfr_c.get('ambiguous', 0),
             non_cfr_c.get('benign', 0), non_cfr_total],
            ['Total', cfr_c.get('pathogenic', 0) + non_cfr_c.get('pathogenic', 0),
             cfr_c.get('ambiguous', 0) + non_cfr_c.get('ambiguous', 0),
             cfr_c.get('benign', 0) + non_cfr_c.get('benign', 0),
             cfr_total + non_cfr_total],
        ]
        f.write(format_table(headers, rows, aligns))


def write_pca_variance(f, pca_variance):
    """Section 11: PCA explained variance."""
    f.write(section_header(11, 'PCA VARIANCE EXPLAINED'))

    if not pca_variance:
        f.write('  (PCA not computed — insufficient data)\n')
        return

    f.write('  PCA on per-GPCR TM-domain feature vectors (14 features: 7 TM + H8 + 3 ICL + 3 ECL)\n\n')
    headers = ['Component', 'Variance Explained', 'Cumulative']
    aligns = ['<', '>', '>']
    cumulative = 0
    rows = []
    for i, var in enumerate(pca_variance, 1):
        cumulative += var
        rows.append([
            f'PC{i}',
            f'{var:.1%}',
            f'{cumulative:.1%}',
        ])
    f.write(format_table(headers, rows, aligns))


def write_ligand_family_analysis(f, info_df):
    """Section 12: Breakdown by ligand type and receptor family."""
    f.write(section_header(12, 'ANALYSIS BY LIGAND TYPE & RECEPTOR FAMILY'))

    # By Ligand Type
    f.write('  --- By Ligand Type ---\n\n')
    lt_groups = info_df.groupby('ligand_type')['sum_abs_delta']
    lt_stats = lt_groups.agg(['count', 'mean', 'median', 'std', 'min', 'max'])
    lt_stats = lt_stats.sort_values('mean', ascending=False)
    headers = ['Ligand Type', 'N', 'Mean Sum|dR|', 'Median', 'Std', 'Min', 'Max']
    aligns = ['<', '>', '>', '>', '>', '>', '>']
    rows = []
    for lt, r in lt_stats.iterrows():
        rows.append([
            str(lt)[:25],
            int(r['count']),
            f'{r["mean"]:.2f}',
            f'{r["median"]:.2f}',
            f'{r["std"]:.2f}' if not pd.isna(r['std']) else '',
            f'{r["min"]:.2f}',
            f'{r["max"]:.2f}',
        ])
    f.write(format_table(headers, rows, aligns))

    # By Receptor Family
    f.write('\n  --- By Receptor Family (families with >= 2 members) ---\n\n')
    fam_groups = info_df.groupby('receptor_family')
    fam_stats = fam_groups['sum_abs_delta'].agg(['count', 'mean', 'median'])
    fam_stats = fam_stats[fam_stats['count'] >= 2].sort_values('mean', ascending=False)

    # Find top GPCR per family
    top_per_family = {}
    for fam in fam_stats.index:
        subset = info_df[info_df['receptor_family'] == fam]
        top_row = subset.nlargest(1, 'sum_abs_delta').iloc[0]
        top_per_family[fam] = (top_row['uniprot_name'], top_row['sum_abs_delta'])

    headers2 = ['Receptor Family', 'N', 'Mean Sum|dR|', 'Median Sum|dR|', 'Top GPCR', 'Its Sum|dR|']
    aligns2 = ['<', '>', '>', '>', '<', '>']
    rows2 = []
    for fam, r in fam_stats.iterrows():
        top_name, top_val = top_per_family.get(fam, ('', 0))
        rows2.append([
            str(fam)[:30],
            int(r['count']),
            f'{r["mean"]:.2f}',
            f'{r["median"]:.2f}',
            str(top_name),
            f'{top_val:.2f}',
        ])
    f.write(format_table(headers2, rows2, aligns2))


def write_interpretation_guide(f):
    """Section 13: Guide to interpreting the data."""
    f.write(section_header(13, 'INTERPRETATION GUIDE'))
    f.write("""  Metric Definitions:
  --------------------

  delta_rrcs ............. RRCS(active) - RRCS(inactive) for each residue-residue contact.
                           Positive = contact strengthens upon activation.
                           Negative = contact weakens upon activation.

  abs_delta (|dR|) ...... Absolute value of delta_rrcs. Measures magnitude of change
                           regardless of direction.

  Sum |dR| .............. Sum of abs_delta across all contacts for a GPCR. Measures
                           total conformational rearrangement.

  Mean |dR| ............. Average abs_delta per contact. Measures intensity of change
                           independent of receptor size.

  Max |dR| .............. Largest single abs_delta. Identifies the most dramatically
                           changed contact.

  Significant change ..... A contact where |delta_rrcs| exceeds the threshold (typically
                           2 standard deviations above the mean for that GPCR).

  CFR Score .............. Composite of normalized frequency (how many GPCRs show this
                           position as significant) and normalized mean |delta|. Range 0-1.

  CV (Coeff. of Var.) ... std / mean of |delta| at a generic number position across GPCRs.
                           Low CV = consistently functional. High CV = variable.

  AM Score .............. AlphaMissense pathogenicity score (0-1). Higher = more likely
                           pathogenic.

  Conservation .......... Sequence conservation score from variant analysis.


  Key Biological Insights:
  ------------------------

  - DRY motif (3.50x50) as CFR #1 confirms the ionic lock mechanism as the most
    universally important activation switch across Class A GPCRs.

  - TM6 positions (6.48, 6.44) in top CFRs reflect the toggle switch / rotamer toggle
    mechanism critical for G-protein coupling.

  - ICL2/ICL3 showing highest mean conformational change is expected, as these
    intracellular loops directly mediate G-protein and arrestin binding.

  - Lower allele frequency at CFR positions (if significant) indicates stronger
    purifying selection at functionally critical residues.

  - Higher pathogenicity enrichment at CFRs confirms that mutations at structurally
    important activation positions are more likely to cause disease.
""")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive text statistics report')
    parser.add_argument('--output', type=str, default=None,
                        help='Output file path (default: GPCompReports/output/comprehensive_statistics_report.txt)')
    parser.add_argument('--batch-dir', type=str, default=None, help='Batch analysis directory')
    parser.add_argument('--metadata', type=str, default=None, help='Path to class_A_all.csv')
    args = parser.parse_args()

    output_path = Path(args.output) if args.output else (
        Path(__file__).resolve().parent / 'output' / 'comprehensive_statistics_report.txt'
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    start = time.time()
    print('=' * 60)
    print('GPCompReports — Comprehensive Text Report Generator')
    print('=' * 60)

    # Step 1: Load data
    print('\n[1/6] Loading data...')
    store = GPCRDataStore(batch_dir=args.batch_dir, metadata_csv=args.metadata)
    store.load_all()
    info_df = store.get_all_info_df()

    # Step 2: Cross-GPCR rankings
    print('[2/6] Computing rankings...')
    rankings = compute_rankings(info_df)
    info_df_ranked = rankings['info_df_with_ranks']

    # Step 3: TM domain analysis
    print('[3/6] Computing TM domain statistics...')
    domain_stats = compute_domain_stats(store)
    interaction_matrix = compute_domain_interactions(store)
    cv_data = compute_generic_number_variation(store)

    # Step 4: CFR analysis
    print('[4/6] Identifying Core Functional Residues...')
    cfr_table = identify_cfrs(store)
    cfr_network = build_cfr_network(store, cfr_table)

    # Step 5: Variant analysis
    print('[5/6] Running variant correlation analysis...')
    variant_results = run_variant_analysis(store, cfr_table)

    # Step 6: PCA (we only need variance_explained)
    print('[6/6] Running PCA analysis...')
    pca_result = run_pca_analysis(store)
    pca_variance = pca_result['variance_explained']

    # Write report
    elapsed_analysis = time.time() - start
    print(f'\nAnalysis complete in {elapsed_analysis:.1f}s. Writing report...')

    with open(output_path, 'w', encoding='utf-8') as f:
        write_header_and_toc(f, store)
        write_dataset_overview(f, store, info_df_ranked)
        write_top100_gpcrs(f, info_df_ranked)
        write_complete_rankings(f, info_df_ranked)
        write_tm_domain_stats(f, domain_stats)
        write_interaction_matrix(f, interaction_matrix)
        write_position_variation(f, cv_data)
        write_cfr_table(f, cfr_table)
        write_cfr_network(f, cfr_network, len(store.gpcr_ids))
        write_high_impact_variants(f, variant_results.get('high_impact_variants', pd.DataFrame()))
        write_statistical_tests(f, variant_results)
        write_pca_variance(f, pca_variance)
        write_ligand_family_analysis(f, info_df_ranked)
        write_interpretation_guide(f)

        # Footer
        total_time = time.time() - start
        f.write('\n' + '=' * 80 + '\n')
        f.write(f'  End of report. Generated in {total_time:.1f}s.\n')
        f.write('=' * 80 + '\n')

    total_time = time.time() - start
    print(f'\nReport written to: {output_path}')
    print(f'Total time: {total_time:.1f}s')


if __name__ == '__main__':
    main()
