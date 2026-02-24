"""
Simple Report Generator for Batch GPCR Analysis
Generates clean, professional HTML reports with interactive plots for individual GPCRs.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional
from datetime import datetime

# Try to import plotly for interactive plots
try:
    import plotly.express as px
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: plotly not available, plots will be disabled")


class SimpleReportGenerator:
    """Generate streamlined HTML reports for GPCR analysis."""

    def __init__(self, gpcr_name: str, receptor: str, species: str,
                 uniprot_id: Optional[str], delta_matrix: pd.DataFrame,
                 variants_df: pd.DataFrame, tm_domains: Optional[Dict],
                 gpcrdb_numbers: Dict, lab_logo_path: Optional[str] = None,
                 snake_plot_svg: Optional[str] = None,
                 snake_plot_json: Optional[str] = None):
        """
        Initialize report generator.

        Args:
            gpcr_name: Full GPCR name
            receptor: Receptor code
            species: Species name
            uniprot_id: UniProt identifier
            delta_matrix: RRCS delta matrix
            variants_df: Variants data
            tm_domains: TM domain definitions
            gpcrdb_numbers: GPCRdb numbering
            lab_logo_path: Optional path to lab logo image file
            snake_plot_svg: Optional base SVG string for snake plot
            snake_plot_json: Optional JSON string with view data for snake plot
        """
        self.gpcr_name = gpcr_name
        self.receptor = receptor
        self.species = species
        self.uniprot_id = uniprot_id
        self.delta_matrix = delta_matrix
        self.variants_df = variants_df
        self.tm_domains = tm_domains
        self.gpcrdb_numbers = gpcrdb_numbers
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.lab_logo_base64 = self._load_lab_logo(lab_logo_path) if lab_logo_path else None
        self.snake_plot_svg = snake_plot_svg
        self.snake_plot_json = snake_plot_json

        # Calculate statistically-based significance threshold
        self.significance_threshold = self._calculate_significance_threshold()

    def _load_lab_logo(self, logo_path: str) -> Optional[str]:
        """Load and encode lab logo as base64."""
        try:
            import base64
            if Path(logo_path).exists():
                with open(logo_path, 'rb') as f:
                    logo_data = base64.b64encode(f.read()).decode('utf-8')
                    print(f"Loaded lab logo from {logo_path}")
                    return logo_data
            else:
                print(f"Lab logo file not found: {logo_path}")
                return None
        except Exception as e:
            print(f"Warning: Could not load lab logo: {e}")
            return None
    
    def _calculate_significance_threshold(self) -> float:
        """
        Calculate significance threshold based on statistical distribution.
        Uses mean + 2*std to identify outliers (top ~2.5% of changes).
        """
        abs_deltas = self.delta_matrix['delta_rrcs'].abs()
        mean_abs_delta = abs_deltas.mean()
        std_abs_delta = abs_deltas.std()

        # Mean + 2*std captures ~95% of data, beyond that is significant
        threshold = mean_abs_delta + 2 * std_abs_delta

        # Ensure minimum threshold of 1.0
        threshold = max(threshold, 1.0)

        print(f"Calculated significance threshold: {threshold:.2f} (mean: {mean_abs_delta:.2f}, std: {std_abs_delta:.2f})")
        return threshold

    def _generate_tm_rrcs_summary(self):
        """Generate TM domain-wise RRCS summary (active vs inactive preferences)."""
        if not self.tm_domains:
            return []

        tm_summary = []
        threshold = 3.0  # Standard significance threshold

        # Process each TM domain
        for tm_name, tm_range in self.tm_domains.items():
            start, end = tm_range

            # Find all residues in this TM domain with significant RRCS changes
            active_residues = []  # Positive ΔRRCS (active-favoring)
            inactive_residues = []  # Negative ΔRRCS (inactive-favoring)

            for position in range(start, end + 1):
                # Get all contacts involving this position
                position_changes = self.delta_matrix[
                    ((self.delta_matrix['res1'] == position) | (self.delta_matrix['res2'] == position))
                ]

                if position_changes.empty:
                    continue

                # Get max absolute ΔRRCS for this position
                max_change = position_changes.loc[position_changes['delta_rrcs'].abs().idxmax()]
                delta = float(max_change['delta_rrcs'])

                # Categorize if significant
                if abs(delta) >= threshold:
                    if delta > 0:
                        active_residues.append({
                            'position': position,
                            'delta_rrcs': delta
                        })
                    else:
                        inactive_residues.append({
                            'position': position,
                            'delta_rrcs': delta
                        })

            # Sort residues by magnitude of change
            active_residues.sort(key=lambda x: x['delta_rrcs'], reverse=True)
            inactive_residues.sort(key=lambda x: x['delta_rrcs'])  # Most negative first

            tm_summary.append({
                'domain': tm_name,
                'range': f"{start}-{end}",
                'active_residues': active_residues,
                'inactive_residues': inactive_residues,
                'total_active': len(active_residues),
                'total_inactive': len(inactive_residues)
            })

        return tm_summary

    def generate(self, output_path: Path):
        """
        Generate the HTML report.

        Args:
            output_path: Path to save the report
        """
        # Generate plots if plotly is available
        if PLOTLY_AVAILABLE:
            self._generate_plots()

        html_content = self._build_html()

        with open(output_path, 'w') as f:
            f.write(html_content)

    def _generate_plots(self):
        """Generate interactive plotly plots."""
        self.plots = {}

        try:
            # 1. RRCS Distribution
            fig_dist = px.histogram(
                self.delta_matrix,
                x='delta_rrcs',
                nbins=50,
                title='Distribution of ΔRRCS Values',
                labels={'delta_rrcs': 'ΔRRCS', 'count': 'Frequency'},
                template='plotly_white'
            )
            fig_dist.update_layout(height=400, font=dict(size=12))
            self.plots['distribution'] = fig_dist.to_html(full_html=False, include_plotlyjs='cdn')
        except Exception as e:
            print(f"Could not generate distribution plot: {e}")
            self.plots['distribution'] = None

        try:
            # 2. Active vs Inactive Scatter
            # Add formatted hover text with amino acid names
            self.delta_matrix['hover_text'] = self.delta_matrix.apply(
                lambda row: f"Residue 1: {row['res1_name']}{row['res1']}<br>" +
                           f"Residue 2: {row['res2_name']}{row['res2']}<br>" +
                           f"Active: {row['active_rrcs']:.2f}<br>" +
                           f"Inactive: {row['inactive_rrcs']:.2f}<br>" +
                           f"ΔRRCS: {row['delta_rrcs']:.2f}",
                axis=1
            )

            fig_scatter = px.scatter(
                self.delta_matrix,
                x='inactive_rrcs',
                y='active_rrcs',
                color='delta_rrcs',
                title='RRCS Comparison: Active vs Inactive State',
                labels={'active_rrcs': 'Active RRCS', 'inactive_rrcs': 'Inactive RRCS'},
                color_continuous_scale='RdBu_r',
                template='plotly_white',
                hover_name='hover_text'
            )
            fig_scatter.update_traces(hovertemplate='%{hovertext}<extra></extra>')
            # Add diagonal line
            min_val = min(self.delta_matrix['active_rrcs'].min(), self.delta_matrix['inactive_rrcs'].min())
            max_val = max(self.delta_matrix['active_rrcs'].max(), self.delta_matrix['inactive_rrcs'].max())
            fig_scatter.add_shape(type='line', x0=min_val, y0=min_val, x1=max_val, y1=max_val,
                                line=dict(color='gray', dash='dash', width=1))
            fig_scatter.update_layout(height=500, font=dict(size=12))
            self.plots['scatter'] = fig_scatter.to_html(full_html=False, include_plotlyjs='cdn')
        except Exception as e:
            print(f"Could not generate scatter plot: {e}")
            self.plots['scatter'] = None

        try:
            # 3. Residue-wise Changes with detailed info
            residue_changes = []
            all_residues = set(self.delta_matrix['res1']) | set(self.delta_matrix['res2'])

            for res in sorted(all_residues):
                changes = self.delta_matrix[
                    (self.delta_matrix['res1'] == res) | (self.delta_matrix['res2'] == res)
                ]['delta_rrcs']

                if not changes.empty:
                    # Get amino acid name
                    aa_name = self.delta_matrix[
                        (self.delta_matrix['res1'] == res)
                    ]['res1_name'].iloc[0] if len(self.delta_matrix[self.delta_matrix['res1'] == res]) > 0 else \
                             self.delta_matrix[self.delta_matrix['res2'] == res]['res2_name'].iloc[0]

                    # Get GPCRdb info if available
                    gpcrdb_info = ""
                    tm_info = ""
                    if self.gpcrdb_numbers and res in self.gpcrdb_numbers:
                        generic = self.gpcrdb_numbers[res].get('display_generic_number', '')
                        segment = self.gpcrdb_numbers[res].get('protein_segment', '')
                        if generic:
                            gpcrdb_info = f" ({generic})"
                        if segment:
                            tm_info = f" [{segment}]"

                    # Check TM domain
                    if not tm_info and self.tm_domains:
                        for tm_name, (start, end) in self.tm_domains.items():
                            if start <= res <= end:
                                tm_info = f" [{tm_name}]"
                                break

                    hover_text = f"{aa_name}{res}{gpcrdb_info}{tm_info}<br>Mean Δ: {changes.mean():.2f}<br>Max |Δ|: {changes.abs().max():.2f}"

                    residue_changes.append({
                        'residue': int(res),
                        'aa_name': aa_name,
                        'mean_change': float(changes.mean()),
                        'max_change': float(changes.abs().max()),
                        'hover_text': hover_text
                    })

            df_changes = pd.DataFrame(residue_changes)
            fig_res = px.scatter(
                df_changes, x='residue', y='mean_change', size='max_change',
                title='Residue-wise RRCS Changes',
                labels={'residue': 'Residue Position', 'mean_change': 'Mean ΔRRCS'},
                color='mean_change', color_continuous_scale='RdBu_r',
                template='plotly_white',
                hover_name='hover_text'
            )
            fig_res.update_traces(hovertemplate='%{hovertext}<extra></extra>')
            fig_res.update_layout(height=400, font=dict(size=12))
            self.plots['residues'] = fig_res.to_html(full_html=False, include_plotlyjs='cdn')
        except Exception as e:
            print(f"Could not generate residue plot: {e}")
            import traceback
            traceback.print_exc()
            self.plots['residues'] = None

        try:
            # 4. TM Domain Analysis
            if self.tm_domains:
                tm_changes = []
                for tm_name, (start, end) in self.tm_domains.items():
                    tm_data = self.delta_matrix[
                        ((self.delta_matrix['res1'] >= start) & (self.delta_matrix['res1'] <= end)) |
                        ((self.delta_matrix['res2'] >= start) & (self.delta_matrix['res2'] <= end))
                    ]
                    if not tm_data.empty:
                        tm_changes.append({
                            'domain': tm_name,
                            'mean_change': float(tm_data['delta_rrcs'].mean()),
                            'significant': int(len(tm_data[abs(tm_data['delta_rrcs']) >= self.significance_threshold]))
                        })

                df_tm = pd.DataFrame(tm_changes)
                fig_tm = px.bar(
                    df_tm, x='domain', y='mean_change',
                    title='TM Domain RRCS Changes',
                    labels={'domain': 'TM Domain', 'mean_change': 'Mean ΔRRCS'},
                    color='mean_change', color_continuous_scale='RdBu_r',
                    template='plotly_white'
                )
                fig_tm.update_layout(height=400, font=dict(size=12))
                self.plots['tm_domains'] = fig_tm.to_html(full_html=False, include_plotlyjs='cdn')
            else:
                self.plots['tm_domains'] = None
        except Exception as e:
            print(f"Could not generate TM plot: {e}")
            self.plots['tm_domains'] = None

    def _build_html(self) -> str:
        """Build the complete HTML report."""
        return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{self.receptor} - GPCR Conformational Analysis</title>
    {self._get_styles()}
</head>
<body>
    <div class="container">
        {self._get_header()}
        {self._get_summary_section()}
        {self._get_methods_section()}
        {self._get_top_changes_section()}
        {self._get_tm_rrcs_summary_section()}
        {self._get_snake_plot_section()}
        {self._get_tm_domain_section()}
        {self._get_plots_section() if PLOTLY_AVAILABLE and hasattr(self, 'plots') else ''}
        {self._get_variants_section()}
        {self._get_complete_rrcs_table()}
        {self._get_rrcs_section()}
        {self._get_interpretation_guide_section()}
        {self._get_footer()}
    </div>
    {self._get_table_scripts()}
</body>
</html>
        """
    
    def _get_styles(self) -> str:
        """Get CSS styles."""
        return """
<style>
    * {
        margin: 0;
        padding: 0;
        box-sizing: border-box;
    }
    
    body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background-color: #f5f7fa;
        color: #2c3e50;
        line-height: 1.6;
    }
    
    .container {
        max-width: 1200px;
        margin: 0 auto;
        padding: 20px;
    }
    
    .header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 40px;
        border-radius: 12px;
        margin-bottom: 30px;
        box-shadow: 0 8px 16px rgba(0,0,0,0.15);
    }

    .header h1 {
        font-size: 2.8em;
        margin-bottom: 10px;
        font-weight: 400;
        letter-spacing: -0.5px;
    }

    .header .subtitle {
        font-size: 1.3em;
        opacity: 0.95;
        margin-bottom: 15px;
    }

    .header .meta {
        margin-top: 20px;
        font-size: 0.95em;
        opacity: 0.9;
    }
    
    .section {
        background: white;
        padding: 35px;
        margin-bottom: 25px;
        border-radius: 12px;
        box-shadow: 0 4px 12px rgba(0,0,0,0.08);
        border: 1px solid #f0f0f0;
    }

    .section h2 {
        color: #667eea;
        margin-bottom: 25px;
        padding-bottom: 12px;
        border-bottom: 3px solid #e8eaf0;
        font-size: 2em;
        font-weight: 500;
        letter-spacing: -0.5px;
    }

    .section h3 {
        color: #4a5568;
        margin-top: 30px;
        margin-bottom: 15px;
        font-size: 1.4em;
        font-weight: 500;
    }
    
    .stats-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
        gap: 20px;
        margin: 20px 0;
    }
    
    .stat-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 25px;
        border-radius: 10px;
        text-align: center;
        border: none;
        box-shadow: 0 4px 8px rgba(102, 126, 234, 0.25);
        transition: transform 0.2s, box-shadow 0.2s;
    }

    .stat-card:hover {
        transform: translateY(-3px);
        box-shadow: 0 6px 12px rgba(102, 126, 234, 0.35);
    }

    .stat-value {
        font-size: 3em;
        font-weight: 600;
        color: white;
        display: block;
        margin-bottom: 8px;
    }

    .stat-label {
        color: rgba(255, 255, 255, 0.95);
        font-size: 0.9em;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        font-weight: 500;
    }
    
    table {
        width: 100%;
        border-collapse: collapse;
        margin-top: 20px;
        background: white;
        border-radius: 8px;
        overflow: hidden;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }

    th {
        padding: 12px;
        text-align: left;
        border-bottom: 1px solid #e8eaf0;
        background-color: #f7fafc;
        color: #4a5568;
        font-weight: 600;
        text-transform: uppercase;
        font-size: 0.85em;
        letter-spacing: 0.05em;
        cursor: pointer;
        user-select: none;
        position: relative;
    }

    th:hover {
        background-color: #e2e8f0;
    }

    th::after {
        content: ' ⇅';
        opacity: 0.5;
        font-size: 0.8em;
    }

    th.asc::after {
        content: ' ↑';
        opacity: 1;
    }

    th.desc::after {
        content: ' ↓';
        opacity: 1;
    }

    td {
        padding: 12px;
        text-align: left;
        border-bottom: 1px solid #e8eaf0;
    }

    tbody tr {
        transition: background-color 0.2s;
    }

    tr:hover {
        background-color: #f0f4ff !important;
    }

    .significant-row {
        background-color: #fff3e0 !important;
        font-weight: 500;
    }

    .positive {
        color: #48bb78;
        font-weight: 600;
    }

    .negative {
        color: #f56565;
        font-weight: 600;
    }

    .positive-delta {
        color: #d32f2f;
        font-weight: 600;
    }

    .negative-delta {
        color: #1976d2;
        font-weight: 600;
    }

    .controls {
        display: flex;
        gap: 15px;
        margin: 20px 0;
        align-items: center;
        flex-wrap: wrap;
    }

    .search-box {
        padding: 10px 15px;
        border: 2px solid #e0e0e0;
        border-radius: 5px;
        font-size: 1em;
        flex: 1;
        min-width: 200px;
    }

    .search-box:focus {
        outline: none;
        border-color: #667eea;
    }

    .btn {
        padding: 10px 20px;
        background-color: #667eea;
        color: white;
        border: none;
        border-radius: 5px;
        cursor: pointer;
        font-size: 1em;
        transition: background-color 0.3s;
    }

    .btn:hover {
        background-color: #764ba2;
    }

    .btn[disabled] {
        background-color: #cbd5e0;
        cursor: not-allowed;
        opacity: 0.6;
    }

    .pagination-controls {
        display: flex;
        align-items: center;
        gap: 15px;
        margin: 20px 0;
        justify-content: center;
        flex-wrap: wrap;
    }

    .page-info {
        font-weight: 500;
        color: #4a5568;
    }

    .page-size-select {
        padding: 8px 12px;
        border: 2px solid #e0e0e0;
        border-radius: 5px;
        font-size: 0.95em;
        cursor: pointer;
    }

    .page-size-select:focus {
        outline: none;
        border-color: #667eea;
    }

    .pagination-hidden {
        display: none !important;
    }
    
    .tm-domains {
        display: flex;
        flex-wrap: wrap;
        gap: 15px;
        margin: 20px 0;
    }
    
    .tm-domain {
        background: #f7fafc;
        padding: 15px 20px;
        border-radius: 8px;
        border-left: 4px solid #667eea;
    }
    
    .tm-domain strong {
        color: #667eea;
        font-size: 1.1em;
    }
    
    .warning-box {
        background: #fffbeb;
        border-left: 4px solid #f59e0b;
        padding: 15px;
        margin: 20px 0;
        border-radius: 4px;
    }
    
    .info-box {
        background: linear-gradient(to right, #eff6ff, #f9fafb);
        border-left: 4px solid #667eea;
        padding: 20px;
        margin: 20px 0;
        border-radius: 6px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }

    .info-box strong {
        color: #667eea;
    }
    
    .footer {
        text-align: center;
        color: #a0aec0;
        padding: 30px;
        font-size: 0.9em;
    }
    
    .badge {
        display: inline-block;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.85em;
        font-weight: 600;
    }
    
    .badge-high {
        background: #fed7d7;
        color: #c53030;
    }
    
    .badge-medium {
        background: #feebc8;
        color: #c05621;
    }
    
    .badge-low {
        background: #c6f6d5;
        color: #2f855a;
    }
</style>
        """
    
    def _get_header(self) -> str:
        """Generate enhanced header section with metadata and external links."""
        # Build external links
        links_html = ""

        if self.uniprot_id:
            links_html += f"""
            <a href="https://www.uniprot.org/uniprotkb/{self.uniprot_id}" target="_blank" class="db-link">
                <strong>UniProt</strong> {self.uniprot_id}
            </a>
            """

            # Use the full GPCR name for GPCRdb link (e.g., "5ht1a_human")
            gpcrdb_id = self.gpcr_name.lower().replace(' ', '_').replace('-', '_')
            links_html += f"""
            <a href="https://gpcrdb.org/protein/{gpcrdb_id}/" target="_blank" class="db-link">
                <strong>GPCRdb</strong>
            </a>
            """

            # Add gnomAD gene link
            links_html += f"""
            <a href="https://www.ncbi.nlm.nih.gov/gene/?term={self.uniprot_id}[Accession]" target="_blank" class="db-link">
                <strong>NCBI Gene</strong>
            </a>
            """

        return f"""
<style>
    .db-link {{
        display: inline-block;
        background-color: rgba(255,255,255,0.2);
        padding: 6px 12px;
        border-radius: 4px;
        color: white;
        text-decoration: none;
        margin: 5px;
        font-size: 0.9em;
        transition: background-color 0.3s;
    }}
    .db-link:hover {{
        background-color: rgba(255,255,255,0.3);
    }}
    .external-links {{
        margin-top: 15px;
        padding-top: 15px;
        border-top: 1px solid rgba(255,255,255,0.3);
    }}
</style>
<div class="header">
    <h1>{self.receptor}</h1>
    <div class="subtitle">{self.species} · Conformational State Analysis</div>
    <div class="meta">
        GPCR: {self.gpcr_name} ·
        {'UniProt: ' + self.uniprot_id if self.uniprot_id else 'UniProt: Not Available'} ·
        Generated: {self.timestamp}
    </div>
    {f'<div class="external-links">{links_html}</div>' if links_html else ''}
</div>
        """
    
    def _get_gpcr_metadata_section(self) -> str:
        """Generate GPCR metadata section with detailed information."""
        # Build metadata cards
        metadata_cards = ""

        if self.uniprot_id:
            metadata_cards += f"""
        <div class="metadata-card">
            <div class="metadata-label">UniProt ID</div>
            <div class="metadata-value">
                <a href="https://www.uniprot.org/uniprotkb/{self.uniprot_id}" target="_blank">{self.uniprot_id}</a>
            </div>
        </div>
            """

        # Add TM domains info
        if self.tm_domains:
            tm_count = len(self.tm_domains)
            metadata_cards += f"""
        <div class="metadata-card">
            <div class="metadata-label">TM Domains</div>
            <div class="metadata-value">{tm_count} domains</div>
        </div>
            """

        # Add species
        metadata_cards += f"""
        <div class="metadata-card">
            <div class="metadata-label">Species</div>
            <div class="metadata-value">{self.species}</div>
        </div>
        """

        # Add receptor type
        metadata_cards += f"""
        <div class="metadata-card">
            <div class="metadata-label">Receptor Type</div>
            <div class="metadata-value">GPCR</div>
        </div>
        """

        metadata_css = """
    .metadata-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
        gap: 15px;
        margin: 20px 0;
    }
    .metadata-card {
        background: #f7fafc;
        padding: 15px;
        border-radius: 6px;
        border-left: 3px solid #667eea;
    }
    .metadata-label {
        font-size: 0.85em;
        color: #718096;
        margin-bottom: 5px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    .metadata-value {
        font-size: 1.1em;
        font-weight: 600;
        color: #2d3748;
    }
    .metadata-value a {
        color: #667eea;
        text-decoration: none;
    }
    .metadata-value a:hover {
        text-decoration: underline;
    }
        """

        return f"""
<style>{metadata_css}</style>
<div class="section">
    <h2>GPCR Information</h2>
    <div class="metadata-grid">
        {metadata_cards}
    </div>
</div>
        """

    def _get_summary_section(self) -> str:
        """Generate summary statistics section."""
        total_contacts = len(self.delta_matrix)
        significant = len(self.delta_matrix[abs(self.delta_matrix['delta_rrcs']) >= self.significance_threshold])
        max_increase = self.delta_matrix['delta_rrcs'].max()
        max_decrease = self.delta_matrix['delta_rrcs'].min()

        return f"""
<div class="section">
    <h2>Analysis Summary</h2>
    <div class="stats-grid">
        <div class="stat-card">
            <span class="stat-value">{total_contacts}</span>
            <div class="stat-label">Total Contacts</div>
        </div>
        <div class="stat-card">
            <span class="stat-value">{significant}</span>
            <div class="stat-label">Significant Changes</div>
        </div>
        <div class="stat-card">
            <span class="stat-value">{max_increase:.2f}</span>
            <div class="stat-label">Max Increase</div>
        </div>
        <div class="stat-card">
            <span class="stat-value">{max_decrease:.2f}</span>
            <div class="stat-label">Max Decrease</div>
        </div>
    </div>

    <div class="info-box">
        <strong>Analysis Method:</strong> This analysis compares residue-residue contact scores (RRCS)
        between active and inactive conformational states of the receptor. Changes with |ΔRRCS| ≥ {self.significance_threshold:.1f} are considered
        significant and may indicate functionally important conformational rearrangements.
    </div>
</div>
        """
    
    def _get_plots_section(self) -> str:
        """Generate interactive plots section."""
        if not hasattr(self, 'plots'):
            return ""

        plots_html = '<div class="section"><h2>Interactive Visualizations</h2>'

        if self.plots.get('distribution'):
            plots_html += f'<h3 style="margin-top: 20px;">RRCS Distribution</h3><div class="plot-container">{self.plots["distribution"]}</div>'

        if self.plots.get('scatter'):
            plots_html += f'<h3 style="margin-top: 30px;">Active vs Inactive Comparison</h3><div class="plot-container">{self.plots["scatter"]}</div>'

        if self.plots.get('residues'):
            plots_html += f'<h3 style="margin-top: 30px;">Residue-wise Changes</h3><div class="plot-container">{self.plots["residues"]}</div>'

        if self.plots.get('tm_domains'):
            plots_html += f'<h3 style="margin-top: 30px;">TM Domain Analysis</h3><div class="plot-container">{self.plots["tm_domains"]}</div>'

        plots_html += '</div>'
        return plots_html

    def _get_rrcs_section(self) -> str:
        """Generate RRCS statistics section."""
        # Distribution statistics
        delta_values = self.delta_matrix['delta_rrcs']
        mean_delta = delta_values.mean()
        std_delta = delta_values.std()
        median_delta = delta_values.median()
        
        # Count by magnitude
        high_changes = len(self.delta_matrix[abs(self.delta_matrix['delta_rrcs']) >= 5.0])
        med_changes = len(self.delta_matrix[(abs(self.delta_matrix['delta_rrcs']) >= self.significance_threshold) & 
                                            (abs(self.delta_matrix['delta_rrcs']) < 5.0)])
        
        return f"""
<div class="section">
    <h2>RRCS Change Distribution</h2>
    
    <div class="stats-grid">
        <div class="stat-card">
            <span class="stat-value">{mean_delta:.2f}</span>
            <div class="stat-label">Mean Δ RRCS</div>
        </div>
        <div class="stat-card">
            <span class="stat-value">{std_delta:.2f}</span>
            <div class="stat-label">Std Dev</div>
        </div>
        <div class="stat-card">
            <span class="stat-value">{median_delta:.2f}</span>
            <div class="stat-label">Median</div>
        </div>
    </div>
    
    <h3>Change Magnitude Classification</h3>
    <div class="stats-grid">
        <div class="stat-card">
            <span class="stat-value">{high_changes}</span>
            <div class="stat-label">High (|Δ| ≥ 5.0)</div>
        </div>
        <div class="stat-card">
            <span class="stat-value">{med_changes}</span>
            <div class="stat-label">Medium (self.significance_threshold ≤ |Δ| < 5.0)</div>
        </div>
    </div>
</div>
        """
    
    def _get_snake_plot_section(self) -> str:
        """Generate interactive snake plot section with multi-view recoloring."""
        if not self.snake_plot_svg or not self.snake_plot_json:
            return ""

        # Escape </script> in JSON to prevent premature tag closure
        safe_json = self.snake_plot_json.replace('</script>', '<\\/script>')

        return f"""
<style>
    .snake-toolbar {{
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        margin-bottom: 15px;
        align-items: center;
    }}
    .snake-view-btn {{
        padding: 8px 16px;
        background: #f7fafc;
        border: 2px solid #e0e0e0;
        border-radius: 6px;
        cursor: pointer;
        font-size: 0.9em;
        transition: all 0.2s;
        color: #4a5568;
        font-family: inherit;
    }}
    .snake-view-btn:hover {{
        border-color: #667eea;
        color: #667eea;
    }}
    .snake-view-btn.active {{
        background: #667eea;
        color: white;
        border-color: #667eea;
    }}
    .snake-links-toggle {{
        margin-left: auto;
        display: flex;
        align-items: center;
        gap: 6px;
        font-size: 0.9em;
        color: #4a5568;
        cursor: pointer;
        user-select: none;
    }}
    .snake-links-toggle input {{
        width: 16px;
        height: 16px;
        cursor: pointer;
    }}
    #snake-legend {{
        margin: 15px 0;
        padding: 12px 16px;
        background: #f7fafc;
        border-radius: 6px;
        border: 1px solid #e8eaf0;
        min-height: 40px;
    }}
    .snake-legend-title {{
        font-weight: 600;
        color: #2d3748;
        margin-bottom: 4px;
    }}
    .snake-legend-desc {{
        font-size: 0.85em;
        color: #718096;
        margin-bottom: 8px;
    }}
    .snake-legend-gradient {{
        height: 16px;
        border-radius: 4px;
        border: 1px solid #e0e0e0;
    }}
    .snake-legend-labels {{
        display: flex;
        justify-content: space-between;
        font-size: 0.8em;
        color: #718096;
        margin-top: 4px;
    }}
    .snake-legend-swatches {{
        display: flex;
        flex-wrap: wrap;
        gap: 12px;
    }}
    .snake-legend-swatch {{
        display: flex;
        align-items: center;
        gap: 6px;
        font-size: 0.85em;
        color: #4a5568;
    }}
    .swatch-color {{
        display: inline-block;
        width: 14px;
        height: 14px;
        border-radius: 3px;
        border: 1px solid #ccc;
    }}
    /* Hide expanded loop segments by default (GPCRdb toggle pattern) */
    #snake-plot-container .long {{
        display: none;
    }}
    #snake-plot-container .segment {{
        cursor: pointer;
    }}
    #snake-plot-container svg {{
        width: 100%;
        height: auto;
        max-width: 900px;
        margin: 0 auto;
        display: block;
    }}
    #snake-links path {{
        transition: opacity 0.15s;
    }}
    #snake-links path:hover {{
        opacity: 1 !important;
        stroke-width: 4px;
    }}
</style>

<div class="section">
    <h2>Snake Plot</h2>

    <div class="info-box">
        <p><strong>Interactive GPCRdb-style snake plot.</strong> Click the buttons below to switch between
        different data views. Each residue circle is colored according to the selected view.
        Toggle contact links to see significant residue-residue contacts drawn as arcs.
        Click loop labels (ICL/ECL) to expand or collapse loop regions.</p>
        <p style="margin-top: 8px;"><strong>Color convention:</strong>
        <span style="color: #2166AC; font-weight: 600;">Blue = active-favoring</span>,
        <span style="color: #B2182B; font-weight: 600;">Red = inactive-favoring</span></p>
    </div>

    <div class="snake-toolbar">
        <button class="snake-view-btn active" data-view="delta_rrcs" onclick="switchSnakeView('delta_rrcs')">Delta RRCS</button>
        <button class="snake-view-btn" data-view="active_rrcs" onclick="switchSnakeView('active_rrcs')">Active RRCS</button>
        <button class="snake-view-btn" data-view="inactive_rrcs" onclick="switchSnakeView('inactive_rrcs')">Inactive RRCS</button>
        <button class="snake-view-btn" data-view="variants" onclick="switchSnakeView('variants')">Variants</button>
        <button class="snake-view-btn" data-view="conservation" onclick="switchSnakeView('conservation')">Conservation</button>
        <button class="snake-view-btn" data-view="alphamissense" onclick="switchSnakeView('alphamissense')">AlphaMissense</button>
        <button class="snake-view-btn" data-view="cfr" onclick="switchSnakeView('cfr')">Core Functional</button>
        <button class="snake-view-btn" data-view="aa_properties" onclick="switchSnakeView('aa_properties')">AA Properties</button>
        <label class="snake-links-toggle">
            <input type="checkbox" id="snake-links-checkbox" onchange="toggleSnakeLinks(this.checked)">
            Show Contact Links
        </label>
    </div>

    <div id="snake-legend"></div>

    <div id="snake-plot-container">
        {self.snake_plot_svg}
    </div>
</div>

<script type="application/json" id="snake-plot-data">
{safe_json}
</script>

<script>
(function() {{
    var snakeData = JSON.parse(document.getElementById('snake-plot-data').textContent);
    var views = snakeData.views;
    var links = snakeData.contact_links;
    var positions = snakeData.positions;
    var linksDrawn = false;

    function textColorForFill(hex) {{
        try {{
            var r = parseInt(hex.substr(1, 2), 16);
            var g = parseInt(hex.substr(3, 2), 16);
            var b = parseInt(hex.substr(5, 2), 16);
            var brightness = 0.2126 * r + 0.7152 * g + 0.0722 * b;
            return brightness > 128 ? '#000000' : '#FFFFFF';
        }} catch (e) {{
            return '#000000';
        }}
    }}

    function switchSnakeView(viewName) {{
        var view = views[viewName];
        if (!view) return;
        var colors = view.colors;

        var container = document.getElementById('snake-plot-container');
        var circles = container.querySelectorAll('circle');

        circles.forEach(function(circle) {{
            var id = circle.getAttribute('id');
            if (!id || !/^\\d+$/.test(id)) return;

            var color = colors[id] || '#FFFFFF';
            circle.setAttribute('fill', color);

            var textEl = document.getElementById(id + 't');
            if (textEl) {{
                textEl.setAttribute('fill', textColorForFill(color));
            }}
        }});

        // Update active button
        var buttons = document.querySelectorAll('.snake-view-btn');
        buttons.forEach(function(btn) {{ btn.classList.remove('active'); }});
        var activeBtn = document.querySelector('[data-view="' + viewName + '"]');
        if (activeBtn) activeBtn.classList.add('active');

        // Update legend
        updateLegend(view);
    }}

    function updateLegend(view) {{
        var el = document.getElementById('snake-legend');
        var html = '<div class="snake-legend-title">' + view.label + '</div>';
        html += '<div class="snake-legend-desc">' + view.description + '</div>';

        if (view.legend_type === 'diverging' || view.legend_type === 'sequential') {{
            var cols = view.legend_colors;
            var labs = view.legend_labels;
            html += '<div class="snake-legend-gradient" style="background:linear-gradient(to right,' + cols.join(',') + ');"></div>';
            html += '<div class="snake-legend-labels">';
            for (var i = 0; i < labs.length; i++) {{
                html += '<span>' + labs[i] + '</span>';
            }}
            html += '</div>';
        }} else if (view.legend_type === 'categorical') {{
            var items = view.legend_items || [];
            html += '<div class="snake-legend-swatches">';
            for (var i = 0; i < items.length; i++) {{
                html += '<span class="snake-legend-swatch"><span class="swatch-color" style="background:' + items[i].color + ';"></span>' + items[i].label + '</span>';
            }}
            html += '</div>';
        }}

        if (Object.keys(view.colors).length === 0 && view.legend_type !== 'categorical') {{
            html += '<div style="color:#a0aec0; font-size:0.85em; margin-top:6px;">No data available for this view.</div>';
        }}

        el.innerHTML = html;
    }}

    function toggleSnakeLinks(show) {{
        var container = document.getElementById('snake-plot-container');
        // Append links inside <g id='snake'> so they share its transform
        var snakeGroup = container.querySelector('#snake');
        if (!snakeGroup) return;

        if (!linksDrawn && show) {{
            var ns = 'http://www.w3.org/2000/svg';
            var g = document.createElementNS(ns, 'g');
            g.setAttribute('id', 'snake-links');

            links.forEach(function(link) {{
                var p1 = positions[String(link.from_seq)];
                var p2 = positions[String(link.to_seq)];
                if (!p1 || !p2) return;

                var midX = (p1[0] + p2[0]) / 2;
                var midY = (p1[1] + p2[1]) / 2;
                var dist = Math.abs(link.to_seq - link.from_seq);
                var arcH = 30 + dist * 0.5;
                var direction = dist < 20 ? -1 : 1;
                var ctrlY = midY + direction * arcH;

                var path = document.createElementNS(ns, 'path');
                path.setAttribute('d', 'M ' + p1[0] + ',' + p1[1] + ' Q' + midX + ',' + ctrlY + ' ' + p2[0] + ',' + p2[1]);
                path.setAttribute('stroke', link.color);
                path.setAttribute('stroke-width', String(link.width));
                path.setAttribute('fill', 'none');
                path.setAttribute('opacity', '0.6');

                var title = document.createElementNS(ns, 'title');
                title.textContent = 'Res ' + link.from_seq + ' \\u2194 ' + link.to_seq + ' | \\u0394RRCS: ' + link.delta;
                path.appendChild(title);

                g.appendChild(path);
            }});

            snakeGroup.appendChild(g);
            linksDrawn = true;
        }}

        var linksGroup = document.getElementById('snake-links');
        if (linksGroup) {{
            linksGroup.style.display = show ? '' : 'none';
        }}
    }}

    // GPCRdb loop toggle: switch between collapsed (short) and expanded (long)
    function toggleLoop(selector, currentMode, expanded, element) {{
        var container = document.getElementById('snake-plot-container');
        var otherMode = currentMode === 'short' ? 'long' : 'short';

        // Hide elements in current mode
        container.querySelectorAll(selector + '.' + currentMode).forEach(function(el) {{
            el.style.display = 'none';
        }});

        // Show elements in other mode — must use explicit value to
        // override the CSS rule `.long {{ display: none }}`
        container.querySelectorAll(selector + '.' + otherMode).forEach(function(el) {{
            el.style.display = 'inline';
        }});
    }}

    // Expose globally for onclick handlers
    window.switchSnakeView = switchSnakeView;
    window.toggleSnakeLinks = toggleSnakeLinks;
    window.toggleLoop = toggleLoop;

    // Initialize with Delta RRCS view on load
    document.addEventListener('DOMContentLoaded', function() {{
        switchSnakeView('delta_rrcs');
    }});
}})();
</script>
"""

    def _get_tm_domain_section(self) -> str:
        """Generate basic TM domain section."""
        if not self.tm_domains:
            return """
<div class="section">
    <h2>Transmembrane Domains</h2>
    <div class="warning-box">
        <strong>No TM domain data available.</strong> UniProt ID mapping may be missing or GPCRdb data unavailable.
    </div>
</div>
            """

        tm_html = """
<div class="section">
    <h2>Transmembrane Domains</h2>
    <div class="tm-domains">
        """

        for tm_name, (start, end) in sorted(self.tm_domains.items()):
            tm_html += f"""
        <div class="tm-domain">
            <strong>{tm_name}</strong><br>
            Residues {start}-{end}
        </div>
            """

        tm_html += """
    </div>
</div>
        """

        return tm_html

    def _get_tm_rrcs_summary_section(self) -> str:
        """Generate TM Domain RRCS Summary table."""
        tm_summary = self._generate_tm_rrcs_summary()

        if not tm_summary:
            return """
<div class="section">
    <h2>Transmembrane Domain Analysis</h2>
    <div class="info-box">
        TM domain information not available for this analysis.
    </div>
</div>
            """

        # Additional CSS for active/inactive coloring
        tm_css = """
    .active-residues { color: #00796b; }
    .inactive-residues { color: #f57c00; }
        """

        table_rows = ""
        for domain in tm_summary:
            # Format active-favoring residues
            active_html = "—"
            if domain['active_residues']:
                active_parts = []
                for res in domain['active_residues']:
                    active_parts.append(f"<strong>{res['position']}</strong> ({res['delta_rrcs']:.1f})")
                active_html = ", ".join(active_parts)

            # Format inactive-favoring residues
            inactive_html = "—"
            if domain['inactive_residues']:
                inactive_parts = []
                for res in domain['inactive_residues']:
                    inactive_parts.append(f"<strong>{res['position']}</strong> ({res['delta_rrcs']:.1f})")
                inactive_html = ", ".join(inactive_parts)

            table_rows += f"""
        <tr>
            <td><strong>{domain['domain']}</strong></td>
            <td style="text-align: center;">{domain['range']}</td>
            <td class="active-residues">{active_html}</td>
            <td style="text-align: center;" class="active-residues"><strong>{domain['total_active']}</strong></td>
            <td class="inactive-residues">{inactive_html}</td>
            <td style="text-align: center;" class="inactive-residues"><strong>{domain['total_inactive']}</strong></td>
        </tr>
            """

        return f"""
<style>{tm_css}</style>
<div class="section">
    <h2>Transmembrane Domain Analysis</h2>

    <div class="info-box">
        <p><strong>What this table shows:</strong> Residues within each transmembrane (TM) domain that undergo
        significant conformational changes (|ΔRRCS| ≥ 3.0). Residues are categorized based on whether they form
        stronger contacts in the <strong style="color: #00796b;">active state</strong> (positive ΔRRCS) or
        <strong style="color: #f57c00;">inactive state</strong> (negative ΔRRCS).</p>

        <p style="margin-top: 10px;"><strong>How to read:</strong> Each row shows a TM domain with its residue range.
        Active-favoring residues strengthen contacts during activation, while inactive-favoring residues strengthen
        contacts in the resting state.</p>
    </div>

    <div class="controls">
        <button class="btn" onclick="exportTableToCSV('tm-summary-table', 'tm_domain_summary.csv')">
            📥 Export to CSV
        </button>
    </div>

    <table id="tm-summary-table">
        <thead>
            <tr>
                <th onclick="sortTable('tm-summary-table', 0)" title="Transmembrane domain identifier">TM Domain</th>
                <th onclick="sortTable('tm-summary-table', 1)" title="Residue range of this domain">Domain Range</th>
                <th onclick="sortTable('tm-summary-table', 2)" title="Residues with positive ΔRRCS (stronger in active state)">Active-Favoring Residues (ΔRRCS &gt; 0)</th>
                <th onclick="sortTable('tm-summary-table', 3)" title="Number of active-favoring residues">Count</th>
                <th onclick="sortTable('tm-summary-table', 4)" title="Residues with negative ΔRRCS (stronger in inactive state)">Inactive-Favoring Residues (ΔRRCS &lt; 0)</th>
                <th onclick="sortTable('tm-summary-table', 5)" title="Number of inactive-favoring residues">Count</th>
            </tr>
        </thead>
        <tbody>
            {table_rows}
        </tbody>
    </table>
</div>
        """
    
    def _get_variants_section(self) -> str:
        """Generate comprehensive variants section matching original GUI format."""
        if self.variants_df.empty:
            return """
<div class="section">
    <h2>Variants of Interest</h2>
    <div class="info-box">
        No variant data available. This may occur if population data could not be retrieved or
        no variants were found at positions with significant RRCS changes.
    </div>
</div>
            """

        # Sort by RRCS magnitude (descending) like original GUI
        df = self.variants_df.copy()

        # Add RRCS delta for each variant position
        if 'max_delta_rrcs' not in df.columns:
            df['max_delta_rrcs'] = df.apply(lambda row: self._get_max_delta_for_position(row.get('protein_position')), axis=1)

        df_sorted = df.nlargest(100, 'max_delta_rrcs') if 'max_delta_rrcs' in df.columns else df.head(100)

        # Add pathogenicity color coding CSS
        pathogenic_css = """
    .pathogenic { color: #d32f2f; font-weight: 500; }
    .benign { color: #388e3c; font-weight: 500; }
    .ambiguous { color: #f57c00; }
        """

        table_rows = ""
        for _, var in df_sorted.iterrows():
            # Extract data with fallbacks
            position = var.get('protein_position', var.get('position', 'N/A'))
            hgvsp = var.get('hgvsp', var.get('protein_change', 'N/A'))
            hgvsc = var.get('hgvsc', var.get('dna_change', 'N/A'))
            af = var.get('af', var.get('allele_frequency', 0))
            het_count = var.get('het_count', 0)
            hom_count = var.get('ac_hom', 0)
            delta_rrcs = var.get('max_delta_rrcs', 0)
            am_score = var.get('am_score')
            am_class = var.get('am_class', 'N/A')
            conservation = var.get('conservation')
            rsids = var.get('rsids', [])

            # Format values
            af_str = f"{af:.2e}" if af and af > 0 else "N/A"
            am_score_str = f"{am_score:.3f}" if am_score is not None else "—"
            conservation_str = f"{conservation:.2f}" if conservation is not None else "—"
            rsids_str = ", ".join(rsids) if rsids and isinstance(rsids, list) else (rsids if rsids else "—")

            # Pathogenicity class with color coding
            am_class_html = "—"
            if am_class and isinstance(am_class, str) and am_class != 'N/A':
                css_class = am_class.lower()
                am_class_html = f'<span class="{css_class}">{am_class}</span>'

            table_rows += f"""
        <tr>
            <td><strong>{position}</strong></td>
            <td>{hgvsp}</td>
            <td style="font-family: monospace; font-size: 0.9em;">{hgvsc}</td>
            <td>{af_str}</td>
            <td style="text-align: center;">{het_count} / {hom_count}</td>
            <td style="text-align: right;"><strong>{delta_rrcs:.2f}</strong></td>
            <td style="text-align: center;">{am_score_str}</td>
            <td>{am_class_html}</td>
            <td style="text-align: center;">{conservation_str}</td>
            <td style="font-size: 0.85em;">{rsids_str}</td>
        </tr>
            """

        return f"""
<style>{pathogenic_css}</style>
<div class="section">
    <h2>Variants of Interest</h2>

    <div class="info-box">
        <p><strong>What this table shows:</strong> Genetic variants (mutations) found at residues that undergo significant
        conformational changes between protein states. These variants may affect protein function by disrupting critical
        structural transitions.</p>

        <p style="margin-top: 15px;"><strong>Column Descriptions:</strong></p>
        <ul style="margin-left: 20px; margin-top: 5px;">
            <li><strong>Position:</strong> Residue position in the protein sequence</li>
            <li><strong>Protein Change:</strong> Amino acid substitution (e.g., p.Ala123Val)</li>
            <li><strong>DNA Change:</strong> Corresponding DNA-level mutation</li>
            <li><strong>Allele Frequency:</strong> How common the variant is in the population (lower = rarer)</li>
            <li><strong>Heterozygotes/Homozygotes:</strong> Number of individuals with one or two copies</li>
            <li><strong>ΔRRCS:</strong> Change in contact strength between states (higher magnitude = more impact)</li>
            <li><strong>AlphaMissense Score:</strong> AI prediction of deleteriousness (0 = benign, 1 = pathogenic)</li>
            <li><strong>Conservation Score:</strong> How preserved this position is across species (1 = highly conserved)</li>
            <li><strong>dbSNP ID:</strong> Reference SNP identifier(s) from dbSNP database</li>
        </ul>
    </div>

    <div class="controls">
        <input type="text" id="search-variants-table" class="search-box"
               placeholder="Search variants..." onkeyup="searchTable('variants-table')">
        <button class="btn" onclick="exportTableToCSV('variants-table', 'variants_of_interest.csv')">
            📥 Export to CSV
        </button>
    </div>

    <table id="variants-table">
        <thead>
            <tr>
                <th onclick="sortTable('variants-table', 0)" title="Residue position in protein sequence">Position</th>
                <th onclick="sortTable('variants-table', 1)" title="Amino acid substitution">Protein Change</th>
                <th onclick="sortTable('variants-table', 2)" title="DNA-level change">DNA Change</th>
                <th onclick="sortTable('variants-table', 3)" title="Minor Allele Frequency from gnomAD">Allele Frequency</th>
                <th onclick="sortTable('variants-table', 4)" title="Number of heterozygous / homozygous carriers">Heterozygotes / Homozygotes</th>
                <th onclick="sortTable('variants-table', 5)" title="Change in residue-residue contact score">ΔRRCS</th>
                <th onclick="sortTable('variants-table', 6)" title="AlphaMissense pathogenicity score (0-1)">AlphaMissense Score</th>
                <th onclick="sortTable('variants-table', 7)" title="AlphaMissense classification">Pathogenicity Class</th>
                <th onclick="sortTable('variants-table', 8)" title="Sequence conservation score (0-1)">Conservation</th>
                <th onclick="sortTable('variants-table', 9)" title="dbSNP reference SNP ID">dbSNP ID (rs#)</th>
            </tr>
        </thead>
        <tbody>
            {table_rows}
        </tbody>
    </table>

    <p style="margin-top: 15px; color: #718096; font-size: 0.9em;">
        Total variants found: {len(self.variants_df)}
    </p>
</div>
        """

    def _get_max_delta_for_position(self, position) -> float:
        """Get maximum absolute ΔRRCS for a given residue position."""
        if position is None:
            return 0.0
        try:
            position_changes = self.delta_matrix[
                (self.delta_matrix['res1'] == position) |
                (self.delta_matrix['res2'] == position)
            ]
            if not position_changes.empty:
                return float(position_changes['delta_rrcs'].abs().max())
        except:
            pass
        return 0.0
    
    def _get_top_changes_section(self) -> str:
        """Generate top changes section."""
        # Top 100 most significant changes
        top_changes = self.delta_matrix.nlargest(100, 'abs_delta')

        table_rows = ""
        for idx, (_, row) in enumerate(top_changes.iterrows(), 1):
            delta = row['delta_rrcs']
            delta_class = "positive-delta" if delta > 0 else "negative-delta"
            significant = "significant-row" if abs(delta) >= self.significance_threshold else ""

            # Magnitude badge
            abs_delta = abs(delta)
            if abs_delta >= 5.0:
                badge = '<span class="badge badge-high">HIGH</span>'
            elif abs_delta >= self.significance_threshold:
                badge = '<span class="badge badge-medium">MED</span>'
            else:
                badge = '<span class="badge badge-low">LOW</span>'

            # Get GPCRdb numbers if available
            res1_gpcrdb = self.gpcrdb_numbers.get(row['res1'], {}).get('display_generic_number', '')
            res2_gpcrdb = self.gpcrdb_numbers.get(row['res2'], {}).get('display_generic_number', '')

            res1_display = f"{row['res1']}{row['res1_name']}"
            res2_display = f"{row['res2']}{row['res2_name']}"
            if res1_gpcrdb:
                res1_display += f" <small>({res1_gpcrdb})</small>"
            if res2_gpcrdb:
                res2_display += f" <small>({res2_gpcrdb})</small>"

            table_rows += f"""
        <tr class="{significant}">
            <td>{idx}</td>
            <td>{res1_display}</td>
            <td>{res2_display}</td>
            <td>{row['active_rrcs']:.3f}</td>
            <td>{row['inactive_rrcs']:.3f}</td>
            <td class="{delta_class}"><strong>{delta:+.3f}</strong></td>
            <td>{badge}</td>
        </tr>
            """

        return f"""
<div class="section">
    <h2>Top 100 Conformational Changes</h2>
    <p>Residue pairs with the largest RRCS changes between active and inactive states. GPCRdb numbering in parentheses.</p>

    <div class="controls">
        <input type="text" id="search-top-changes-table" class="search-box"
               placeholder="Search residues..." onkeyup="searchTable('top-changes-table')">
        <button class="btn" onclick="exportTableToCSV('top-changes-table', 'top_changes.csv')">
            📥 Export to CSV
        </button>
    </div>

    <table id="top-changes-table">
        <thead>
            <tr>
                <th onclick="sortTable('top-changes-table', 0)">Rank</th>
                <th onclick="sortTable('top-changes-table', 1)">Residue 1</th>
                <th onclick="sortTable('top-changes-table', 2)">Residue 2</th>
                <th onclick="sortTable('top-changes-table', 3)">Active RRCS</th>
                <th onclick="sortTable('top-changes-table', 4)">Inactive RRCS</th>
                <th onclick="sortTable('top-changes-table', 5)">Δ RRCS</th>
                <th onclick="sortTable('top-changes-table', 6)">Magnitude</th>
            </tr>
        </thead>
        <tbody>
            {table_rows}
        </tbody>
    </table>

    <div id="pagination-top-changes-table"></div>
</div>
        """
    
    def _get_complete_rrcs_table(self) -> str:
        """Generate complete RRCS results table with all contacts."""
        # Show all data, pagination will handle display
        display_data = self.delta_matrix.nlargest(1000, 'abs_delta')

        rows = ""
        for _, row in display_data.iterrows():
            # Color code based on delta direction
            delta_class = "positive-delta" if row['delta_rrcs'] > 0 else "negative-delta"
            significant = "significant-row" if abs(row['delta_rrcs']) >= self.significance_threshold else ""

            rows += f"""
            <tr class="{significant}">
                <td><strong>{row['res1']}</strong></td>
                <td>{row['res1_name']}</td>
                <td><strong>{row['res2']}</strong></td>
                <td>{row['res2_name']}</td>
                <td>{row['active_rrcs']:.3f}</td>
                <td>{row['inactive_rrcs']:.3f}</td>
                <td class="{delta_class}"><strong>{row['delta_rrcs']:.3f}</strong></td>
                <td>{row['abs_delta']:.3f}</td>
            </tr>
            """

        return f"""
<div class="section">
    <h2>Complete RRCS Results (Top 1000 by |ΔRRCS|)</h2>

    <div class="info-box">
        <strong>Significance Threshold:</strong> |ΔRRCS| ≥ {self.significance_threshold:.2f}
        (calculated as mean + 2×std of absolute deltas)
        <br><strong>Highlighted rows</strong> indicate statistically significant changes.
        <br><strong>Color coding:</strong> <span class="positive-delta">Red = Increased in active</span>,
        <span class="negative-delta">Blue = Increased in inactive</span>
    </div>

    <div class="controls">
        <input type="text" id="search-rrcs-table" class="search-box"
               placeholder="Search residues, amino acids..." onkeyup="searchTable('rrcs-table')">
        <button class="btn" onclick="exportTableToCSV('rrcs-table', 'complete_rrcs_results.csv')">
            📥 Export to CSV
        </button>
    </div>

    <table id="rrcs-table" class="sortable-table">
        <thead>
            <tr>
                <th onclick="sortTable('rrcs-table', 0)" title="First residue position">Res1</th>
                <th onclick="sortTable('rrcs-table', 1)" title="First residue amino acid">AA1</th>
                <th onclick="sortTable('rrcs-table', 2)" title="Second residue position">Res2</th>
                <th onclick="sortTable('rrcs-table', 3)" title="Second residue amino acid">AA2</th>
                <th onclick="sortTable('rrcs-table', 4)" title="RRCS in active state">Active RRCS</th>
                <th onclick="sortTable('rrcs-table', 5)" title="RRCS in inactive state">Inactive RRCS</th>
                <th onclick="sortTable('rrcs-table', 6)" title="Change in RRCS (Active - Inactive)">ΔRRCS</th>
                <th onclick="sortTable('rrcs-table', 7)" title="Absolute value of ΔRRCS">|ΔRRCS|</th>
            </tr>
        </thead>
        <tbody>
            {rows}
        </tbody>
    </table>

    <div id="pagination-rrcs-table"></div>
</div>
        """

    def _get_methods_section(self) -> str:
        """Generate Methods section explaining the analysis."""
        return f"""
<div class="section">
    <h2>Methods</h2>
    <div class="info-box">
        <p><strong>RRCS Analysis:</strong> Residue-Residue Contact Scores were calculated for each conformational state
        based on inter-atomic distances. Changes in contact patterns (ΔRRCS) were computed by comparing states.</p>

        <p style="margin-top: 10px;"><strong>RRCS Methodology:</strong>
        <a href="https://github.com/zhaolabSHT/RRCS" target="_blank" style="color: #667eea;">github.com/zhaolabSHT/RRCS</a></p>

        <p style="margin-top: 10px;"><strong>Significance Threshold:</strong> |ΔRRCS| ≥ {self.significance_threshold:.1f}
        (calculated as mean + 2×std of absolute deltas)</p>

        <p style="margin-top: 10px;"><strong>Data Sources:</strong></p>
        <ul style="margin-left: 20px; margin-top: 5px;">
            <li>Structural data: AlphaFold multistate models in PDB format</li>
            <li>Variant data: gnomAD v4 (exomes + genomes)</li>
            <li>Pathogenicity predictions: AlphaMissense</li>
            <li>Conservation scores: ProtVar/UniProt</li>
            <li>TM domain definitions: GPCRdb</li>
        </ul>
    </div>
</div>
        """

    def _get_interpretation_guide_section(self) -> str:
        """Generate Understanding the Results section."""
        return """
<div class="section">
    <h2>Understanding the Results</h2>

    <h3 style="margin-top: 20px; margin-bottom: 15px; color: #667eea;">What is RRCS?</h3>
    <div class="info-box">
        <p><strong>Residue-Residue Contact Score (RRCS)</strong> measures how strongly two amino acids interact in a
        protein structure. When a protein changes shape (e.g., from inactive to active), these contact strengths change.</p>

        <p style="margin-top: 10px;"><strong>ΔRRCS (Delta RRCS)</strong> shows the difference in contact strength between states:</p>
        <ul style="margin-left: 20px; margin-top: 5px;">
            <li><strong>Positive ΔRRCS:</strong> Contact is stronger in active state</li>
            <li><strong>Negative ΔRRCS:</strong> Contact is stronger in inactive state</li>
            <li><strong>Large |ΔRRCS| (≥3.0):</strong> Significant structural rearrangement</li>
        </ul>

        <p style="margin-top: 10px;"><strong>Why it matters:</strong> Residues with large RRCS changes are critical for protein
        function. Mutations at these positions may disrupt the protein's ability to change shape properly.</p>
    </div>

    <h3 style="margin-top: 30px; margin-bottom: 15px; color: #667eea;">Pathogenicity Predictions</h3>
    <div class="info-box">
        <p><strong>AlphaMissense</strong> uses artificial intelligence to predict whether a mutation will harm protein function:</p>
        <ul style="margin-left: 20px; margin-top: 5px;">
            <li><strong style="color: #d32f2f;">Pathogenic (score ≥ 0.564):</strong> Mutation likely damages protein function</li>
            <li><strong style="color: #f57c00;">Ambiguous (score 0.34-0.563):</strong> Effect uncertain</li>
            <li><strong style="color: #388e3c;">Benign (score &lt; 0.34):</strong> Mutation likely tolerated</li>
        </ul>
    </div>

    <h3 style="margin-top: 30px; margin-bottom: 15px; color: #667eea;">Conservation & Population Data</h3>
    <div class="info-box">
        <p><strong>Conservation Score:</strong> Measures how well-preserved a position is across different species.
        Scores range from 0 (highly variable) to 1 (perfectly conserved). High conservation suggests the position is
        critical for protein function.</p>

        <p style="margin-top: 10px;"><strong>Allele Frequency:</strong> How often a variant appears in the population.
        Rare variants (low frequency) may be more likely to be harmful, as natural selection removes common deleterious mutations.</p>

        <p style="margin-top: 10px;"><strong>Heterozygotes/Homozygotes:</strong> Number of individuals with one copy
        (heterozygous) or two copies (homozygous) of the variant. Absence of homozygotes for a common variant may suggest lethality.</p>
    </div>
</div>
        """

    def _get_footer(self) -> str:
        """Generate footer with optional lab logo."""
        logo_html = ""
        if self.lab_logo_base64:
            logo_html = f"""
    <img src="data:image/png;base64,{self.lab_logo_base64}" alt="Lab Logo" style="max-width: 150px; margin-bottom: 15px;">
            """

        return f"""
<div class="footer">
    {logo_html}
    <p><strong>GPCR Batch Analysis Pipeline</strong></p>
    <p>RRCS methodology: Residue-Residue Contact Score analysis of AlphaFold multistate models</p>
    <p style="margin-top: 10px; font-size: 0.85em;">Report generated on {self.timestamp}</p>
</div>
        """

    def _get_table_scripts(self) -> str:
        """Get JavaScript for table sorting, searching, and pagination."""
        return """
<script>
// Pagination state for each table
const paginationState = {};

// Initialize pagination for a table
function initPagination(tableId, pageSize = 25) {
    if (!paginationState[tableId]) {
        paginationState[tableId] = {
            currentPage: 1,
            pageSize: pageSize,
            totalRows: 0
        };
    }
    updatePagination(tableId);
}

// Update pagination display
function updatePagination(tableId) {
    const table = document.getElementById(tableId);
    const tbody = table.getElementsByTagName('tbody')[0];
    const rows = Array.from(tbody.getElementsByTagName('tr'));

    // Filter visible rows (after search)
    const visibleRows = rows.filter(row => row.style.display !== 'none');

    const state = paginationState[tableId];
    state.totalRows = visibleRows.length;
    const totalPages = Math.ceil(state.totalRows / state.pageSize);

    // Hide all rows first
    rows.forEach(row => row.classList.add('pagination-hidden'));

    // Show only current page rows
    const start = (state.currentPage - 1) * state.pageSize;
    const end = start + state.pageSize;
    visibleRows.slice(start, end).forEach(row => row.classList.remove('pagination-hidden'));

    // Update pagination controls
    const paginationDiv = document.getElementById('pagination-' + tableId);
    if (paginationDiv && totalPages > 1) {
        const pageInfo = `Page ${state.currentPage} of ${totalPages} (${state.totalRows} rows)`;
        const prevDisabled = state.currentPage <= 1 ? 'disabled' : '';
        const nextDisabled = state.currentPage >= totalPages ? 'disabled' : '';

        paginationDiv.innerHTML = `
            <div class="pagination-controls">
                <button onclick="changePage('${tableId}', ${state.currentPage - 1})" ${prevDisabled} class="btn">← Prev</button>
                <span class="page-info">${pageInfo}</span>
                <button onclick="changePage('${tableId}', ${state.currentPage + 1})" ${nextDisabled} class="btn">Next →</button>
                <select onchange="changePageSize('${tableId}', this.value)" class="page-size-select">
                    <option value="25" ${state.pageSize === 25 ? 'selected' : ''}>25 per page</option>
                    <option value="50" ${state.pageSize === 50 ? 'selected' : ''}>50 per page</option>
                    <option value="100" ${state.pageSize === 100 ? 'selected' : ''}>100 per page</option>
                    <option value="99999" ${state.pageSize === 99999 ? 'selected' : ''}>All</option>
                </select>
            </div>
        `;
    } else if (paginationDiv) {
        paginationDiv.innerHTML = `<div class="pagination-controls"><span class="page-info">${state.totalRows} rows</span></div>`;
    }
}

// Change page
function changePage(tableId, newPage) {
    const state = paginationState[tableId];
    const totalPages = Math.ceil(state.totalRows / state.pageSize);

    if (newPage >= 1 && newPage <= totalPages) {
        state.currentPage = newPage;
        updatePagination(tableId);
    }
}

// Change page size
function changePageSize(tableId, newSize) {
    const state = paginationState[tableId];
    state.pageSize = parseInt(newSize);
    state.currentPage = 1;
    updatePagination(tableId);
}

// Table search functionality
function searchTable(tableId) {
    const input = document.getElementById('search-' + tableId);
    const filter = input.value.toUpperCase();
    const table = document.getElementById(tableId);
    const tr = table.getElementsByTagName('tr');

    for (let i = 1; i < tr.length; i++) {
        let visible = false;
        const td = tr[i].getElementsByTagName('td');

        for (let j = 0; j < td.length; j++) {
            const cell = td[j];
            if (cell) {
                const text = cell.textContent || cell.innerText;
                if (text.toUpperCase().indexOf(filter) > -1) {
                    visible = true;
                    break;
                }
            }
        }
        tr[i].style.display = visible ? '' : 'none';
    }

    // Reset to page 1 and update pagination after search
    if (paginationState[tableId]) {
        paginationState[tableId].currentPage = 1;
        updatePagination(tableId);
    }
}

// Table sorting functionality
function sortTable(tableId, columnIndex) {
    const table = document.getElementById(tableId);
    const th = table.getElementsByTagName('th')[columnIndex];
    const tbody = table.getElementsByTagName('tbody')[0];
    const rows = Array.from(tbody.getElementsByTagName('tr'));

    // Determine sort direction
    const isAsc = th.classList.contains('asc');

    // Clear all sorting indicators
    Array.from(table.getElementsByTagName('th')).forEach(header => {
        header.classList.remove('asc', 'desc');
    });

    // Set new sorting indicator
    th.classList.add(isAsc ? 'desc' : 'asc');

    // Sort rows
    rows.sort((a, b) => {
        const aValue = a.getElementsByTagName('td')[columnIndex].textContent.trim();
        const bValue = b.getElementsByTagName('td')[columnIndex].textContent.trim();

        // Try numeric comparison first
        const aNum = parseFloat(aValue.replace(/[^0-9.-]/g, ''));
        const bNum = parseFloat(bValue.replace(/[^0-9.-]/g, ''));

        if (!isNaN(aNum) && !isNaN(bNum)) {
            return isAsc ? bNum - aNum : aNum - bNum;
        } else {
            return isAsc ?
                bValue.localeCompare(aValue) :
                aValue.localeCompare(bValue);
        }
    });

    // Reappend sorted rows
    rows.forEach(row => tbody.appendChild(row));
}

// CSV export functionality
function exportTableToCSV(tableId, filename) {
    const table = document.getElementById(tableId);
    const rows = table.getElementsByTagName('tr');
    const csv = [];

    for (let i = 0; i < rows.length; i++) {
        // Skip hidden rows
        if (rows[i].style.display === 'none') continue;

        const row = [];
        const cols = rows[i].querySelectorAll('td, th');

        for (let j = 0; j < cols.length; j++) {
            let text = cols[j].innerText.replace(/"/g, '""');
            row.push('"' + text + '"');
        }
        csv.push(row.join(','));
    }

    const csvFile = new Blob([csv.join('\\n')], {type: 'text/csv'});
    const downloadLink = document.createElement('a');
    downloadLink.download = filename;
    downloadLink.href = window.URL.createObjectURL(csvFile);
    downloadLink.style.display = 'none';
    document.body.appendChild(downloadLink);
    downloadLink.click();
    document.body.removeChild(downloadLink);
}

// Initialize all tables on page load
document.addEventListener('DOMContentLoaded', function() {
    // Find all tables with pagination
    const tables = document.querySelectorAll('table[id]');
    tables.forEach(table => {
        const tableId = table.id;
        // Initialize pagination if pagination div exists
        const paginationDiv = document.getElementById('pagination-' + tableId);
        if (paginationDiv) {
            initPagination(tableId, 25);
        }
    });
});
</script>
        """
