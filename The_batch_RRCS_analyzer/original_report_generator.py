"""
GPCompaRe Report Generator
Generates comprehensive HTML reports for GPCR conformational state analysis.
"""

import pandas as pd
import numpy as np
from jinja2 import Template
import plotly.graph_objects as go
import plotly.express as px
from datetime import datetime
import os
import webbrowser
import requests
import json
import re


class ReportGenerator:
    """Generate academic-style HTML reports for RRCS analysis."""
    
    def __init__(self, analyzer):
        """
        Initialize the report generator.
        
        Args:
            analyzer: RRCSAnalyzer instance with completed analysis
        """
        self.analyzer = analyzer
        self.report_data = {}
        
    def collect_data(self):
        """Collect all necessary data for the report with error handling."""
        if not hasattr(self.analyzer, 'delta_matrix') or self.analyzer.delta_matrix is None:
            raise ValueError("No analysis results available. Please run analysis first.")
            
        # Basic information
        self.report_data['timestamp'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.report_data['uniprot_id'] = getattr(self.analyzer, 'uniprot_id', 'N/A')
        
        # Get gene name
        try:
            self.report_data['gene_name'] = self.analyzer.get_gene_from_uniprot(self.analyzer.uniprot_id) or 'Unknown'
        except Exception as e:
            print(f"Warning: Could not fetch gene name: {e}")
            self.report_data['gene_name'] = 'Unknown'
        
        # Get actual gnomAD dataset used
        self.report_data['gnomad_version'] = getattr(self.analyzer, 'gnomad_dataset_used', 'gnomAD (version detection failed)')
        
        # Calculate RRCS statistics
        delta_matrix = self.analyzer.delta_matrix
        threshold = 3.0  # Standard threshold for significant changes
        
        significant_changes = delta_matrix[abs(delta_matrix['delta_rrcs']) >= threshold]
        
        self.report_data['rrcs_stats'] = {
            'total_contacts': len(delta_matrix),
            'significant_changes': len(significant_changes),
            'max_increase': float(delta_matrix['delta_rrcs'].max()),
            'max_decrease': float(delta_matrix['delta_rrcs'].min()),
            'mean_change': float(delta_matrix['delta_rrcs'].mean()),
            'std_change': float(delta_matrix['delta_rrcs'].std()),
            'threshold': threshold
        }
        
        # TM domain information
        if hasattr(self.analyzer, 'tm_definitions') and self.analyzer.tm_definitions:
            self.report_data['tm_domains'] = self.analyzer.tm_definitions
        else:
            self.report_data['tm_domains'] = None
        
        # Collect variant data
        self.collect_variant_data()
        
        # Generate TM domain RRCS summary
        self.generate_tm_rrcs_summary()
        
        # Generate plots
        self.generate_plots()
        
    def collect_variant_data(self):
        """Identify variants of interest based on multiple criteria."""
        variants_of_interest = []
        
        # Amino acid code conversion
        aa_codes = {
            'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E',
            'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N',
            'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S',
            'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
        }
        
        # Get positions with significant RRCS changes
        threshold = self.report_data['rrcs_stats']['threshold']
        significant_positions = set()
        
        for _, row in self.analyzer.delta_matrix.iterrows():
            if abs(row['delta_rrcs']) >= threshold:
                significant_positions.add(int(row['res1']))
                significant_positions.add(int(row['res2']))
        
        print(f"Analyzing {len(significant_positions)} positions with significant RRCS changes...")
        
        # Get gnomAD data
        try:
            gene_name = self.report_data.get('gene_name')
            if gene_name and gene_name != 'Unknown':
                gnomad_data = self.analyzer.get_gnomad_data(gene_name)
            else:
                gnomad_data = []
        except Exception as e:
            print(f"Warning: Could not fetch gnomAD data: {e}")
            gnomad_data = []
        
        # Create protein position to variant mapping
        protein_variants = {}
        for variant in gnomad_data:
            match = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', variant.get('hgvsp', ''))
            if match:
                ref_aa, pos, alt_aa = match.groups()
                pos = int(pos)
                if pos not in protein_variants:
                    protein_variants[pos] = []
                protein_variants[pos].append(variant)
        
        # Process each significant position
        for position in sorted(significant_positions):
            # Get RRCS context
            position_changes = self.analyzer.delta_matrix[
                (self.analyzer.delta_matrix['res1'] == position) |
                (self.analyzer.delta_matrix['res2'] == position)
            ]
            
            if position_changes.empty:
                continue
                
            max_delta = float(position_changes['delta_rrcs'].abs().max())
            mean_delta = float(position_changes['delta_rrcs'].mean())
            
            # Get ProtVar data
            try:
                protvar_data = self.analyzer.get_population_data(position)
            except Exception as e:
                print(f"Warning: Could not fetch ProtVar data for position {position}: {e}")
                protvar_data = None
            
            # Get variants at this position
            position_variants = protein_variants.get(position, [])
            
            # Process each variant
            for variant in position_variants:
                match = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', variant.get('hgvsp', ''))
                if not match:
                    continue
                    
                ref_aa, pos, alt_aa = match.groups()
                alt_aa_one = aa_codes.get(alt_aa, alt_aa)
                
                # Get prediction scores
                am_score = None
                am_class = 'N/A'
                conservation = None
                
                if protvar_data and 'score_data' in protvar_data:
                    for score in protvar_data['score_data']:
                        if score.get('name') == 'AM' and score.get('mt') == alt_aa_one:
                            am_score = score.get('amPathogenicity')
                            am_class = score.get('amClass', 'N/A')
                        elif score.get('name') == 'CONSERV' and conservation is None:
                            conservation = score.get('score')
                
                # Calculate priority score (DISABLED - not used currently)
                # priority = self._calculate_priority(max_delta, variant.get('af', 0), am_score, conservation)
                
                # Compile variant information
                variant_info = {
                    'position': position,
                    'ref_aa': ref_aa,
                    'alt_aa': alt_aa,
                    'protein_change': variant.get('hgvsp', 'N/A'),
                    'dna_change': variant.get('hgvsc', 'N/A'),
                    'allele_frequency': float(variant.get('af', 0)),
                    'het_count': int(variant.get('ac', 0) - 2 * variant.get('ac_hom', 0)),
                    'hom_count': int(variant.get('ac_hom', 0)),
                    'max_delta_rrcs': max_delta,
                    'mean_delta_rrcs': mean_delta,
                    'am_score': float(am_score) if am_score is not None else None,
                    'am_class': am_class,
                    'conservation': float(conservation) if conservation is not None else None,
                    'rsids': variant.get('rsids', []),
                    # 'priority': priority,  # DISABLED
                    'population_freqs': variant.get('population_frequencies', {})
                }
                
                variants_of_interest.append(variant_info)
        
        # Sort by RRCS magnitude (descending)
        variants_of_interest.sort(key=lambda x: abs(x['max_delta_rrcs']), reverse=True)
        self.report_data['variants'] = variants_of_interest
        print(f"Identified {len(variants_of_interest)} variants of interest")
    
    def _calculate_priority(self, delta_rrcs, allele_freq, am_score, conservation):
        """Calculate priority score for variant (0-100). CURRENTLY DISABLED."""
        priority = 0
        
        # RRCS change contribution (0-40 points)
        if abs(delta_rrcs) >= 5.0:
            priority += 40
        elif abs(delta_rrcs) >= 3.0:
            priority += 30
        elif abs(delta_rrcs) >= 2.0:
            priority += 20
        else:
            priority += 10
        
        # Rarity contribution (0-25 points)
        if allele_freq < 0.0001:
            priority += 25
        elif allele_freq < 0.001:
            priority += 15
        elif allele_freq < 0.01:
            priority += 5
        
        # Pathogenicity contribution (0-25 points)
        if am_score is not None:
            if am_score >= 0.8:
                priority += 25
            elif am_score >= 0.6:
                priority += 15
            elif am_score >= 0.4:
                priority += 5
        
        # Conservation contribution (0-10 points)
        if conservation is not None:
            if conservation >= 0.9:
                priority += 10
            elif conservation >= 0.7:
                priority += 5
        
        return priority
    
    def generate_tm_rrcs_summary(self):
        """Generate TM domain-wise RRCS summary (active vs inactive preferences)."""
        tm_summary = []
        
        # Check if we have TM domain definitions
        if not self.report_data.get('tm_domains'):
            self.report_data['tm_summary'] = []
            print("No TM domain definitions available for summary")
            return
        
        delta_matrix = self.analyzer.delta_matrix
        threshold = 3.0  # Significant RRCS change threshold
        
        # Process each TM domain
        for tm_name, tm_range in self.report_data['tm_domains'].items():
            start, end = tm_range
            
            # Find all residues in this TM domain with significant RRCS changes
            active_residues = []  # Positive ΔRRCS (active-favoring)
            inactive_residues = []  # Negative ΔRRCS (inactive-favoring)
            
            for position in range(start, end + 1):
                # Get all contacts involving this position
                position_changes = delta_matrix[
                    ((delta_matrix['res1'] == position) | (delta_matrix['res2'] == position))
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
        
        self.report_data['tm_summary'] = tm_summary
        print(f"Generated TM domain summary for {len(tm_summary)} domains")
        
    def generate_plots(self):
        """Generate all plots for the report with error handling."""
        delta_matrix = self.analyzer.delta_matrix
        
        # 1. RRCS Distribution
        try:
            fig_dist = px.histogram(
                delta_matrix,
            x='delta_rrcs',
            nbins=50,
                title='Distribution of ΔRRCS Values',
            labels={'delta_rrcs': 'ΔRRCS', 'count': 'Frequency'},
                template='plotly_white'
            )
            fig_dist.update_layout(
                font=dict(family="Arial, sans-serif", size=12),
                showlegend=False,
                height=400
            )
            self.report_data['rrcs_dist_plot'] = fig_dist.to_html(full_html=False, include_plotlyjs='cdn')
        except Exception as e:
            print(f"Warning: Could not generate distribution plot: {e}")
            self.report_data['rrcs_dist_plot'] = None
        
        # 2. Residue-wise Changes
        try:
            residue_changes = []
            all_residues = set(delta_matrix['res1'].unique()) | set(delta_matrix['res2'].unique())
            
            for res in sorted(all_residues):
                changes = delta_matrix[
                    (delta_matrix['res1'] == res) | (delta_matrix['res2'] == res)
            ]['delta_rrcs']
                
                if not changes.empty:
                    residue_changes.append({
                        'residue': int(res),
                        'mean_change': float(changes.mean()),
                        'max_change': float(changes.abs().max())
            })
        
            df_changes = pd.DataFrame(residue_changes)
            fig_res = px.scatter(
                df_changes,
                x='residue',
            y='mean_change',
            size='max_change',
            title='Residue-wise RRCS Changes',
            labels={
                'residue': 'Residue Position',
                'mean_change': 'Mean ΔRRCS',
                    'max_change': 'Max |ΔRRCS|'
            },
            color='mean_change',
                color_continuous_scale='RdBu_r',
                template='plotly_white'
            )
            fig_res.update_layout(
                font=dict(family="Arial, sans-serif", size=12),
                height=400
            )
            self.report_data['residue_plot'] = fig_res.to_html(full_html=False, include_plotlyjs='cdn')
        except Exception as e:
            print(f"Warning: Could not generate residue plot: {e}")
            self.report_data['residue_plot'] = None
        
        # 3. TM Domain Analysis
        if self.report_data['tm_domains']:
            try:
                tm_changes = []
                for tm, (start, end) in self.report_data['tm_domains'].items():
                    tm_data = delta_matrix[
                        ((delta_matrix['res1'] >= start) & (delta_matrix['res1'] <= end)) |
                        ((delta_matrix['res2'] >= start) & (delta_matrix['res2'] <= end))
                    ]
                    if not tm_data.empty:
                        tm_changes.append({
                            'domain': tm,
                            'mean_change': float(tm_data['delta_rrcs'].mean()),
                            'significant_changes': int(len(tm_data[abs(tm_data['delta_rrcs']) >= 3.0]))
                        })
                
                df_tm = pd.DataFrame(tm_changes)
                fig_tm = px.bar(
                    df_tm,
                    x='domain',
                    y='mean_change',
                    title='TM Domain RRCS Changes',
                    labels={'domain': 'TM Domain', 'mean_change': 'Mean ΔRRCS'},
                    color='mean_change',
                    color_continuous_scale='RdBu_r',
                    template='plotly_white'
                )
                fig_tm.update_layout(
                    font=dict(family="Arial, sans-serif", size=12),
                    height=400
                )
                self.report_data['tm_plot'] = fig_tm.to_html(full_html=False, include_plotlyjs='cdn')
            except Exception as e:
                print(f"Warning: Could not generate TM plot: {e}")
                self.report_data['tm_plot'] = None
        else:
            self.report_data['tm_plot'] = None
        
        # 4. RRCS Scatter Plot
        try:
            fig_scatter = px.scatter(
                delta_matrix,
                x='rrcs_active',
                y='rrcs_inactive',
                color='delta_rrcs',
                title='RRCS Comparison: Active vs Inactive State',
                labels={
                    'rrcs_active': 'Active State RRCS',
                    'rrcs_inactive': 'Inactive State RRCS',
                    'delta_rrcs': 'ΔRRCS'
                },
                color_continuous_scale='RdBu_r',
                template='plotly_white'
            )
            
            # Add diagonal line
            min_val = min(delta_matrix['rrcs_active'].min(), delta_matrix['rrcs_inactive'].min())
            max_val = max(delta_matrix['rrcs_active'].max(), delta_matrix['rrcs_inactive'].max())
            fig_scatter.add_shape(
                type='line',
                x0=min_val, y0=min_val,
                x1=max_val, y1=max_val,
                line=dict(color='gray', dash='dash', width=1)
            )
            
            fig_scatter.update_layout(
                font=dict(family="Arial, sans-serif", size=12),
                height=500
            )
            self.report_data['scatter_plot'] = fig_scatter.to_html(full_html=False, include_plotlyjs='cdn')
        except Exception as e:
            print(f"Warning: Could not generate scatter plot: {e}")
            self.report_data['scatter_plot'] = None
    
    def generate_html_report(self, output_dir='.'):
        """
        Generate comprehensive HTML report.
        
        Args:
            output_dir: Directory to save the report
            
        Returns:
            str: Path to generated report
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        report_path = os.path.join(output_dir, f"gpcompare_report_{timestamp}.html")
        
        # Encode lab logo as base64
        try:
            import base64
            logo_path = 'lablogo.png'
            if os.path.exists(logo_path):
                with open(logo_path, 'rb') as f:
                    logo_data = base64.b64encode(f.read()).decode('utf-8')
                    self.report_data['lab_logo_base64'] = logo_data
            else:
                self.report_data['lab_logo_base64'] = ''
        except Exception as e:
            print(f"Warning: Could not load lab logo: {e}")
            self.report_data['lab_logo_base64'] = ''
        
        # HTML template with academic styling
        template = Template("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GPCompaRe Analysis Report</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #2c3e50;
            background-color: #ecf0f1;
            padding: 20px;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }
        
        .header {
            background: linear-gradient(135deg, #00897b 0%, #00695c 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
            font-weight: 300;
            letter-spacing: 2px;
        }
        
        .header .subtitle {
            font-size: 1.2em;
            opacity: 0.9;
            margin-bottom: 20px;
        }
        
        .metadata {
            display: flex;
            justify-content: center;
            gap: 30px;
            margin-top: 20px;
            flex-wrap: wrap;
        }
        
        .metadata-item {
            background-color: rgba(255,255,255,0.1);
            padding: 10px 20px;
            border-radius: 5px;
        }
        
        .content {
            padding: 40px;
        }
        
        .section {
            margin-bottom: 50px;
        }
        
        .section-title {
            font-size: 1.8em;
            color: #00695c;
            border-bottom: 3px solid #00897b;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }
        
        .info-box {
            background-color: #e0f2f1;
            border-left: 4px solid #00897b;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 4px;
        }
        
        .warning-box {
            background-color: #fff3e0;
            border-left: 4px solid #ff6f00;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 4px;
        }
        
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        
        .stat-card {
            background: linear-gradient(135deg, #00897b 0%, #00695c 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        
        .stat-card h3 {
            font-size: 0.9em;
            font-weight: 400;
            margin-bottom: 10px;
            opacity: 0.9;
        }
        
        .stat-card .value {
            font-size: 2em;
            font-weight: 600;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
        }
        
        thead {
            background-color: #00695c;
            color: white;
        }
        
        th {
            padding: 12px;
            text-align: left;
            font-weight: 500;
            cursor: pointer;
            user-select: none;
            position: relative;
        }
        
        th:hover {
            background-color: #00897b;
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
            padding: 10px 12px;
            border-bottom: 1px solid #e0e0e0;
        }
        
        tbody tr {
            transition: background-color 0.2s;
        }
        
        tbody tr:hover {
            background-color: #b2dfdb !important;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
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
            border-color: #00897b;
        }
        
        .btn {
            padding: 10px 20px;
            background-color: #ff6f00;
            color: white;
            border: none;
            border-radius: 5px;
            cursor: pointer;
            font-size: 1em;
            transition: background-color 0.3s;
        }
        
        .btn:hover {
            background-color: #e65100;
        }
        
        .plot-container {
            margin: 30px 0;
            padding: 20px;
            background-color: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }
        
        .no-data {
            text-align: center;
            padding: 40px;
            color: #999;
            font-style: italic;
        }
        
        .footer {
            background-color: #f5f5f5;
            padding: 20px 40px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #e0e0e0;
        }
        
        
        /* Pathogenicity color coding */
        .pathogenic {
            color: #d32f2f;
            font-weight: 500;
        }
        
        .benign {
            color: #388e3c;
            font-weight: 500;
        }
        
        .ambiguous {
            color: #f57c00;
        }
        
        /* Tooltips for table headers */
        th[title] {
            cursor: help;
        }
        
        th[title]:hover {
            position: relative;
        }
        
        .footer {
            background-color: #f5f5f5;
            padding: 30px 40px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #e0e0e0;
        }
        
        .footer-logo {
            max-width: 150px;
            margin-bottom: 15px;
        }
        
        @media print {
            body {
                background-color: white;
                padding: 0;
            }
            .container {
                box-shadow: none;
            }
            .btn {
                display: none;
            }
        }
    </style>
        </head>
        <body>
            <div class="container">
        <!-- Header -->
                <div class="header">
            <h1>GPCompaRe</h1>
            <div class="subtitle">GPCR Conformational State Comparison Analysis</div>
            <div class="metadata">
                <div class="metadata-item">
                    <strong>Protein:</strong> {{ gene_name }} ({{ uniprot_id }})
                </div>
                <div class="metadata-item">
                    <strong>Generated:</strong> {{ timestamp }}
                </div>
            </div>
                </div>
        
        <!-- Content -->
        <div class="content">
            <!-- Analysis Summary -->
            <div class="section">
                <h2 class="section-title">Analysis Summary</h2>

                <div class="info-box">
                    <strong>Analysis Overview:</strong> This report presents the results of Residue-Residue Contact Score (RRCS) 
                    analysis comparing different conformational states of {{ gene_name }}. The analysis identifies residue contacts 
                    that undergo significant changes between states, potentially indicating functionally important regions.
                </div>

                <div class="stats-grid">
                    <div class="stat-card">
                            <h3>Total Contacts</h3>
                        <div class="value">{{ rrcs_stats.total_contacts }}</div>
                        </div>
                    <div class="stat-card">
                            <h3>Significant Changes</h3>
                        <div class="value">{{ rrcs_stats.significant_changes }}</div>
                        </div>
                    <div class="stat-card">
                            <h3>Max Increase</h3>
                        <div class="value">{{ "%.2f"|format(rrcs_stats.max_increase) }}</div>
                        </div>
                    <div class="stat-card">
                            <h3>Max Decrease</h3>
                        <div class="value">{{ "%.2f"|format(rrcs_stats.max_decrease) }}</div>
                        </div>
                    <div class="stat-card">
                            <h3>Mean Change</h3>
                        <div class="value">{{ "%.2f"|format(rrcs_stats.mean_change) }}</div>
                        </div>
                    </div>
                </div>

            <!-- Methods -->
                <div class="section">
                <h2 class="section-title">Methods</h2>
                <div class="info-box">
                    <p><strong>RRCS Analysis:</strong> Residue-Residue Contact Scores were calculated for each conformational state 
                    based on inter-atomic distances. Changes in contact patterns (ΔRRCS) were computed by comparing states.</p>
                    
                    <p style="margin-top: 10px;"><strong>Significance Threshold:</strong> |ΔRRCS| ≥ {{ "%.1f"|format(rrcs_stats.threshold) }}</p>
                    
                    <p style="margin-top: 10px;"><strong>Data Sources:</strong></p>
                    <ul style="margin-left: 20px; margin-top: 5px;">
                        <li>Structural data: PDB format files</li>
                        <li>Variant data: {{ gnomad_version }}</li>
                        <li>Pathogenicity predictions: AlphaMissense</li>
                        <li>Conservation scores: ProtVar/UniProt</li>
                        <li>TM domain definitions: GPCRdb</li>
                    </ul>
                </div>
            </div>
            
            <!-- TM Domain RRCS Summary -->
            {% if tm_summary and tm_summary|length > 0 %}
            <div class="section">
                <h2 class="section-title">Transmembrane Domain Analysis</h2>
                
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
                        Export to CSV
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
                        {% for domain in tm_summary %}
                        <tr>
                            <td><strong>{{ domain.domain }}</strong></td>
                            <td style="text-align: center;">{{ domain.range }}</td>
                            <td style="color: #00796b;">
                                {% if domain.active_residues|length > 0 %}
                                    {% for res in domain.active_residues %}
                                        <span style="display: inline-block; margin: 2px;">
                                            <strong>{{ res.position }}</strong> ({{ "%.1f"|format(res.delta_rrcs) }}){% if not loop.last %},{% endif %}
                                        </span>
                                    {% endfor %}
                                {% else %}
                                    —
                                {% endif %}
                            </td>
                            <td style="text-align: center; color: #00796b;"><strong>{{ domain.total_active }}</strong></td>
                            <td style="color: #f57c00;">
                                {% if domain.inactive_residues|length > 0 %}
                                    {% for res in domain.inactive_residues %}
                                        <span style="display: inline-block; margin: 2px;">
                                            <strong>{{ res.position }}</strong> ({{ "%.1f"|format(res.delta_rrcs) }}){% if not loop.last %},{% endif %}
                                        </span>
                                    {% endfor %}
                                {% else %}
                                    —
                                {% endif %}
                            </td>
                            <td style="text-align: center; color: #f57c00;"><strong>{{ domain.total_inactive }}</strong></td>
                        </tr>
                        {% endfor %}
                    </tbody>
                    </table>
                </div>
            {% else %}
                <div class="section">
                <h2 class="section-title">Transmembrane Domain Analysis</h2>
                <div class="no-data">
                    TM domain information not available for this analysis.
                </div>
            </div>
            {% endif %}

            <!-- Variants of Interest -->
            {% if variants and variants|length > 0 %}
                <div class="section">
                <h2 class="section-title">Variants of Interest</h2>
                
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
                    <input type="text" id="search-variants" class="search-box" 
                           placeholder="Search variants..." onkeyup="searchTable('variants-table')">
                    <button class="btn" onclick="exportTableToCSV('variants-table', 'variants.csv')">
                        Export to CSV
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
                            <th onclick="sortTable('variants-table', 6)" title="AlphaMissense pathogenicity score (0-1, higher = more pathogenic)">AlphaMissense Score</th>
                            <th onclick="sortTable('variants-table', 7)" title="AlphaMissense classification">Pathogenicity Class</th>
                            <th onclick="sortTable('variants-table', 8)" title="Sequence conservation score (0-1, higher = more conserved)">Conservation Score</th>
                            <th onclick="sortTable('variants-table', 9)" title="dbSNP reference SNP ID">dbSNP ID (rs#)</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for variant in variants %}
                        <tr>
                            <td><strong>{{ variant.position }}</strong></td>
                            <td>{{ variant.protein_change }}</td>
                            <td style="font-family: monospace; font-size: 0.9em;">{{ variant.dna_change }}</td>
                            <td>{{ "%.2e"|format(variant.allele_frequency) if variant.allele_frequency else "N/A" }}</td>
                            <td style="text-align: center;">{{ variant.het_count }} / {{ variant.hom_count }}</td>
                            <td style="text-align: right;"><strong>{{ "%.2f"|format(variant.max_delta_rrcs) }}</strong></td>
                            <td style="text-align: center;">{{ "%.3f"|format(variant.am_score) if variant.am_score is not none else "—" }}</td>
                            <td class="{% if variant.am_class == 'PATHOGENIC' %}pathogenic{% elif variant.am_class == 'BENIGN' %}benign{% else %}ambiguous{% endif %}">
                                {{ variant.am_class if variant.am_class != 'N/A' else '—' }}
                            </td>
                            <td style="text-align: center;">{{ "%.2f"|format(variant.conservation) if variant.conservation is not none else "—" }}</td>
                            <td style="font-size: 0.85em;">{{ variant.rsids|join(", ") if variant.rsids else "—" }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
                </div>
            {% else %}
                <div class="section">
                <h2 class="section-title">Variants of Interest</h2>
                <div class="no-data">
                    No variant data available. This may occur if population data could not be retrieved or 
                    no variants were found at positions with significant RRCS changes.
                </div>
                </div>
                {% endif %}

            <!-- Visualizations -->
                <div class="section">
                <h2 class="section-title">Data Visualizations</h2>
                
                {% if rrcs_dist_plot %}
                <h3 style="margin-top: 30px; margin-bottom: 15px; color: #00695c;">RRCS Distribution</h3>
                <div class="plot-container">
                    {{ rrcs_dist_plot|safe }}
                </div>
                {% endif %}
                
                {% if scatter_plot %}
                <h3 style="margin-top: 30px; margin-bottom: 15px; color: #00695c;">State Comparison</h3>
                <div class="plot-container">
                    {{ scatter_plot|safe }}
                </div>
                {% endif %}

                {% if residue_plot %}
                <h3 style="margin-top: 30px; margin-bottom: 15px; color: #00695c;">Residue-wise Analysis</h3>
                <div class="plot-container">
                    {{ residue_plot|safe }}
                </div>
                {% endif %}
                
                {% if tm_plot %}
                <h3 style="margin-top: 30px; margin-bottom: 15px; color: #00695c;">Transmembrane Domain Analysis</h3>
                <div class="plot-container">
                    {{ tm_plot|safe }}
                </div>
                {% endif %}

                {% if not rrcs_dist_plot and not scatter_plot and not residue_plot and not tm_plot %}
                <div class="no-data">
                    No visualization data available.
                </div>
                {% endif %}
            </div>
            
            <!-- Interpretation Guide -->
                <div class="section">
                <h2 class="section-title">Understanding the Results</h2>
                
                <h3 style="margin-top: 20px; margin-bottom: 15px; color: #00695c;">What is RRCS?</h3>
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
                
                <h3 style="margin-top: 30px; margin-bottom: 15px; color: #00695c;">Pathogenicity Predictions</h3>
                <div class="info-box">
                    <p><strong>AlphaMissense</strong> uses artificial intelligence to predict whether a mutation will harm protein function:</p>
                    <ul style="margin-left: 20px; margin-top: 5px;">
                        <li><strong style="color: #d32f2f;">Pathogenic (score ≥ 0.564):</strong> Mutation likely damages protein function</li>
                        <li><strong style="color: #f57c00;">Ambiguous (score 0.34-0.563):</strong> Effect uncertain</li>
                        <li><strong style="color: #388e3c;">Benign (score &lt; 0.34):</strong> Mutation likely tolerated</li>
                        </ul>
                    </div>
                
                <h3 style="margin-top: 30px; margin-bottom: 15px; color: #00695c;">Conservation & Population Data</h3>
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
            </div>
        
        <!-- Footer -->
        <div class="footer">
            <img src="data:image/png;base64,{{ lab_logo_base64 }}" alt="Lab Logo" class="footer-logo">
            <p><strong>GPCompaRe</strong> - GPCR Conformational State Comparison Tool</p>
            <p>Report generated on {{ timestamp }}</p>
            <p style="margin-top: 10px; font-size: 0.85em;">
                For questions or issues, please refer to the documentation.
            </p>
        </div>
    </div>
    
    <script>
        // Table search functionality
        function searchTable(tableId) {
            const input = document.getElementById('search-' + tableId.replace('-table', ''));
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
    </script>
        </body>
        </html>
        """)
        
        # Render template
        try:
            report_html = template.render(**self.report_data)
        except Exception as e:
            print(f"Error rendering template: {e}")
            raise
        
        # Save report
        try:
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write(report_html)
                print(f"Report saved to: {report_path}")
        except Exception as e:
            print(f"Error saving report: {e}")
            raise
        
        # Open in browser
        try:
            webbrowser.open(f'file:///{os.path.abspath(report_path)}')
        except Exception as e:
            print(f"Could not open browser: {e}")
        
        return report_path 
