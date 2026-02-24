#!/usr/bin/env python3
"""
Batch GPCR Analysis Pipeline
Systematically processes all GPCR active/inactive pairs and generates comprehensive reports.
"""

import os
import sys
import json
import logging
import traceback
from pathlib import Path
from datetime import datetime
import pandas as pd
from typing import Dict, List, Tuple, Optional
import time
import argparse


class BatchGPCRAnalyzer:
    """Orchestrates batch processing of GPCR conformational state analysis."""
    
    def __init__(self, data_dir: str, output_dir: str, log_level: str = "INFO"):
        """
        Initialize the batch analyzer.
        
        Args:
            data_dir: Directory containing active/inactive PDB files
            output_dir: Root directory for all outputs
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create output structure
        self.setup_output_directories()
        
        # Setup logging
        self.setup_logging(log_level)
        
        # Processing statistics
        self.stats = {
            'total_gpcrs': 0,
            'successful': 0,
            'failed': 0,
            'skipped': 0,
            'start_time': None,
            'end_time': None
        }
        
        self.logger.info(f"Batch GPCR Analyzer initialized")
        self.logger.info(f"Data directory: {self.data_dir}")
        self.logger.info(f"Output directory: {self.output_dir}")
    
    def setup_output_directories(self):
        """Create organized directory structure for outputs."""
        self.dirs = {
            'root': self.output_dir / f"batch_analysis_{self.timestamp}",
            'reports': self.output_dir / f"batch_analysis_{self.timestamp}" / "reports",
            'csv_data': self.output_dir / f"batch_analysis_{self.timestamp}" / "csv_data",
            'plots': self.output_dir / f"batch_analysis_{self.timestamp}" / "plots",
            'logs': self.output_dir / f"batch_analysis_{self.timestamp}" / "logs",
            'summary': self.output_dir / f"batch_analysis_{self.timestamp}" / "summary",
            'rrcs_matrices': self.output_dir / f"batch_analysis_{self.timestamp}" / "rrcs_matrices",
            'variants': self.output_dir / f"batch_analysis_{self.timestamp}" / "variants"
        }
        
        for dir_path in self.dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def setup_logging(self, log_level: str):
        """Configure logging to both file and console."""
        log_file = self.dirs['logs'] / f"batch_analysis_{self.timestamp}.log"
        
        # Create logger
        self.logger = logging.getLogger('BatchGPCRAnalyzer')
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        file_handler.setFormatter(file_formatter)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(getattr(logging, log_level.upper()))
        console_formatter = logging.Formatter(
            '%(levelname)s - %(message)s'
        )
        console_handler.setFormatter(console_formatter)
        
        self.logger.addHandler(file_handler)
        self.logger.addHandler(console_handler)
    
    def discover_gpcr_pairs(self) -> List[Dict[str, str]]:
        """
        Discover all GPCR active/inactive pairs in the data directory.
        
        Returns:
            List of dictionaries containing GPCR information
        """
        self.logger.info("Discovering GPCR pairs...")
        
        active_dir = self.data_dir / "active"
        inactive_dir = self.data_dir / "inactive"
        
        if not active_dir.exists() or not inactive_dir.exists():
            self.logger.error(f"Required directories not found: {active_dir}, {inactive_dir}")
            return []
        
        # Get all active PDB files
        active_files = {f.stem.replace('_active', ''): f 
                       for f in active_dir.glob("*_active.pdb")}
        
        # Get all inactive PDB files
        inactive_files = {f.stem.replace('_inactive', ''): f 
                         for f in inactive_dir.glob("*_inactive.pdb")}
        
        # Find matching pairs
        gpcr_pairs = []
        common_names = set(active_files.keys()) & set(inactive_files.keys())
        
        for gpcr_name in sorted(common_names):
            # Extract receptor code and species
            # Format: receptor_species (e.g., "5ht1a_human")
            parts = gpcr_name.rsplit('_', 1)
            receptor = parts[0] if len(parts) > 1 else gpcr_name
            species = parts[1] if len(parts) > 1 else "unknown"
            
            gpcr_info = {
                'name': gpcr_name,
                'receptor': receptor.upper(),
                'species': species.title(),
                'active_pdb': str(active_files[gpcr_name]),
                'inactive_pdb': str(inactive_files[gpcr_name])
            }
            gpcr_pairs.append(gpcr_info)
        
        self.logger.info(f"Found {len(gpcr_pairs)} GPCR pairs")
        self.stats['total_gpcrs'] = len(gpcr_pairs)
        
        return gpcr_pairs
    
    def process_single_gpcr(self, gpcr_info: Dict[str, str]) -> Dict:
        """
        Process a single GPCR pair.
        
        Args:
            gpcr_info: Dictionary containing GPCR information
            
        Returns:
            Dictionary with processing results
        """
        gpcr_name = gpcr_info['name']
        self.logger.info(f"Processing {gpcr_name}...")
        
        result = {
            'gpcr': gpcr_name,
            'status': 'unknown',
            'error': None,
            'files_generated': []
        }
        
        try:
            # Import the processor (we'll create this next)
            from single_gpcr_processor import SingleGPCRProcessor
            
            processor = SingleGPCRProcessor(
                gpcr_info=gpcr_info,
                output_dirs=self.dirs,
                logger=self.logger
            )
            
            # Run the complete analysis
            processor_result = processor.run_complete_analysis()
            
            result['status'] = 'success'
            result['files_generated'] = processor_result.get('files', [])
            result['summary'] = processor_result.get('summary', {})
            
            self.stats['successful'] += 1
            self.logger.info(f"✓ Successfully processed {gpcr_name}")
            
        except Exception as e:
            result['status'] = 'failed'
            result['error'] = str(e)
            result['traceback'] = traceback.format_exc()
            
            self.stats['failed'] += 1
            self.logger.error(f"✗ Failed to process {gpcr_name}: {str(e)}")
            self.logger.debug(traceback.format_exc())
        
        return result
    
    def process_all_gpcrs(self, gpcr_list: List[Dict[str, str]], 
                         skip_existing: bool = True,
                         delay_between: float = 1.0) -> pd.DataFrame:
        """
        Process all GPCRs in the list.
        
        Args:
            gpcr_list: List of GPCR information dictionaries
            skip_existing: Skip GPCRs with existing reports
            delay_between: Delay in seconds between API calls (rate limiting)
            
        Returns:
            DataFrame with processing results
        """
        self.stats['start_time'] = datetime.now()
        self.logger.info(f"Starting batch processing of {len(gpcr_list)} GPCRs...")
        
        results = []
        
        for idx, gpcr_info in enumerate(gpcr_list, 1):
            self.logger.info(f"\n{'='*60}")
            self.logger.info(f"Progress: {idx}/{len(gpcr_list)} - {gpcr_info['name']}")
            self.logger.info(f"{'='*60}")
            
            # Check if already processed
            if skip_existing:
                report_file = self.dirs['reports'] / f"{gpcr_info['name']}_report.html"
                if report_file.exists():
                    self.logger.info(f"Skipping {gpcr_info['name']} - report already exists")
                    self.stats['skipped'] += 1
                    continue
            
            # Process the GPCR
            result = self.process_single_gpcr(gpcr_info)
            results.append(result)
            
            # Rate limiting delay
            if idx < len(gpcr_list) and delay_between > 0:
                time.sleep(delay_between)
        
        self.stats['end_time'] = datetime.now()
        
        # Convert results to DataFrame
        results_df = pd.DataFrame(results)
        
        # Save results
        results_file = self.dirs['summary'] / "processing_results.csv"
        results_df.to_csv(results_file, index=False)
        self.logger.info(f"Results saved to {results_file}")
        
        return results_df
    
    def generate_master_summary(self, results_df: pd.DataFrame):
        """
        Generate master summary report across all GPCRs.
        
        Args:
            results_df: DataFrame with processing results for all GPCRs
        """
        self.logger.info("Generating master summary...")
        
        # Create summary statistics
        summary = {
            'Analysis Metadata': {
                'Timestamp': self.timestamp,
                'Total GPCRs': self.stats['total_gpcrs'],
                'Successful': self.stats['successful'],
                'Failed': self.stats['failed'],
                'Skipped': self.stats['skipped'],
                'Start Time': self.stats['start_time'].strftime("%Y-%m-%d %H:%M:%S"),
                'End Time': self.stats['end_time'].strftime("%Y-%m-%d %H:%M:%S"),
                'Duration': str(self.stats['end_time'] - self.stats['start_time'])
            }
        }
        
        # Save summary JSON
        summary_file = self.dirs['summary'] / "batch_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        
        self.logger.info(f"Master summary saved to {summary_file}")
        
        # Generate HTML index
        self.generate_html_index(results_df)
        
        # Print final statistics
        self.print_final_stats()
    
    def generate_html_index(self, results_df: pd.DataFrame):
        """Generate a professional database website-style HTML index of all processed GPCRs."""
        from jinja2 import Template

        template = Template("""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GPCR Conformational Analysis Database</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }

        .container {
            max-width: 1600px;
            margin: 0 auto;
        }

        .header {
            background: white;
            padding: 40px;
            border-radius: 12px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
            text-align: center;
        }

        .header h1 {
            font-size: 2.5em;
            color: #2c3e50;
            margin-bottom: 10px;
            font-weight: 300;
        }

        .header .subtitle {
            font-size: 1.2em;
            color: #7f8c8d;
            margin-bottom: 20px;
        }

        .header .timestamp {
            font-size: 0.9em;
            color: #95a5a6;
        }

        .stats-section {
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-bottom: 30px;
        }

        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }

        .stat-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 25px;
            border-radius: 10px;
            color: white;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            transition: transform 0.3s;
        }

        .stat-card:hover {
            transform: translateY(-5px);
        }

        .stat-value {
            font-size: 3em;
            font-weight: bold;
            margin-bottom: 10px;
        }

        .stat-label {
            font-size: 0.95em;
            opacity: 0.9;
            text-transform: uppercase;
            letter-spacing: 1px;
        }

        .database-section {
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
        }

        .section-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            flex-wrap: wrap;
            gap: 15px;
        }

        .section-header h2 {
            color: #2c3e50;
            font-size: 1.8em;
            font-weight: 300;
        }

        .controls {
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
        }

        .search-box {
            padding: 10px 15px;
            border: 2px solid #e0e0e0;
            border-radius: 6px;
            font-size: 1em;
            min-width: 250px;
            transition: border-color 0.3s;
        }

        .search-box:focus {
            outline: none;
            border-color: #667eea;
        }

        .filter-select {
            padding: 10px 15px;
            border: 2px solid #e0e0e0;
            border-radius: 6px;
            font-size: 1em;
            background: white;
            cursor: pointer;
        }

        .btn {
            padding: 10px 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-size: 1em;
            transition: opacity 0.3s;
        }

        .btn:hover {
            opacity: 0.9;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }

        th {
            padding: 15px 12px;
            text-align: left;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: 500;
            cursor: pointer;
            user-select: none;
            position: sticky;
            top: 0;
            z-index: 10;
        }

        th:hover {
            opacity: 0.9;
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
            border-bottom: 1px solid #ecf0f1;
        }

        tbody tr {
            transition: background-color 0.2s;
        }

        tbody tr:hover {
            background-color: #f0f4ff;
        }

        .status-badge {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
            text-transform: uppercase;
        }

        .status-success {
            background: #d4edda;
            color: #155724;
        }

        .status-failed {
            background: #f8d7da;
            color: #721c24;
        }

        .status-skipped {
            background: #e2e3e5;
            color: #6c757d;
        }

        .gpcr-link {
            color: #667eea;
            text-decoration: none;
            font-weight: 600;
            font-size: 1.05em;
        }

        .gpcr-link:hover {
            text-decoration: underline;
        }

        .action-links a {
            color: #667eea;
            text-decoration: none;
            margin-right: 10px;
            font-size: 0.9em;
        }

        .action-links a:hover {
            text-decoration: underline;
        }

        .external-link {
            color: #95a5a6;
            font-size: 0.85em;
            margin-left: 5px;
        }

        .footer {
            background: white;
            padding: 30px;
            border-radius: 12px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            margin-top: 30px;
            text-align: center;
            color: #7f8c8d;
        }

        .footer a {
            color: #667eea;
            text-decoration: none;
        }

        .footer a:hover {
            text-decoration: underline;
        }

        .info-box {
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }

        @media (max-width: 768px) {
            .header h1 {
                font-size: 1.8em;
            }
            .stats-grid {
                grid-template-columns: 1fr;
            }
            .controls {
                width: 100%;
            }
            .search-box {
                width: 100%;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <div class="header">
            <h1>🧬 GPCR Conformational Analysis Database</h1>
            <div class="subtitle">Multistate RRCS Analysis Collection</div>
            <div class="timestamp">Last Updated: {{ timestamp }}</div>
        </div>

        <!-- Statistics Section -->
        <div class="stats-section">
            <h2 style="margin-bottom: 20px; color: #2c3e50;">Database Statistics</h2>
            <div class="stats-grid">
                <div class="stat-card">
                    <div class="stat-value">{{ total }}</div>
                    <div class="stat-label">Total GPCRs</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{{ successful }}</div>
                    <div class="stat-label">Analyzed</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{{ failed }}</div>
                    <div class="stat-label">Failed</div>
                </div>
                <div class="stat-card">
                    <div class="stat-value">{{ skipped }}</div>
                    <div class="stat-label">Skipped</div>
                </div>
            </div>
        </div>

        <!-- Database Table Section -->
        <div class="database-section">
            <div class="section-header">
                <h2>GPCR Receptors</h2>
                <div class="controls">
                    <input type="text" id="searchInput" class="search-box"
                           placeholder="🔍 Search GPCRs..." onkeyup="searchTable()">
                    <select id="filterStatus" class="filter-select" onchange="filterTable()">
                        <option value="all">All Status</option>
                        <option value="success">Success Only</option>
                        <option value="failed">Failed Only</option>
                        <option value="skipped">Skipped Only</option>
                    </select>
                    <button class="btn" onclick="exportToCSV()">📥 Export CSV</button>
                </div>
            </div>

            <div class="info-box">
                <strong>About this database:</strong> This collection contains conformational state analysis for GPCR receptors,
                comparing active and inactive states using Residue-Residue Contact Score (RRCS) methodology. Click on any receptor
                name to view its detailed analysis report.
            </div>

            <table id="gpcrTable">
                <thead>
                    <tr>
                        <th onclick="sortTable(0)">GPCR Name</th>
                        <th onclick="sortTable(1)">Status</th>
                        <th onclick="sortTable(2)">Report</th>
                        <th onclick="sortTable(3)">Data Files</th>
                        <th onclick="sortTable(4)">External Links</th>
                        <th onclick="sortTable(5)">Notes</th>
                    </tr>
                </thead>
                <tbody>
                    {% for result in results %}
                    <tr data-status="{{ result.status }}">
                        <td>
                            {% if result.status == 'success' %}
                            <a href="../reports/{{ result.gpcr }}_report.html" class="gpcr-link">{{ result.gpcr }}</a>
                            {% else %}
                            <strong>{{ result.gpcr }}</strong>
                            {% endif %}
                        </td>
                        <td>
                            <span class="status-badge status-{{ result.status }}">{{ result.status }}</span>
                        </td>
                        <td>
                            {% if result.status == 'success' %}
                            <a href="../reports/{{ result.gpcr }}_report.html" target="_blank">📊 View Report</a>
                            {% else %}
                            <span style="color: #95a5a6;">—</span>
                            {% endif %}
                        </td>
                        <td class="action-links">
                            {% if result.status == 'success' %}
                            <a href="../csv_data/{{ result.gpcr }}_rrcs_delta.csv" download>RRCS</a>
                            <a href="../csv_data/{{ result.gpcr }}_significant_changes.csv" download>Significant</a>
                            <a href="../variants/{{ result.gpcr }}_variants.csv" download>Variants</a>
                            {% else %}
                            <span style="color: #95a5a6;">—</span>
                            {% endif %}
                        </td>
                        <td class="action-links">
                            {% if result.status == 'success' %}
                            <a href="https://gpcrdb.org/protein/{{ result.gpcr|lower }}" target="_blank" title="View in GPCRdb">GPCRdb</a>
                            <a href="https://www.uniprot.org/uniprotkb/?query={{ result.gpcr }}" target="_blank" title="Search UniProt">UniProt</a>
                            {% else %}
                            <span style="color: #95a5a6;">—</span>
                            {% endif %}
                        </td>
                        <td>
                            {% if result.status == 'failed' %}
                            <span style="color: #e74c3c; font-size: 0.9em;">{{ result.error }}</span>
                            {% else %}
                            <span style="color: #95a5a6;">—</span>
                            {% endif %}
                        </td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>

            <p style="margin-top: 20px; color: #7f8c8d; font-size: 0.9em;">
                Showing {{ total }} receptor{{ 's' if total != 1 else '' }}
            </p>
        </div>

        <!-- Footer -->
        <div class="footer">
            <p><strong>GPCR Conformational Analysis Database</strong></p>
            <p style="margin-top: 10px;">
                Powered by <a href="https://github.com/anthropics/claude-code" target="_blank">GPCR Batch Analysis Pipeline</a>
            </p>
            <p style="margin-top: 10px; font-size: 0.85em;">
                RRCS methodology: Residue-Residue Contact Score analysis of AlphaFold multistate models
            </p>
            <p style="margin-top: 15px; font-size: 0.85em;">
                Data sources: AlphaFold, GPCRdb, gnomAD, AlphaMissense, ProtVar, UniProt
            </p>
        </div>
    </div>

    <script>
        // Table sorting
        function sortTable(columnIndex) {
            const table = document.getElementById('gpcrTable');
            const tbody = table.getElementsByTagName('tbody')[0];
            const rows = Array.from(tbody.getElementsByTagName('tr'));
            const th = table.getElementsByTagName('th')[columnIndex];

            const isAsc = th.classList.contains('asc');

            // Clear all sort indicators
            Array.from(table.getElementsByTagName('th')).forEach(header => {
                header.classList.remove('asc', 'desc');
            });

            // Set new sort indicator
            th.classList.add(isAsc ? 'desc' : 'asc');

            // Sort rows
            rows.sort((a, b) => {
                const aValue = a.getElementsByTagName('td')[columnIndex].textContent.trim();
                const bValue = b.getElementsByTagName('td')[columnIndex].textContent.trim();

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

        // Table search
        function searchTable() {
            const input = document.getElementById('searchInput');
            const filter = input.value.toUpperCase();
            const table = document.getElementById('gpcrTable');
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

        // Filter by status
        function filterTable() {
            const select = document.getElementById('filterStatus');
            const filter = select.value;
            const table = document.getElementById('gpcrTable');
            const tr = table.getElementsByTagName('tr');

            for (let i = 1; i < tr.length; i++) {
                const status = tr[i].getAttribute('data-status');
                if (filter === 'all' || status === filter) {
                    tr[i].style.display = '';
                } else {
                    tr[i].style.display = 'none';
                }
            }
        }

        // Export to CSV
        function exportToCSV() {
            const table = document.getElementById('gpcrTable');
            const rows = table.getElementsByTagName('tr');
            const csv = [];

            for (let i = 0; i < rows.length; i++) {
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
            downloadLink.download = 'gpcr_database_{{ timestamp.replace(":", "-").replace(" ", "_") }}.csv';
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
        
        html_content = template.render(
            timestamp=self.timestamp,
            total=self.stats['total_gpcrs'],
            successful=self.stats['successful'],
            failed=self.stats['failed'],
            skipped=self.stats['skipped'],
            results=results_df.to_dict('records')
        )
        
        index_file = self.dirs['summary'] / "index.html"
        with open(index_file, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"HTML index saved to {index_file}")
    
    def print_final_stats(self):
        """Print final processing statistics."""
        duration = self.stats['end_time'] - self.stats['start_time']
        
        self.logger.info("\n" + "="*60)
        self.logger.info("BATCH PROCESSING COMPLETE")
        self.logger.info("="*60)
        self.logger.info(f"Total GPCRs: {self.stats['total_gpcrs']}")
        self.logger.info(f"Successful: {self.stats['successful']}")
        self.logger.info(f"Failed: {self.stats['failed']}")
        self.logger.info(f"Skipped: {self.stats['skipped']}")
        self.logger.info(f"Total Duration: {duration}")
        self.logger.info(f"Average time per GPCR: {duration / self.stats['total_gpcrs'] if self.stats['total_gpcrs'] > 0 else 'N/A'}")
        self.logger.info(f"Output directory: {self.dirs['root']}")
        self.logger.info("="*60)


def main():
    """Main entry point for batch analysis."""
    parser = argparse.ArgumentParser(
        description='Batch GPCR Conformational State Analysis'
    )
    parser.add_argument(
        '--data-dir',
        required=True,
        help='Directory containing active/inactive PDB subdirectories'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Root directory for all outputs'
    )
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level'
    )
    parser.add_argument(
        '--skip-existing',
        action='store_true',
        help='Skip GPCRs with existing reports'
    )
    parser.add_argument(
        '--delay',
        type=float,
        default=1.0,
        help='Delay between API calls in seconds (rate limiting)'
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help='Test mode - process only first 5 GPCRs'
    )
    
    args = parser.parse_args()
    
    # Initialize the batch analyzer
    analyzer = BatchGPCRAnalyzer(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        log_level=args.log_level
    )
    
    # Discover GPCR pairs
    gpcr_list = analyzer.discover_gpcr_pairs()
    
    if not gpcr_list:
        print("No GPCR pairs found!")
        return 1
    
    # Test mode - limit to first 5
    if args.test:
        print("TEST MODE: Processing only first 5 GPCRs")
        gpcr_list = gpcr_list[:5]
    
    # Process all GPCRs
    results_df = analyzer.process_all_gpcrs(
        gpcr_list=gpcr_list,
        skip_existing=args.skip_existing,
        delay_between=args.delay
    )
    
    # Generate master summary
    analyzer.generate_master_summary(results_df)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
