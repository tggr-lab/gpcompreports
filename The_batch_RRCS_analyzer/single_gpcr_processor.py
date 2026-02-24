"""
Single GPCR Processor
Handles complete analysis pipeline for one GPCR active/inactive pair.
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
import json

# Add project directory to path to import existing modules
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

from RRCS import calc_contact
from rrcs_analysis import RRCSAnalyzer

# Optional ouroboros import for snake plot generation
try:
    from ouroboros import fetch_protein_data, SnakePlotRenderer, RenderConfig
    OUROBOROS_AVAILABLE = True
except ImportError:
    OUROBOROS_AVAILABLE = False


class SingleGPCRProcessor:
    """Process a single GPCR active/inactive pair."""
    
    def __init__(self, gpcr_info: Dict, output_dirs: Dict, logger):
        """
        Initialize processor for one GPCR.
        
        Args:
            gpcr_info: Dictionary with GPCR information
            output_dirs: Dictionary of output directories
            logger: Logger instance
        """
        self.gpcr_info = gpcr_info
        self.gpcr_name = gpcr_info['name']
        self.output_dirs = output_dirs
        self.logger = logger
        
        # Initialize RRCS analyzer
        self.analyzer = RRCSAnalyzer()
        
        # Results storage
        self.results = {
            'gpcr_name': self.gpcr_name,
            'receptor': gpcr_info['receptor'],
            'species': gpcr_info['species'],
            'files_generated': []
        }
    
    def run_complete_analysis(self) -> Dict:
        """
        Run the complete analysis pipeline.
        
        Returns:
            Dictionary with analysis results and file paths
        """
        self.logger.info(f"Starting analysis for {self.gpcr_name}")
        
        # Step 1: Calculate RRCS for both states
        self.logger.info("Step 1/7: Calculating RRCS scores...")
        active_rrcs = self.calculate_rrcs(
            self.gpcr_info['active_pdb'], 
            'active'
        )
        inactive_rrcs = self.calculate_rrcs(
            self.gpcr_info['inactive_pdb'], 
            'inactive'
        )
        
        # Step 2: Calculate delta RRCS
        self.logger.info("Step 2/7: Computing RRCS differences...")
        delta_matrix = self.compute_delta_rrcs(active_rrcs, inactive_rrcs)
        
        # Step 3: Get residue information
        self.logger.info("Step 3/7: Extracting residue information...")
        residue_info = self.analyzer.get_residue_info(self.gpcr_info['active_pdb'])
        
        # Step 4: Try to get UniProt ID and fetch additional data
        self.logger.info("Step 4/7: Fetching metadata...")
        uniprot_id = self.try_get_uniprot_id()
        
        tm_domains = None
        gpcrdb_numbers = {}
        if uniprot_id:
            self.analyzer.uniprot_id = uniprot_id
            tm_domains = self.analyzer.fetch_tm_domains(uniprot_id)
            
            # Get GPCRdb numbers for significant residues
            significant_positions = self.get_significant_positions(delta_matrix)
            if significant_positions:
                gpcrdb_numbers = self.analyzer.get_gpcrdb_numbers(
                    uniprot_id, 
                    significant_positions
                )
        
        # Step 5: Fetch variant data (if UniProt ID available)
        self.logger.info("Step 5/7: Fetching variant data...")
        variants_df = self.fetch_variant_data(uniprot_id, significant_positions) if uniprot_id else pd.DataFrame()

        # Step 5.5: Prepare snake plot data (optional, requires ouroboros)
        snake_plot_svg, snake_plot_json = self._prepare_snake_plot_data(
            delta_matrix, variants_df
        )

        # Step 6: Save CSV files
        self.logger.info("Step 6/7: Saving CSV files...")
        self.save_csv_files(
            delta_matrix, 
            variants_df, 
            residue_info, 
            gpcrdb_numbers
        )
        
        # Step 7: Generate report
        self.logger.info("Step 7/7: Generating HTML report...")
        report_path = self.generate_report(
            delta_matrix,
            variants_df,
            residue_info,
            tm_domains,
            gpcrdb_numbers,
            uniprot_id,
            snake_plot_svg=snake_plot_svg,
            snake_plot_json=snake_plot_json
        )
        
        # Compile results summary
        self.results['summary'] = {
            'total_contacts': len(delta_matrix),
            'significant_changes': len(delta_matrix[abs(delta_matrix['delta_rrcs']) >= 3.0]),
            'variants_found': len(variants_df) if not variants_df.empty else 0,
            'uniprot_id': uniprot_id,
            'tm_domains_found': len(tm_domains) if tm_domains else 0
        }
        
        self.logger.info(f"Analysis complete for {self.gpcr_name}")
        self.logger.info(f"  - Total contacts: {self.results['summary']['total_contacts']}")
        self.logger.info(f"  - Significant changes: {self.results['summary']['significant_changes']}")
        self.logger.info(f"  - Variants found: {self.results['summary']['variants_found']}")
        
        return self.results
    
    def calculate_rrcs(self, pdb_file: str, state: str) -> Dict:
        """
        Calculate RRCS for a PDB file.
        
        Args:
            pdb_file: Path to PDB file
            state: 'active' or 'inactive'
            
        Returns:
            Dictionary of RRCS scores
        """
        self.logger.debug(f"Calculating RRCS for {state} state: {pdb_file}")
        
        try:
            contact_scores = calc_contact(pdb_file)
            
            # Save RRCS matrix
            matrix_file = self.output_dirs['rrcs_matrices'] / f"{self.gpcr_name}_{state}_rrcs.csv"
            self.save_rrcs_matrix(contact_scores, matrix_file)
            self.results['files_generated'].append(str(matrix_file))
            
            return contact_scores
            
        except Exception as e:
            self.logger.error(f"Failed to calculate RRCS for {state}: {str(e)}")
            raise
    
    def save_rrcs_matrix(self, contact_scores: Dict, output_file: Path):
        """Save RRCS matrix to CSV."""
        rows = []
        for res1 in contact_scores:
            for res2, score in contact_scores[res1].items():
                if score > 0:
                    # Parse residue identifiers
                    res1_num = int(res1.split(':')[1].split('_')[0])
                    res1_name = res1.split('_')[1]
                    res2_num = int(res2.split(':')[1].split('_')[0])
                    res2_name = res2.split('_')[1]
                    
                    rows.append({
                        'res1': res1_num,
                        'res1_name': res1_name,
                        'res2': res2_num,
                        'res2_name': res2_name,
                        'rrcs_score': score
                    })
        
        df = pd.DataFrame(rows)
        df.to_csv(output_file, index=False)
    
    def compute_delta_rrcs(self, active_rrcs: Dict, inactive_rrcs: Dict) -> pd.DataFrame:
        """
        Compute delta RRCS (active - inactive).
        
        Args:
            active_rrcs: RRCS scores for active state
            inactive_rrcs: RRCS scores for inactive state
            
        Returns:
            DataFrame with delta RRCS values
        """
        delta_data = []
        
        # Collect all unique residue pairs
        all_pairs = set()
        for res1 in active_rrcs:
            for res2 in active_rrcs[res1]:
                all_pairs.add((res1, res2))
        for res1 in inactive_rrcs:
            for res2 in inactive_rrcs[res1]:
                all_pairs.add((res1, res2))
        
        # Calculate delta for each pair
        for res1, res2 in all_pairs:
            active_score = active_rrcs.get(res1, {}).get(res2, 0)
            inactive_score = inactive_rrcs.get(res1, {}).get(res2, 0)
            
            if active_score > 0 or inactive_score > 0:
                # Parse residue identifiers
                res1_num = int(res1.split(':')[1].split('_')[0])
                res1_name = res1.split('_')[1]
                res2_num = int(res2.split(':')[1].split('_')[0])
                res2_name = res2.split('_')[1]
                
                delta_data.append({
                    'res1': res1_num,
                    'res1_name': res1_name,
                    'res2': res2_num,
                    'res2_name': res2_name,
                    'active_rrcs': active_score,
                    'inactive_rrcs': inactive_score,
                    'delta_rrcs': active_score - inactive_score,
                    'abs_delta': abs(active_score - inactive_score)
                })
        
        df = pd.DataFrame(delta_data)
        
        # Sort by absolute delta (most significant changes first)
        df = df.sort_values('abs_delta', ascending=False).reset_index(drop=True)
        
        # Store in analyzer for compatibility with existing code
        self.analyzer.delta_matrix = df
        
        return df
    
    def get_significant_positions(self, delta_matrix: pd.DataFrame, 
                                  threshold: float = 3.0) -> List[int]:
        """
        Get positions with significant RRCS changes.
        
        Args:
            delta_matrix: DataFrame with delta RRCS values
            threshold: Threshold for significance
            
        Returns:
            List of significant positions
        """
        significant = delta_matrix[abs(delta_matrix['delta_rrcs']) >= threshold]
        positions = set()
        positions.update(significant['res1'].tolist())
        positions.update(significant['res2'].tolist())
        return sorted(list(positions))
    
    def try_get_uniprot_id(self) -> Optional[str]:
        """
        Try to determine UniProt ID from GPCR name.
        
        Returns:
            UniProt ID if found, None otherwise
        """
        # Check if there's a mapping file
        if hasattr(self.analyzer, 'uniprot_mapping'):
            # Try to find matching entry
            gpcr_name_upper = self.gpcr_info['receptor'].upper()
            
            for uniprot_id, gpcrdb_id in self.analyzer.uniprot_mapping.items():
                if gpcr_name_upper in gpcrdb_id.upper():
                    self.logger.info(f"Found UniProt ID: {uniprot_id} (GPCRdb: {gpcrdb_id})")
                    return uniprot_id
        
        self.logger.warning(f"No UniProt ID found for {self.gpcr_name}")
        return None
    
    def fetch_variant_data(self, uniprot_id: str, 
                          significant_positions: List[int]) -> pd.DataFrame:
        """
        Fetch variant data from gnomAD for significant positions.
        
        Args:
            uniprot_id: UniProt identifier
            significant_positions: List of positions to query
            
        Returns:
            DataFrame with variant information
        """
        if not uniprot_id or not significant_positions:
            return pd.DataFrame()
        
        try:
            # Get gene name from UniProt
            gene_name = self.analyzer.get_gene_from_uniprot(uniprot_id)
            
            if not gene_name:
                self.logger.warning(f"Could not get gene name for {uniprot_id}")
                return pd.DataFrame()
            
            self.logger.info(f"Fetching variants for gene: {gene_name}")
            
            # Fetch gnomAD data
            gnomad_variants = self.analyzer.get_gnomad_data(gene_name)
            
            if not gnomad_variants:
                self.logger.warning(f"No gnomAD variants found for {gene_name}")
                return pd.DataFrame()
            
            # Filter to significant positions
            variants_list = []
            for variant in gnomad_variants:
                if variant.get('protein_position') in significant_positions:
                    variants_list.append(variant)

            if not variants_list:
                self.logger.info(f"No variants found at significant positions")
                return pd.DataFrame()

            # Enrich variants with AlphaMissense and conservation data
            import re
            aa_codes = {
                'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E',
                'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N',
                'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S',
                'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
            }

            enriched_variants = []
            positions_fetched = set()

            self.logger.info(f"Enriching {len(variants_list)} variants with AlphaMissense data...")
            for variant in variants_list:
                position = variant.get('protein_position')
                hgvsp = variant.get('hgvsp', '')

                # Fetch ProtVar data for this position (only once per position)
                protvar_data = None
                if position and position not in positions_fetched:
                    try:
                        protvar_data = self.analyzer.get_population_data(position)
                        positions_fetched.add(position)
                    except Exception as e:
                        self.logger.debug(f"Could not fetch ProtVar data for position {position}: {e}")

                # Extract amino acid change
                match = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', hgvsp)
                alt_aa_one = None
                if match:
                    ref_aa, pos, alt_aa = match.groups()
                    alt_aa_one = aa_codes.get(alt_aa, alt_aa)

                # Get AlphaMissense score and conservation
                am_score = None
                am_class = 'N/A'
                conservation = None

                if protvar_data and 'score_data' in protvar_data and alt_aa_one:
                    for score in protvar_data['score_data']:
                        if score.get('name') == 'AM' and score.get('mt') == alt_aa_one:
                            am_score = score.get('amPathogenicity')
                            am_class = score.get('amClass', 'N/A')
                        elif score.get('name') == 'CONSERV' and conservation is None:
                            conservation = score.get('score')

                # Add enriched data to variant
                enriched_variant = variant.copy()
                enriched_variant['am_score'] = am_score
                enriched_variant['am_class'] = am_class
                enriched_variant['conservation'] = conservation

                # Calculate heterozygotes from allele counts
                ac = variant.get('ac', 0)
                ac_hom = variant.get('ac_hom', 0)
                enriched_variant['het_count'] = ac - 2 * ac_hom

                enriched_variants.append(enriched_variant)

            df = pd.DataFrame(enriched_variants)
            self.logger.info(f"Found {len(df)} variants at significant positions (enriched with AlphaMissense data)")

            return df
            
        except Exception as e:
            self.logger.error(f"Error fetching variant data: {str(e)}")
            return pd.DataFrame()
    
    def save_csv_files(self, delta_matrix: pd.DataFrame, 
                      variants_df: pd.DataFrame,
                      residue_info: Dict,
                      gpcrdb_numbers: Dict):
        """
        Save all CSV data files.
        
        Args:
            delta_matrix: RRCS delta matrix
            variants_df: Variants DataFrame
            residue_info: Residue information
            gpcrdb_numbers: GPCRdb numbering
        """
        # 1. Delta RRCS matrix
        delta_file = self.output_dirs['csv_data'] / f"{self.gpcr_name}_rrcs_delta.csv"
        delta_matrix.to_csv(delta_file, index=False)
        self.results['files_generated'].append(str(delta_file))
        self.logger.debug(f"Saved delta RRCS to {delta_file}")
        
        # 2. Significant changes (threshold >= 3.0)
        significant = delta_matrix[abs(delta_matrix['delta_rrcs']) >= 3.0]
        sig_file = self.output_dirs['csv_data'] / f"{self.gpcr_name}_significant_changes.csv"
        significant.to_csv(sig_file, index=False)
        self.results['files_generated'].append(str(sig_file))
        
        # 3. Variants (if available)
        if not variants_df.empty:
            var_file = self.output_dirs['variants'] / f"{self.gpcr_name}_variants.csv"
            variants_df.to_csv(var_file, index=False)
            self.results['files_generated'].append(str(var_file))
            self.logger.debug(f"Saved {len(variants_df)} variants to {var_file}")
        
        # 4. Residue annotations
        if gpcrdb_numbers:
            annotations = []
            for pos, info in gpcrdb_numbers.items():
                annotations.append({
                    'position': pos,
                    'amino_acid': info.get('amino_acid', 'X'),
                    'generic_number': info.get('generic_number', ''),
                    'display_number': info.get('display_generic_number', ''),
                    'protein_segment': info.get('protein_segment', ''),
                    'gpcrdb_number': info.get('gpcrdb_number', '')
                })
            
            if annotations:
                ann_df = pd.DataFrame(annotations)
                ann_file = self.output_dirs['csv_data'] / f"{self.gpcr_name}_annotations.csv"
                ann_df.to_csv(ann_file, index=False)
                self.results['files_generated'].append(str(ann_file))
    
    def _prepare_snake_plot_data(self, delta_matrix: pd.DataFrame,
                                  variants_df: pd.DataFrame) -> tuple:
        """
        Prepare snake plot SVG and interactive view data.

        Returns:
            (svg_string, json_data_string) or (None, None) on failure.
        """
        if not OUROBOROS_AVAILABLE:
            self.logger.debug("ouroboros not available, skipping snake plot")
            return None, None

        try:
            from snake_plot_data_builder import SnakePlotDataBuilder

            # Fetch protein data from GPCRdb
            self.logger.debug(f"Fetching GPCRdb data for {self.gpcr_name}...")
            protein_data = fetch_protein_data(self.gpcr_name)

            # Render base SVG (white circles)
            renderer = SnakePlotRenderer(config=RenderConfig())
            svg = renderer.render(protein_data)

            # Extract circle positions for JS arc drawing
            positions = renderer._extract_positions(svg)

            # Calculate significance threshold (same formula as report generator)
            abs_deltas = delta_matrix['delta_rrcs'].abs()
            threshold = max(abs_deltas.mean() + 2 * abs_deltas.std(), 1.0)

            # Build view data
            builder = SnakePlotDataBuilder(
                delta_matrix=delta_matrix,
                variants_df=variants_df,
                protein_data=protein_data,
                significance_threshold=threshold,
            )
            json_data = builder.to_json(positions)

            self.logger.info(f"Snake plot prepared ({len(positions)} residues, "
                             f"{len(protein_data.residues)} total)")
            return svg, json_data

        except Exception as e:
            self.logger.warning(f"Could not prepare snake plot: {e}")
            return None, None

    def generate_report(self, delta_matrix: pd.DataFrame,
                       variants_df: pd.DataFrame,
                       residue_info: Dict,
                       tm_domains: Optional[Dict],
                       gpcrdb_numbers: Dict,
                       uniprot_id: Optional[str],
                       snake_plot_svg: Optional[str] = None,
                       snake_plot_json: Optional[str] = None) -> Path:
        """
        Generate HTML report.

        Args:
            delta_matrix: RRCS delta matrix
            variants_df: Variants DataFrame
            residue_info: Residue information
            tm_domains: TM domain definitions
            gpcrdb_numbers: GPCRdb numbering
            uniprot_id: UniProt identifier

        Returns:
            Path to generated report
        """
        # Try to use improved report generator with plots, fallback to simple if needed
        try:
            from improved_report_generator import ImprovedReportGenerator

            report_gen = ImprovedReportGenerator(
                gpcr_name=self.gpcr_name,
                receptor=self.gpcr_info['receptor'],
                species=self.gpcr_info['species'],
                uniprot_id=uniprot_id,
                delta_matrix=delta_matrix,
                variants_df=variants_df,
                tm_domains=tm_domains,
                gpcrdb_numbers=gpcrdb_numbers
            )

            report_path = self.output_dirs['reports'] / f"{self.gpcr_name}_report.html"
            report_gen.generate(report_path)

        except Exception as e:
            self.logger.warning(f"Could not use improved report generator: {e}")
            self.logger.info("Falling back to simple report generator")

            from simple_report_generator import SimpleReportGenerator

            # Check for lab logo (optional)
            lab_logo_path = None
            potential_logo_paths = [
                'lablogo.png',
                '../lablogo.png',
                str(Path.cwd() / 'lablogo.png'),
                str(Path.home() / 'lablogo.png')
            ]
            for path in potential_logo_paths:
                if Path(path).exists():
                    lab_logo_path = path
                    break

            report_gen = SimpleReportGenerator(
                gpcr_name=self.gpcr_name,
                receptor=self.gpcr_info['receptor'],
                species=self.gpcr_info['species'],
                uniprot_id=uniprot_id,
                delta_matrix=delta_matrix,
                variants_df=variants_df,
                tm_domains=tm_domains,
                gpcrdb_numbers=gpcrdb_numbers,
                lab_logo_path=lab_logo_path,
                snake_plot_svg=snake_plot_svg,
                snake_plot_json=snake_plot_json
            )

            report_path = self.output_dirs['reports'] / f"{self.gpcr_name}_report.html"
            report_gen.generate(report_path)

        self.results['files_generated'].append(str(report_path))
        self.logger.info(f"Report saved to {report_path}")

        return report_path
