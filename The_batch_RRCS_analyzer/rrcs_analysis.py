import pandas as pd
import numpy as np
import os
import sys
from Bio import PDB
from Bio.Data.IUPACData import protein_letters_3to1
import requests
import time
import json
import traceback
import re

# Add the RRCS script directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
rrcs_path = os.path.join(current_dir, 'RRCS-master', 'RRCS-master')
sys.path.append(rrcs_path)

from RRCS import calc_contact

class RRCSAnalyzer:
    def __init__(self):
        self.results = {}
        self.delta_matrix = None
        self.residue_averages = {}
        self.uniprot_mapping = {}
        self.parser = PDB.PDBParser(QUIET=True)
        self.load_uniprot_mapping()
        
    def load_uniprot_mapping(self):
        """Load Uniprot to GPCRdb mapping from file."""
        try:
            with open('uniprot_mapping.txt', 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        uniprot_id, gpcrdb_id = line.strip().split()
                        self.uniprot_mapping[uniprot_id] = gpcrdb_id
        except FileNotFoundError:
            print("Warning: uniprot_mapping.txt not found")
            
    def get_residue_info(self, pdb_file):
        """Extract residue information from PDB file."""
        residue_info = {}
        current_resid = None
        current_resname = None
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    resid = int(line[22:26].strip())  # Remove chain ID
                    resname = line[17:20].strip()
                    
                    if resid != current_resid:
                        current_resid = resid
                        current_resname = resname
                        residue_info[resid] = {
                            'name': resname,
                            'one_letter': self.three_to_one(resname)
                        }
        
        return residue_info
    
    def three_to_one(self, three_letter):
        """Convert three letter amino acid code to one letter."""
        aa_dict = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
            'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
            'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
            'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        return aa_dict.get(three_letter, 'X')
    
    def get_gpcrdb_link(self, uniprot_id, residue_number):
        """Generate GPCRdb link for a residue."""
        if uniprot_id in self.uniprot_mapping:
            gpcrdb_id = self.uniprot_mapping[uniprot_id]
            return f"https://gpcrdb.org/protein/{gpcrdb_id}/residue/{residue_number}"
        return None
    
    def fetch_tm_domains(self, uniprot_id):
        """Fetch transmembrane domain definitions from GPCRdb."""
        if not uniprot_id:
            print("Warning: No Uniprot ID provided")
            return None
            
        if uniprot_id not in self.uniprot_mapping:
            print(f"Warning: No GPCRdb mapping found for Uniprot ID {uniprot_id}")
            print("Please check if the Uniprot ID is correct and exists in uniprot_mapping.txt")
            return None
        
        gpcrdb_id = self.uniprot_mapping[uniprot_id]
        # Use the residues API endpoint instead of topology
        url = f"https://gpcrdb.org/services/residues/{gpcrdb_id}"
        
        try:
            import requests
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                if not data:
                    print(f"Warning: No residue data returned for {gpcrdb_id}")
                    return None
                
                # Extract TM domains from residue data
                tm_domains = {}
                current_tm = None
                start_pos = None
                
                for residue in data:
                    if 'protein_segment' in residue and residue['protein_segment'].startswith('TM'):
                        tm_num = residue['protein_segment'].replace('TM', '')
                        position = residue['sequence_number']
                        
                        if current_tm != tm_num:
                            if current_tm and start_pos:
                                tm_domains[f'TM{current_tm}'] = (start_pos, prev_pos)
                            current_tm = tm_num
                            start_pos = position
                        prev_pos = position
                
                # Add the last TM domain
                if current_tm and start_pos:
                    tm_domains[f'TM{current_tm}'] = (start_pos, prev_pos)
                
                if not tm_domains:
                    print(f"Warning: No TM domains found in residue data for {gpcrdb_id}")
                    return None
                
                print(f"Successfully fetched {len(tm_domains)} TM domains for {gpcrdb_id}")
                for tm, (start, end) in sorted(tm_domains.items()):
                    print(f"{tm}: {start}-{end}")
                return tm_domains
                
            elif response.status_code == 404:
                print(f"Warning: No data found for GPCRdb ID {gpcrdb_id}")
                print("Please check if the protein exists in GPCRdb")
                return None
            else:
                print(f"Warning: Failed to fetch data from GPCRdb (Status code: {response.status_code})")
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Warning: Network error while fetching TM domains: {str(e)}")
            print("Please check your internet connection")
            return None
        except Exception as e:
            print(f"Warning: Unexpected error while fetching TM domains: {str(e)}")
            return None
    
    def get_gpcrdb_numbers(self, uniprot_id, residue_numbers):
        """Get GPCRdb generic residue numbers for a list of residues."""
        if not uniprot_id or uniprot_id not in self.uniprot_mapping:
            return {}
            
        gpcrdb_id = self.uniprot_mapping[uniprot_id]
        url = f"https://gpcrdb.org/services/residues/{gpcrdb_id}"
        
        try:
            import requests
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                if not data:
                    print(f"Warning: No residue data returned for {gpcrdb_id}")
                    return {}
                
                gpcrdb_numbers = {}
                for residue in data:
                    if 'sequence_number' in residue:
                        seq_num = residue['sequence_number']
                        if seq_num in residue_numbers:
                            gpcrdb_numbers[seq_num] = {
                                'generic_number': residue.get('generic_number', ''),
                                'protein_segment': residue.get('protein_segment', 'Other'),
                                'display_generic_number': residue.get('display_generic_number', ''),
                                'amino_acid': residue.get('amino_acid', 'X'),
                                'gpcrdb_number': residue.get('gpcrdb_number', '')
                            }
                
                if not gpcrdb_numbers:
                    print(f"Warning: No matching residues found for {gpcrdb_id}")
                return gpcrdb_numbers
                
            elif response.status_code == 404:
                print(f"Warning: No data found for GPCRdb ID {gpcrdb_id}")
                return {}
            else:
                print(f"Warning: Failed to fetch data from GPCRdb (Status code: {response.status_code})")
                return {}
                
        except Exception as e:
            print(f"Warning: Failed to fetch GPCRdb numbers: {str(e)}")
            return {}
    
    def get_structure(self, pdb_file):
        """Load and return a PDB structure."""
        return self.parser.get_structure('structure', pdb_file)

    def calculate_rrcs(self, pdb_file, state, uniprot_id, chain_id=None):
        """Calculate RRCS scores for a PDB file."""
        try:
            print(f"\n=== Starting RRCS calculation for {pdb_file} ===")
            
            # Extract atom lines from PDB file
            with open(pdb_file, 'r') as f:
                all_lines = f.readlines()
            atom_lines = [l for l in all_lines if l[0:6] == 'ATOM  ']
            print(f"Found {len(atom_lines)} ATOM lines in PDB file")
            
            # Initialize dictionaries
            dict_coord = {}  # dict to store coordinates. dict_coord[res][atom] = (x,y,z,occ)
            atomnum_2_name = {}  # map atom number to atom name
            
            # Process each atom line
            for line in atom_lines:
                # Extract info from atom line
                atom_num = int(line[6:11].strip())
                atom_name = line[12:16].replace(' ', '_')
                res_name = line[17:20]
                res_num = int(line[22:26].strip())
                chain_id_from_pdb = line[21:22]
                
                # Skip if chain_id is specified and doesn't match
                if chain_id and chain_id != chain_id_from_pdb:
                    continue
                    
                # Get coordinates and occupancy
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occ = float(line[54:60].strip())
                
                # Create residue identifier
                res = f"{chain_id_from_pdb}:{res_num}_{res_name}"
                
                # Store atom information
                atomnum_2_name[atom_num] = atom_name
                if res_num <= 0:
                    continue
                if res not in dict_coord:
                    dict_coord[res] = {}
                dict_coord[res][atom_num] = (x, y, z, occ)
            
            print(f"\nProcessed {len(dict_coord)} residues")
            
            # Calculate RRCS scores
            rrcs_scores = []
            res_list = sorted(dict_coord.keys())
            
            # Process all residue pairs
            for i, res1 in enumerate(res_list):
                res1_num = int(res1.split(':')[1].split('_')[0])
                atoms_in_res1 = dict_coord[res1]
                
                for res2 in res_list[i+1:]:
                    res2_num = int(res2.split(':')[1].split('_')[0])
                    atoms_in_res2 = dict_coord[res2]
                    
                    # Quick distance check first
                    contact_found = False
                    for a1 in atoms_in_res1.values():
                        if contact_found:
                            break
                        for a2 in atoms_in_res2.values():
                            dx = abs(a1[0] - a2[0])
                            dy = abs(a1[1] - a2[1])
                            dz = abs(a1[2] - a2[2])
                            if dx < 4.63 and dy < 4.63 and dz < 4.63:
                                contact_found = True
                                break
                    
                    if not contact_found:
                        continue
                    
                    # Calculate detailed RRCS score
                    total_score = 0
                    close_residues = abs(res1_num - res2_num) < 5
                    
                    if close_residues:
                        # For close residues, exclude backbone atoms
                        for a1_num, a1_coord in atoms_in_res1.items():
                            a1_name = atomnum_2_name[a1_num]
                            if a1_name in ['_N__', '_CA_', '_C__', '_O__']:
                                continue
                            for a2_num, a2_coord in atoms_in_res2.items():
                                a2_name = atomnum_2_name[a2_num]
                                if a2_name in ['_N__', '_CA_', '_C__', '_O__']:
                                    continue
                                d2 = sum((c1-c2)**2 for c1, c2 in zip(a1_coord[:3], a2_coord[:3]))
                                if d2 >= 21.4369:  # 4.63^2
                                    score = 0
                                elif d2 <= 10.4329:  # 3.23^2
                                    score = 1.0 * a1_coord[3] * a2_coord[3]
                                else:
                                    score = (1 - (d2**0.5 - 3.23)/1.4) * a1_coord[3] * a2_coord[3]
                                total_score += score
                    else:
                        # For distant residues, include all atoms
                        for a1_coord in atoms_in_res1.values():
                            for a2_coord in atoms_in_res2.values():
                                d2 = sum((c1-c2)**2 for c1, c2 in zip(a1_coord[:3], a2_coord[:3]))
                                if d2 >= 21.4369:  # 4.63^2
                                    score = 0
                                elif d2 <= 10.4329:  # 3.23^2
                                    score = 1.0 * a1_coord[3] * a2_coord[3]
                                else:
                                    score = (1 - (d2**0.5 - 3.23)/1.4) * a1_coord[3] * a2_coord[3]
                                total_score += score
                    
                    if total_score > 0:
                        rrcs_scores.append({
                            'res1': res1_num,
                            'res2': res2_num,
                            'res1_name': res1.split('_')[1],
                            'res2_name': res2.split('_')[1],
                            'rrcs': total_score,
                            'state': state
                        })
            
            # Convert to DataFrame
            df = pd.DataFrame(rrcs_scores)
            if len(df) == 0:
                print(f"Warning: No RRCS scores calculated for {pdb_file}")
                return None
            
            print(f"\n=== Final Results ===")
            print(f"Total contacts found: {len(df)}")
            
            return df
            
        except Exception as e:
            print(f"Error calculating RRCS for {pdb_file}: {str(e)}")
            import traceback
            traceback.print_exc()
            return None
    
    def calculate_deltas(self, reference_state='active'):
        """Calculate differences in RRCS scores between states using a reference state."""
        if not self.results or reference_state not in self.results:
            raise ValueError(f"Missing results for reference state: {reference_state}")
        
        # Get reference state data
        ref_df = self.results[reference_state]
        if isinstance(ref_df, dict):
            return  # Skip if it's the statistics dictionary
        
        # Calculate deltas for all other states
        delta_dfs = []
        for state, df in self.results.items():
            if state != reference_state and state != 'statistics' and isinstance(df, pd.DataFrame):
                # Merge with reference state
                merged = pd.merge(
                    ref_df, df,
                    on=['residue1', 'residue2'],
                    suffixes=(f'_{reference_state}', f'_{state}')
                )
                
                # Calculate delta
                merged['delta'] = merged[f'score_{reference_state}'] - merged[f'score_{state}']
                merged['comparison'] = f"{reference_state}_vs_{state}"
                delta_dfs.append(merged)
        
        # Combine all delta DataFrames
        if delta_dfs:
            self.delta_matrix = pd.concat(delta_dfs, ignore_index=True)
            return self.delta_matrix
        else:
            raise ValueError("No other states to compare with reference state")
    
    def calculate_residue_averages(self, state):
        """Calculate average RRCS scores for each residue."""
        if state not in self.results:
            raise ValueError(f"No results for state: {state}")
        
        df = self.results[state]
        if isinstance(df, dict):
            return  # Skip if it's the statistics dictionary
        
        # Calculate averages for residue1
        avg1 = df.groupby('residue1')['score'].mean()
        # Calculate averages for residue2
        avg2 = df.groupby('residue2')['score'].mean()
        
        # Combine and average where residue appears in both positions
        combined = pd.concat([avg1, avg2]).groupby(level=0).mean()
        
        self.residue_averages[state] = combined
        return combined
    
    def get_significant_changes(self, threshold=0.5):
        """Get residue pairs with significant RRCS score changes."""
        if self.delta_matrix is None:
            raise ValueError("Please calculate deltas first")
        
        return self.delta_matrix[abs(self.delta_matrix['delta']) >= threshold]
    
    def get_transmembrane_changes(self):
        """Get changes in RRCS scores between transmembrane domains."""
        if not hasattr(self, 'delta_matrix') or self.delta_matrix is None:
            return None
            
        if not hasattr(self, 'tm_definitions') or self.tm_definitions is None:
            return None
            
        # Function to get TM domain for a residue number
        def get_tm_domain(res_num):
            try:
                num = int(str(res_num))
                for domain, (start, end) in self.tm_definitions.items():
                    if start <= num <= end:
                        return domain
                return 'Other'
            except ValueError:
                return 'Other'
        
        # Add domain information
        result = self.delta_matrix.copy()
        result['domain1'] = result['res1'].apply(get_tm_domain)
        result['domain2'] = result['res2'].apply(get_tm_domain)
        
        return result
    
    def get_contact_statistics(self):
        """Calculate detailed contact statistics."""
        stats = {}
        
        for state, df in self.results.items():
            if state != 'statistics':  # Skip the statistics dictionary
                state_stats = {
                    'total_contacts': len(df),
                    'unique_residues': len(set(df['res1_num'].tolist() + df['res2_num'].tolist())),
                    'score_distribution': {
                        'mean': df['score'].mean(),
                        'std': df['score'].std(),
                        'median': df['score'].median(),
                        'q1': df['score'].quantile(0.25),
                        'q3': df['score'].quantile(0.75),
                        'min': df['score'].min(),
                        'max': df['score'].max()
                    }
                }
                
                # Segment analysis
                if 'res1_segment' in df.columns:
                    segment_contacts = df.groupby(['res1_segment', 'res2_segment']).agg({
                        'score': ['count', 'mean', 'std']
                    }).round(3)
                    state_stats['segment_analysis'] = segment_contacts.to_dict()
                
                # Residue type analysis
                residue_contacts = df.groupby(['res1_aa']).size().to_dict()
                state_stats['residue_composition'] = residue_contacts
                
                # Score thresholds
                for threshold in [0.3, 0.5, 0.7]:
                    state_stats[f'contacts_above_{threshold}'] = len(df[df['score'] >= threshold])
                
                stats[state] = state_stats
        
        if self.delta_matrix is not None:
            delta_stats = {
                'total_comparisons': len(self.delta_matrix),
                'mean_change': self.delta_matrix['delta'].mean(),
                'std_change': self.delta_matrix['delta'].std(),
                'max_increase': self.delta_matrix['delta'].max(),
                'max_decrease': self.delta_matrix['delta'].min(),
                'significant_changes': len(self.delta_matrix[abs(self.delta_matrix['delta']) >= 0.5])
            }
            stats['conformational_changes'] = delta_stats
        
        return stats 
    
    def calculate_changes(self, active_rrcs, inactive_rrcs):
        """Calculate differences in RRCS scores between active and inactive states."""
        print("\nActive RRCS columns:")
        print(active_rrcs.columns.tolist())
        print("\nInactive RRCS columns:")
        print(inactive_rrcs.columns.tolist())
        
        # Merge the active and inactive state data
        merged = pd.merge(
            active_rrcs, inactive_rrcs,
            on=['res1', 'res2'],  # Only merge on residue numbers
            suffixes=('_active', '_inactive')
        )
        
        print("\nMerged columns:")
        print(merged.columns.tolist())
        print("\nFirst few rows of merged data:")
        print(merged.head())
        
        # Calculate delta (active - inactive)
        merged['delta_rrcs'] = merged['rrcs_active'] - merged['rrcs_inactive']
        
        # Add additional information
        merged['active_rrcs'] = merged['rrcs_active']
        merged['inactive_rrcs'] = merged['rrcs_inactive']
        
        # Store the delta matrix in the class instance
        self.delta_matrix = merged
        
        # Print column names for debugging
        print("\nFinal delta matrix columns:")
        print(self.delta_matrix.columns.tolist())
        
        return merged
    
    def set_uniprot_id(self, uniprot_id):
        """Set the Uniprot ID for the analysis and fetch TM definitions."""
        self.uniprot_id = uniprot_id
        if uniprot_id not in self.uniprot_mapping:
            print(f"Warning: No GPCRdb mapping found for Uniprot ID {uniprot_id}")
            print("Please check if the Uniprot ID is correct and exists in uniprot_mapping.txt")
            return
            
        # Fetch TM definitions
        print("\nFetching TM definitions...")
        self.tm_definitions = self.fetch_tm_domains(uniprot_id)
        if self.tm_definitions:
            print("Successfully loaded TM definitions")
        else:
            print("Failed to load TM definitions")

    def get_protvar_data(self, threshold=3.5):
        """
        Get ProtVar data for residues with significant RRCS changes.
        
        Args:
            threshold (float): ΔRRCS threshold for significant changes (default: 3.5)
        """
        if self.delta_matrix is None:
            raise ValueError("Please calculate deltas first using calculate_deltas()")
            
        print("\n=== Starting ProtVar Analysis ===")
        print(f"Using threshold: {threshold}")
        
        # Get residues with significant changes
        significant_changes = self.delta_matrix[abs(self.delta_matrix['delta_rrcs']) >= threshold]
        print(f"\nFound {len(significant_changes)} significant changes above threshold")
        
        if len(significant_changes) == 0:
            print(f"No residues found with ΔRRCS >= {threshold}")
            return pd.DataFrame()
            
        # Initialize results dataframe
        results = []
        
        # Get unique residues
        unique_residues = set()
        for _, row in significant_changes.iterrows():
            # Add residue with amino acid name from active state
            unique_residues.add((row['res1'], row['res1_name_active']))
            unique_residues.add((row['res2'], row['res2_name_active']))
        
        print(f"\nFetching ProtVar data for {len(unique_residues)} unique residues...")
        
        # Ensure UniProt ID is uppercase
        uniprot_id = self.uniprot_id.upper() if self.uniprot_id else None
        print("UniProt ID:", uniprot_id)
        
        # Base URL for ProtVar API
        base_url = "https://www.ebi.ac.uk/ProtVar/api"
        print(f"Using ProtVar API base URL: {base_url}")
        
        # Headers for API requests
        headers = {
            'Accept': 'application/json'
        }
        print("Request headers:", headers)
        
        for res_num, res_name in unique_residues:
            print(f"\n--- Processing residue {res_name}{res_num} ---")
            try:
                # Convert three-letter code to one-letter code
                orig_aa = self.three_to_one(res_name)
                if orig_aa == 'X':
                    print(f"Warning: Could not convert residue name {res_name} to one-letter code")
                    continue
                
                print(f"Converted {res_name} to {orig_aa}")
                
                # Get scores (Conservation, EVE, ESM1b, AlphaMissense)
                print("\nFetching scores...")

                # First get conservation and EVE scores (no mt needed)
                score_url = f"{base_url}/score/{self.uniprot_id}/{res_num}"
                print(f"Fetching base scores from: {score_url}")
                score_response = requests.get(score_url, headers=headers)

                print(f"Base score response status: {score_response.status_code}")
                if score_response.status_code != 200:
                    print(f"Error response content: {score_response.text}")
                    print(f"Error fetching score data for residue {res_name}{res_num}")
                    continue

                score_data = score_response.json()
                print("Base score data received:", json.dumps(score_data, indent=2))

                # Process conservation and EVE scores
                conservation_score = None
                eve_score = None

                print("\nProcessing base scores...")
                for score in score_data:
                    if score['name'] == 'CONSERV':
                        conservation_score = score.get('score')
                        print(f"Conservation score from ProtVar: {conservation_score}")
                    elif score['name'] == 'EVE':
                        eve_score = score.get('score')
                        print(f"EVE score: {eve_score}")

                # Get max delta RRCS for this residue
                residue_changes = significant_changes[
                    (significant_changes['res1'] == res_num) |
                    (significant_changes['res2'] == res_num)
                ]
                max_delta = abs(residue_changes['delta_rrcs']).max()
                print(f"Max ΔRRCS: {max_delta}")

                # Initialize result dictionary with basic info
                result = {
                    'residue': f"{res_name}{res_num}",
                    'wt_aa': orig_aa,
                    'position': res_num,
                    'conservation_score': conservation_score,  # Use actual value from API
                    'eve_score': eve_score,  # Use actual value from API
                    'max_delta': max_delta
                }

                # Get FoldX data
                foldx_url = f"{base_url}/foldx/{self.uniprot_id}/{res_num}"
                print(f"\nFetching FoldX data from: {foldx_url}")
                foldx_response = requests.get(foldx_url, headers=headers)

                print(f"FoldX response status: {foldx_response.status_code}")
                if foldx_response.status_code == 200:
                    foldx_data = foldx_response.json()
                    print("FoldX data received:", json.dumps(foldx_data, indent=2))
                    
                    # Process FoldX data
                    for foldx in foldx_data:
                        mt_aa = foldx['mutatedType']
                        result[f'FoldX_{mt_aa}'] = foldx['foldxDdg']
                        if 'plddt' not in result:
                            result['plddt'] = foldx['plddt']

                # Extract AM scores from base response
                print("\nExtracting AM scores from base response...")
                for score in score_data:
                    if score['name'] == 'AM' and score.get('mt'):
                        mt_aa = score.get('mt')
                        result[f'AM_{mt_aa}'] = score.get('amPathogenicity')
                        result[f'AMClass_{mt_aa}'] = score.get('amClass')
                        print(f"Found AM score for {mt_aa}: {score.get('amPathogenicity')} ({score.get('amClass')})")

                results.append(result)
                print(f"\nSuccessfully processed residue {res_name}{res_num}")
                
            except Exception as e:
                print(f"\nError processing residue {res_name}{res_num}:")
                print(f"Exception type: {type(e).__name__}")
                print(f"Exception message: {str(e)}")
                print("Traceback:")
                traceback.print_exc()
                continue
        
        print("\n=== Converting results to DataFrame ===")
        # Convert to DataFrame
        df = pd.DataFrame(results)
        
        if not df.empty:
            print("\nConverting score columns to float...")
            # Convert numeric columns to float
            numeric_columns = [col for col in df.columns if col.startswith(('AM_', 'FoldX_')) or col in ['conservation_score', 'eve_score', 'max_delta', 'plddt']]
            for col in numeric_columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                print(f"Column {col} unique values: {df[col].unique()}")

        print(f"\n=== Analysis Complete ===")
        print(f"Final DataFrame shape: {df.shape}")
        print(f"Columns: {df.columns.tolist()}")
        
        return df 

    def test_protvar_query(self, residue_num, residue_name):
        """
        Test function to debug ProtVar API queries for a single residue.
        
        Args:
            residue_num (int): Residue number
            residue_name (str): Three-letter residue name
        """
        print(f"\n=== Testing ProtVar Query for {residue_name}{residue_num} ===")
        
        # Base URL for ProtVar API
        base_url = "https://www.ebi.ac.uk/ProtVar/api"
        headers = {'Accept': 'application/json'}
        
        # Convert three-letter code to one-letter code
        orig_aa = self.three_to_one(residue_name)
        print(f"\nConverted {residue_name} to {orig_aa}")
        
        # Test just a few mutations
        test_mutations = ['A', 'D', 'K', 'W']  # Testing just 4 mutations
        
        # 1. Get conservation score
        print("\n1. Testing conservation score...")
        score_url = f"{base_url}/score/{self.uniprot_id}/{residue_num}"
        score_response = requests.get(score_url, headers=headers)
        print(f"URL: {score_url}")
        print(f"Response status: {score_response.status_code}")
        
        conservation_score = None
        eve_score = None
        if score_response.status_code == 200:
            score_data = score_response.json()
            print("Score data:")
            print(json.dumps(score_data, indent=2))
            
            for score in score_data:
                if score['name'] == 'CONSERV':
                    conservation_score = score.get('score')
                    print(f"Conservation score from ProtVar: {conservation_score}")
                elif score['name'] == 'EVE':
                    eve_score = score.get('score')
                    print(f"EVE score: {eve_score}")
        
        # 2. Get AM scores for test mutations
        print("\n2. Testing AM scores for specific mutations...")
        am_scores = {}
        for mt_aa in test_mutations:
            if mt_aa != orig_aa:
                am_url = f"{base_url}/score/{self.uniprot_id}/{residue_num}?mt={mt_aa}&name=AM"
                print(f"\nTesting mutation to {mt_aa}")
                print(f"URL: {am_url}")
                
                am_response = requests.get(am_url, headers=headers)
                print(f"Response status: {am_response.status_code}")
                
                if am_response.status_code == 200:
                    am_data = am_response.json()
                    print("Raw AM response:")
                    print(json.dumps(am_data, indent=2))
                    
                    # Extract AM score from amPathogenicity field
                    for score in am_data:
                        if score['name'] == 'AM' and score.get('mt') == mt_aa:
                            am_scores[mt_aa] = {
                                'score': score.get('amPathogenicity'),
                                'class': score.get('amClass')
                            }
                            print(f"Found AM score for {mt_aa}: {am_scores[mt_aa]['score']} ({am_scores[mt_aa]['class']})")
        
        # 3. Get FoldX data for test mutations
        print("\n3. Testing FoldX data...")
        foldx_url = f"{base_url}/foldx/{self.uniprot_id}/{residue_num}"
        foldx_response = requests.get(foldx_url, headers=headers)
        print(f"URL: {foldx_url}")
        print(f"Response status: {foldx_response.status_code}")
        
        foldx_data = {}
        if foldx_response.status_code == 200:
            foldx_raw = foldx_response.json()
            for foldx in foldx_raw:
                mt = foldx['mutatedType']
                if mt in test_mutations:
                    foldx_data[mt] = {
                        'ddg': foldx['foldxDdg'],
                        'plddt': foldx['plddt']
                    }
                    print(f"FoldX for {orig_aa}{residue_num}{mt}: ΔΔG = {foldx_data[mt]['ddg']}, pLDDT = {foldx_data[mt]['plddt']}")
        
        print("\n=== Summary of Results ===")
        print(f"Conservation score: {conservation_score}")
        print(f"EVE score: {eve_score}")
        print("\nAM Scores:")
        for mt, data in am_scores.items():
            print(f"{orig_aa}{residue_num}{mt}: {data['score']} ({data['class']})")
        
        print("\nFoldX Data:")
        for mt, data in foldx_data.items():
            print(f"{orig_aa}{residue_num}{mt}: ΔΔG = {data['ddg']}, pLDDT = {data['plddt']}")
        
        return {
            'conservation': conservation_score,
            'eve': eve_score,
            'am_scores': am_scores,
            'foldx_data': foldx_data
        } 

    def get_population_data(self, position):
        """Get population data for a position from ProtVar API."""
        if not self.uniprot_id:
            return None
            
        data = {}
        
        # Get population data
        pop_url = f"https://www.ebi.ac.uk/ProtVar/api/population/{self.uniprot_id}/{position}"
        try:
            response = requests.get(pop_url)
            if response.status_code == 200:
                data['population_data'] = response.json()
        except Exception as e:
            print(f"Error fetching population data: {e}")
        
        # Get score predictions
        score_url = f"https://www.ebi.ac.uk/ProtVar/api/score/{self.uniprot_id}/{position}"
        try:
            response = requests.get(score_url)
            if response.status_code == 200:
                data['score_data'] = response.json()
        except Exception as e:
            print(f"Error fetching score data: {e}")
        
        # Get function annotations
        function_url = f"https://www.ebi.ac.uk/ProtVar/api/function/{self.uniprot_id}/{position}"
        try:
            response = requests.get(function_url)
            if response.status_code == 200:
                data['function_data'] = response.json()
        except Exception as e:
            print(f"Error fetching function data: {e}")
        
        return data 

    def get_gnomad_data(self, gene_name):
        """Get variant data from gnomAD for a given gene."""
        print(f"Getting gnomAD data for gene: {gene_name}")
        
        # First get the Ensembl ID
        ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"
        response = requests.get(ensembl_url)
        if response.status_code != 200:
            print(f"Error getting Ensembl ID: {response.status_code}")
            return []
        
        ensembl_data = response.json()
        ensembl_id = ensembl_data.get('id')
        if not ensembl_id:
            print("No Ensembl ID found")
            return []
        
        print(f"Found Ensembl ID: {ensembl_id}")
        
        # Try gnomAD v4 first (best coverage with both exomes and genomes)
        # v3 is fallback for genomes-only data
        datasets_to_try = [
            ("gnomad_r4", "gnomAD v4 (exomes + genomes)"),
            ("gnomad_r3", "gnomAD v3 (genomes only)"),
        ]
        
        gnomad_url = "https://gnomad.broadinstitute.org/api"
        query = """
        query GeneVariants($geneId: String!, $datasetId: DatasetId!) {
          gene(gene_id: $geneId, reference_genome: GRCh38) {
            variants(dataset: $datasetId) {
              variant_id
              pos
              ref
              alt
              consequence
              hgvsc
              hgvsp
              rsids
              exome {
                ac
                an
                af
                ac_hom
                filters
                populations {
                  id
                  ac
                  an
                  ac_hom
                }
              }
              genome {
                ac
                an
                af
                ac_hom
                filters
                populations {
                  id
                  ac
                  an
                  ac_hom
                }
              }
              flags
              transcript_consequence {
                lof
                lof_filter
                major_consequence
                hgvsc
                hgvsp
                is_canonical
              }
            }
          }
        }
        """
        
        variants = []
        dataset_used = None
        
        for dataset_id, dataset_name in datasets_to_try:
            print(f"Trying gnomAD {dataset_name}...")
            variables = {
                "geneId": ensembl_id,
                "datasetId": dataset_id
            }
            
            response = requests.post(gnomad_url, json={"query": query, "variables": variables})
            if response.status_code != 200:
                print(f"  Error querying gnomAD {dataset_name}: {response.status_code}")
                continue
            
            data = response.json()
            if "errors" in data:
                print(f"  GraphQL errors for {dataset_name}:")
                print(json.dumps(data["errors"], indent=2))
                continue
            
            variants = data.get("data", {}).get("gene", {}).get("variants", [])
            if variants:
                dataset_used = dataset_name
                print(f"  Successfully retrieved data from gnomAD {dataset_name}")
                break
            else:
                print(f"  No variants found in {dataset_name}")
        
        if not variants:
            print("No variants found in any gnomAD dataset")
            return []
        
        print(f"\nUsing gnomAD {dataset_used}")
        print(f"Found {len(variants)} total variants")
        
        # Store which dataset was actually used for reporting
        self.gnomad_dataset_used = dataset_used
        
        # Filter for missense variants
        missense_variants = [v for v in variants if v.get('consequence') == 'missense_variant']
        print(f"Found {len(missense_variants)} missense variants")
        
        # Process variants into a list of dictionaries with relevant information
        processed_variants = []
        for variant in missense_variants:
            tc = variant.get('transcript_consequence', {})
            if not tc:
                continue
            
            # Extract protein position from HGVS.p notation
            hgvsp = tc.get('hgvsp')
            if not hgvsp:
                continue
            
            # Extract protein position and amino acid changes
            match = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', hgvsp)
            if not match:
                continue
                
            ref_aa, pos, alt_aa = match.groups()
            protein_position = int(pos)
            
            # gnomAD v4 has both exome and genome data - combine them
            # Priority: use exome if available (more samples), otherwise genome
            exome = variant.get("exome", {})
            genome = variant.get("genome", {})
            
            # Choose primary dataset (exome preferred as it has more samples in v4)
            if exome and exome.get("ac"):
                primary_data = exome
                data_source = "exome"
            elif genome and genome.get("ac"):
                primary_data = genome
                data_source = "genome"
            else:
                continue  # No data available
            
            af = primary_data.get("af")
            ac = primary_data.get("ac")
            an = primary_data.get("an")
            ac_hom = primary_data.get("ac_hom", 0)
            
            if not all(x is not None for x in [af, ac, an]):
                continue
            
            # Get population frequencies from primary data source
            pop_freqs = {}
            for pop in primary_data.get('populations', []):
                pop_ac = pop.get('ac', 0)
                pop_an = pop.get('an', 0)
                pop_hom = pop.get('ac_hom', 0)
                if pop_an > 0:
                    pop_af = pop_ac / pop_an
                    pop_freqs[pop['id']] = {
                        'af': pop_af,
                        'ac': pop_ac,
                        'an': pop_an,
                        'ac_hom': pop_hom
                    }
            
            # Get the DNA change (hgvsc) from transcript consequence if available, otherwise from variant
            hgvsc = tc.get('hgvsc') or variant.get('hgvsc', '')
            
            variant_info = {
                'variant_id': variant.get('variant_id'),
                'position': variant['pos'],  # Keep genomic position
                'protein_position': protein_position,  # Add protein position
                'ref': variant['ref'],
                'alt': variant['alt'],
                'ref_aa': ref_aa,  # Add reference amino acid
                'alt_aa': alt_aa,  # Add alternate amino acid
                'hgvsp': hgvsp,
                'hgvsc': hgvsc,  # Add DNA change
                'af': af,
                'ac': ac,
                'an': an,
                'ac_hom': ac_hom,
                'population_frequencies': pop_freqs,
                'rsids': variant.get('rsids', []),
                'data_source': data_source  # Track whether from exome or genome
            }
            
            processed_variants.append(variant_info)
        
        return processed_variants

    def get_gene_from_uniprot(self, uniprot_id):
        """Get gene name from UniProt ID.
        
        Args:
            uniprot_id (str): UniProt ID
            
        Returns:
            str: Gene name/symbol
        """
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            # Get gene name from response
            if 'genes' in data:
                for gene in data['genes']:
                    if 'geneName' in gene:
                        return gene['geneName']['value']
            return None
        except Exception as e:
            print(f"Error fetching gene name from UniProt: {str(e)}")
            return None 

    def get_coordinate_mapping(self, uniprot_id):
        """Get mapping between protein positions and genomic coordinates."""
        try:
            # First get the gene name from UniProt
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
            response = requests.get(url)
            if not response.ok:
                print(f"Error getting UniProt data: {response.status_code}")
                return None
            
            data = response.json()
            
            # Get gene name
            gene_name = None
            if 'genes' in data:
                for gene in data['genes']:
                    if 'geneName' in gene:
                        gene_name = gene['geneName']['value']
                        break
            
            if not gene_name:
                print("Could not find gene name in UniProt data")
                return None
            
            # Get protein sequence length
            sequence_length = len(data.get('sequence', {}).get('value', ''))
            
            # Get genomic coordinates from Ensembl
            ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?"
            headers = {"Content-Type": "application/json"}
            response = requests.get(ensembl_url, headers=headers)
            
            if not response.ok:
                print(f"Error getting Ensembl data: {response.status_code}")
                return None
            
            ensembl_data = response.json()
            
            return {
                'genomic_start': ensembl_data.get('start'),
                'genomic_end': ensembl_data.get('end'),
                'protein_start': 1,  # Assume 1-based protein coordinates
                'sequence_length': sequence_length,
                'is_reverse': ensembl_data.get('strand') == -1
            }
            
        except Exception as e:
            print(f"Error getting coordinate mapping: {str(e)}")
            return None 

    def test_gnomad_datasets(self, gene_name):
        """Test different gnomAD datasets to compare allele counts."""
        print(f"\nTesting gnomAD datasets for {gene_name}...")
        
        # First get the Ensembl ID
        ensembl_url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_name}?content-type=application/json"
        response = requests.get(ensembl_url)
        if response.status_code != 200:
            print(f"Error getting Ensembl ID: {response.status_code}")
            return
        
        ensembl_data = response.json()
        ensembl_id = ensembl_data.get('id')
        if not ensembl_id:
            print("No Ensembl ID found")
            return
        
        print(f"Found Ensembl ID: {ensembl_id}")
        
        # Test different datasets
        datasets = [
            "gnomad_r2_1",  # v2.1
            "gnomad_r3",    # v3
            "gnomad_r4",    # v4 with UKB
            "gnomad_r4_non_ukb"  # v4 without UKB
        ]
        
        query = """
        query GeneVariants($geneId: String!, $datasetId: DatasetId!) {
          gene(gene_id: $geneId, reference_genome: GRCh38) {
            variants(dataset: $datasetId) {
              variantId
              pos
              ref
              alt
              consequence
              genome {
                ac
                an
                af
                ac_hom
              }
            }
          }
        }
        """
        
        results = {}
        for dataset in datasets:
            print(f"\nTesting dataset: {dataset}")
            variables = {
                "geneId": ensembl_id,
                "datasetId": dataset
            }
            
            response = requests.post(
                "https://gnomad.broadinstitute.org/api",
                json={"query": query, "variables": variables}
            )
            
            if response.status_code == 200:
                data = response.json()
                if "errors" in data:
                    print(f"GraphQL errors for {dataset}:")
                    print(json.dumps(data["errors"], indent=2))
                    continue
                    
                variants = data.get("data", {}).get("gene", {}).get("variants", [])
                missense_variants = [v for v in variants if v.get('consequence') == 'missense_variant']
                
                total_ac = sum(v.get('genome', {}).get('ac', 0) for v in missense_variants)
                total_an = sum(v.get('genome', {}).get('an', 0) for v in missense_variants)
                total_hom = sum(v.get('genome', {}).get('ac_hom', 0) for v in missense_variants)
                
                results[dataset] = {
                    'total_variants': len(variants),
                    'missense_variants': len(missense_variants),
                    'total_alleles': total_ac,
                    'total_samples': total_an // 2 if total_an > 0 else 0,
                    'total_hom': total_hom
                }
                
                print(f"Total variants: {len(variants)}")
                print(f"Missense variants: {len(missense_variants)}")
                print(f"Total allele count: {total_ac}")
                print(f"Total samples: {total_an // 2 if total_an > 0 else 0}")
                print(f"Total homozygotes: {total_hom}")
            else:
                print(f"Error querying {dataset}: {response.status_code}")
                print(response.text)
        
        return results 