"""Load all batch RRCS analysis data into memory for cross-GPCR analysis."""

import ast
import os
from pathlib import Path

import pandas as pd


# Default batch results path
BATCH_DIR = Path(__file__).resolve().parent.parent.parent / (
    "The_batch_RRCS_analyzer/batch_analysis_full/batch_analysis_20260202_151051"
)
METADATA_CSV = Path(__file__).resolve().parent.parent.parent / "class_A_all.csv"


class GPCRDataStore:
    """Central data store holding all loaded GPCR data."""

    def __init__(self, batch_dir=None, metadata_csv=None):
        self.batch_dir = Path(batch_dir) if batch_dir else BATCH_DIR
        self.metadata_csv = Path(metadata_csv) if metadata_csv else METADATA_CSV
        self.csv_dir = self.batch_dir / "csv_data"
        self.variant_dir = self.batch_dir / "variants"
        self.summary_csv = self.batch_dir / "summary" / "processing_results.csv"

        # Data containers
        self.metadata = None          # DataFrame: uniprot_name, gene_name, family, ligand_type
        self.processing_summary = None  # DataFrame from processing_results.csv
        self.gpcr_ids = []            # list of filesystem gpcr ids (e.g. '5ht1a_human')
        self.delta_data = {}          # gpcr_id -> DataFrame
        self.annotation_data = {}     # gpcr_id -> DataFrame
        self.significant_data = {}    # gpcr_id -> DataFrame
        self.variant_data = {}        # gpcr_id -> DataFrame
        self.gpcr_info = {}           # gpcr_id -> dict with metadata + summary stats
        self.name_map = {}            # filesystem id -> uniprot name (uppercase)

    def load_all(self):
        """Load all data sources."""
        print("Loading GPCR metadata...")
        self._load_metadata()
        print("Loading processing summary...")
        self._load_processing_summary()
        print("Building name map...")
        self._build_name_map()
        print(f"Loading CSV data for {len(self.gpcr_ids)} GPCRs...")
        self._load_csv_data()
        print(f"Loading variant data...")
        self._load_variant_data()
        print("Building GPCR info index...")
        self._build_gpcr_info()
        print(f"Data loading complete: {len(self.gpcr_ids)} GPCRs loaded")
        return self

    def _load_metadata(self):
        """Load class_A_all.csv with multiline header handling."""
        raw = self.metadata_csv.read_text()
        lines = raw.strip().split('\n')
        # The header spans lines 1-3 due to quoted multiline fields
        # Reconstruct: skip multiline header, use clean column names
        data_lines = []
        i = 0
        in_quotes = False
        header_end = 0
        for idx, line in enumerate(lines):
            if idx == 0:
                in_quotes = line.count('"') % 2 == 1
                if not in_quotes:
                    header_end = idx
                continue
            if in_quotes:
                in_quotes = line.count('"') % 2 == 0
                if not in_quotes:
                    header_end = idx
                continue
            data_lines.append(line)

        self.metadata = pd.DataFrame(
            [self._parse_csv_line(l) for l in data_lines if l.strip()],
            columns=['uniprot_name', 'gene_name', 'receptor_family', 'ligand_type', 'class']
        )

    @staticmethod
    def _parse_csv_line(line):
        """Parse a CSV line handling quoted fields with commas."""
        fields = []
        current = ''
        in_quotes = False
        for ch in line:
            if ch == '"':
                in_quotes = not in_quotes
            elif ch == ',' and not in_quotes:
                fields.append(current.strip())
                current = ''
            else:
                current += ch
        fields.append(current.strip())
        return fields

    def _load_processing_summary(self):
        """Load processing_results.csv and parse summary dicts."""
        self.processing_summary = pd.read_csv(self.summary_csv)
        # Parse the summary dict strings
        summaries = []
        for _, row in self.processing_summary.iterrows():
            try:
                s = ast.literal_eval(row['summary'])
            except (ValueError, SyntaxError):
                s = {}
            summaries.append(s)
        self.processing_summary['parsed_summary'] = summaries

    def _build_name_map(self):
        """Map filesystem IDs to metadata entries.

        Filesystem: '5ht1a_human' -> metadata: '5HT1A'
        Mapping: lowercase(uniprot_name) + '_human' == filesystem_id
        """
        # Discover all GPCRs with delta CSV files
        delta_files = sorted(self.csv_dir.glob("*_rrcs_delta.csv"))
        self.gpcr_ids = [f.stem.replace("_rrcs_delta", "") for f in delta_files]

        # Build forward and reverse maps
        meta_lookup = {}
        for _, row in self.metadata.iterrows():
            key = row['uniprot_name'].lower() + '_human'
            meta_lookup[key] = row['uniprot_name']

        for gid in self.gpcr_ids:
            if gid in meta_lookup:
                self.name_map[gid] = meta_lookup[gid]
            else:
                # Fallback: try matching via processing_summary uniprot_id
                self.name_map[gid] = gid.replace('_human', '').upper()

    def _load_csv_data(self):
        """Load delta, annotation, and significant_changes CSVs."""
        for gid in self.gpcr_ids:
            delta_path = self.csv_dir / f"{gid}_rrcs_delta.csv"
            annot_path = self.csv_dir / f"{gid}_annotations.csv"
            sig_path = self.csv_dir / f"{gid}_significant_changes.csv"

            if delta_path.exists():
                self.delta_data[gid] = pd.read_csv(delta_path)
            if annot_path.exists():
                self.annotation_data[gid] = pd.read_csv(annot_path)
            if sig_path.exists():
                self.significant_data[gid] = pd.read_csv(sig_path)

    def _load_variant_data(self):
        """Load variant CSVs."""
        for gid in self.gpcr_ids:
            var_path = self.variant_dir / f"{gid}_variants.csv"
            if var_path.exists():
                self.variant_data[gid] = pd.read_csv(var_path)

    def _build_gpcr_info(self):
        """Build per-GPCR info dicts combining all data sources."""
        meta_dict = {}
        for _, row in self.metadata.iterrows():
            meta_dict[row['uniprot_name']] = {
                'gene_name': row['gene_name'],
                'receptor_family': row['receptor_family'],
                'ligand_type': row['ligand_type'],
            }

        # Parse processing summary into lookup
        proc_dict = {}
        for _, row in self.processing_summary.iterrows():
            proc_dict[row['gpcr']] = row['parsed_summary']

        for gid in self.gpcr_ids:
            uname = self.name_map.get(gid, gid.replace('_human', '').upper())
            meta = meta_dict.get(uname, {})
            proc = proc_dict.get(gid, {})

            delta_df = self.delta_data.get(gid, pd.DataFrame())
            sig_df = self.significant_data.get(gid, pd.DataFrame())
            var_df = self.variant_data.get(gid, pd.DataFrame())

            info = {
                'gpcr_id': gid,
                'uniprot_name': uname,
                'gene_name': meta.get('gene_name', ''),
                'receptor_family': meta.get('receptor_family', ''),
                'ligand_type': meta.get('ligand_type', ''),
                'uniprot_id': proc.get('uniprot_id', ''),
                'total_contacts': len(delta_df),
                'significant_changes': len(sig_df),
                'variants_found': len(var_df),
                'tm_domains_found': proc.get('tm_domains_found', 0),
                'sum_abs_delta': float(delta_df['abs_delta'].sum()) if len(delta_df) > 0 else 0.0,
                'mean_abs_delta': float(delta_df['abs_delta'].mean()) if len(delta_df) > 0 else 0.0,
                'max_abs_delta': float(delta_df['abs_delta'].max()) if len(delta_df) > 0 else 0.0,
            }
            self.gpcr_info[gid] = info

    def get_all_info_df(self):
        """Return a DataFrame of all GPCR info for ranking/comparison."""
        return pd.DataFrame(list(self.gpcr_info.values()))

    def get_annotation_map(self, gpcr_id):
        """Return position -> annotation dict for a GPCR."""
        df = self.annotation_data.get(gpcr_id, pd.DataFrame())
        if df.empty:
            return {}
        result = {}
        for _, row in df.iterrows():
            pos = int(row['position'])
            result[pos] = {
                'amino_acid': _safe_str(row.get('amino_acid', '')),
                'generic_number': _safe_str(row.get('generic_number', '')),
                'display_number': _safe_str(row.get('display_number', '')),
                'protein_segment': _safe_str(row.get('protein_segment', '')),
            }
        return result


def _safe_str(val):
    """Convert value to string, returning '' for NaN/None."""
    if pd.isna(val):
        return ''
    return str(val).strip()
