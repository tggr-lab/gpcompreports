"""
Snake Plot Data Builder
Prepares multi-view color maps, contact links, and metadata for the
interactive snake plot embedded in RRCS HTML reports.
"""

import json
import pandas as pd
from typing import Dict, List, Optional, Tuple


# ── Amino acid property classification (mirrors ouroboros/coloring.py) ────

_HYDROPHOBIC = set("AILMFWVP")
_POLAR = set("STNQ")
_POSITIVE = set("RHK")
_NEGATIVE = set("DE")
_AROMATIC = set("FWY")
_SPECIAL = set("CGP")

PROPERTY_COLORS = {
    "hydrophobic": "#F0C05A",  # gold
    "polar": "#87CEEB",        # sky blue
    "positive": "#FF7F7F",     # coral
    "negative": "#5F9EA0",     # teal
    "aromatic": "#DDA0DD",     # plum
    "special": "#98FB98",      # mint
    "unknown": "#FFFFFF",
}


def _aa_property(aa: str) -> str:
    """Classify a single-letter amino acid code by physicochemical property."""
    c = aa.upper()
    if c in _AROMATIC:
        return "aromatic"
    if c in _HYDROPHOBIC:
        return "hydrophobic"
    if c in _POSITIVE:
        return "positive"
    if c in _NEGATIVE:
        return "negative"
    if c in _POLAR:
        return "polar"
    if c in _SPECIAL:
        return "special"
    return "unknown"


def _interpolate(value: float, low: str, high: str) -> str:
    """Linearly interpolate between two hex colors for value in [0, 1]."""
    value = max(0.0, min(1.0, value))
    r1, g1, b1 = int(low[1:3], 16), int(low[3:5], 16), int(low[5:7], 16)
    r2, g2, b2 = int(high[1:3], 16), int(high[3:5], 16), int(high[5:7], 16)
    r = int(r1 + (r2 - r1) * value)
    g = int(g1 + (g2 - g1) * value)
    b = int(b1 + (b2 - b1) * value)
    return f"#{r:02x}{g:02x}{b:02x}"


def _diverging_color(value: float, min_val: float, max_val: float,
                     neg_color: str, pos_color: str,
                     mid_color: str = "#FFFFFF") -> str:
    """Map a value to a diverging color scale symmetric around zero."""
    abs_max = max(abs(min_val), abs(max_val))
    if abs_max == 0:
        return mid_color
    normalized = value / abs_max  # range [-1, 1]
    if normalized >= 0:
        return _interpolate(normalized, mid_color, pos_color)
    else:
        return _interpolate(-normalized, mid_color, neg_color)


class SnakePlotDataBuilder:
    """Build multi-view color data for the interactive snake plot.

    Produces 8 view color maps, contact link data, and legend metadata,
    all serializable to JSON for embedding in an HTML report.
    """

    def __init__(
        self,
        delta_matrix: pd.DataFrame,
        variants_df: pd.DataFrame,
        protein_data,  # ouroboros ProteinData
        significance_threshold: float = 1.0,
        max_links: int = 50,
    ):
        self.delta_matrix = delta_matrix
        self.variants_df = variants_df
        self.protein_data = protein_data
        self.significance_threshold = significance_threshold
        self.max_links = max_links

        # Build residue lookups
        self.residue_map: Dict[int, str] = {}  # seq → one-letter AA
        self.segment_map: Dict[int, str] = {}  # seq → protein_segment
        _helix_segments = {'TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8'}
        self.helix_residues: set = set()  # residues always visible in collapsed view
        for r in protein_data.residues:
            self.residue_map[r.sequence_number] = r.one_letter
            if r.protein_segment:
                self.segment_map[r.sequence_number] = r.protein_segment
                if r.protein_segment in _helix_segments:
                    self.helix_residues.add(r.sequence_number)

        # Pre-compute per-residue aggregates from delta matrix
        self._compute_residue_aggregates()

    def _compute_residue_aggregates(self):
        """Compute per-residue mean delta, sum active, and sum inactive."""
        dm = self.delta_matrix
        all_residues = set(dm['res1'].tolist()) | set(dm['res2'].tolist())

        self.residue_mean_delta: Dict[int, float] = {}
        self.residue_sum_active: Dict[int, float] = {}
        self.residue_sum_inactive: Dict[int, float] = {}

        for res in all_residues:
            mask = (dm['res1'] == res) | (dm['res2'] == res)
            subset = dm[mask]
            self.residue_mean_delta[res] = float(subset['delta_rrcs'].mean())
            self.residue_sum_active[res] = float(subset['active_rrcs'].sum())
            self.residue_sum_inactive[res] = float(subset['inactive_rrcs'].sum())

    # ── Individual view builders ─────────────────────────────────────

    def _build_delta_rrcs_view(self) -> Dict[str, str]:
        """Delta RRCS: red (inactive-favoring) ← white → blue (active-favoring)."""
        colors: Dict[str, str] = {}
        if not self.residue_mean_delta:
            return colors

        values = list(self.residue_mean_delta.values())
        min_val = min(values)
        max_val = max(values)

        for seq, val in self.residue_mean_delta.items():
            colors[str(seq)] = _diverging_color(
                val, min_val, max_val,
                neg_color="#B2182B",  # dark red  = inactive-favoring
                pos_color="#2166AC",  # dark blue = active-favoring
            )
        return colors

    def _build_active_rrcs_view(self) -> Dict[str, str]:
        """Active RRCS: white → dark blue."""
        colors: Dict[str, str] = {}
        if not self.residue_sum_active:
            return colors

        max_val = max(self.residue_sum_active.values()) or 1.0
        for seq, val in self.residue_sum_active.items():
            colors[str(seq)] = _interpolate(val / max_val, "#FFFFFF", "#2166AC")
        return colors

    def _build_inactive_rrcs_view(self) -> Dict[str, str]:
        """Inactive RRCS: white → dark red."""
        colors: Dict[str, str] = {}
        if not self.residue_sum_inactive:
            return colors

        max_val = max(self.residue_sum_inactive.values()) or 1.0
        for seq, val in self.residue_sum_inactive.items():
            colors[str(seq)] = _interpolate(val / max_val, "#FFFFFF", "#B2182B")
        return colors

    def _build_variants_view(self) -> Dict[str, str]:
        """Variant positions: colored by variant count (rarity). All positions with
        any variant are shown, independent of AlphaMissense classification."""
        colors: Dict[str, str] = {}
        if self.variants_df is None or self.variants_df.empty:
            return colors
        if 'protein_position' not in self.variants_df.columns:
            return colors

        # Count distinct variants per position
        counts = self.variants_df.groupby('protein_position').size()
        max_count = counts.max() if len(counts) > 0 else 1

        for pos, count in counts.items():
            # Yellow (1 variant) → deep orange (many variants)
            normalized = min(count / max(max_count, 1), 1.0)
            colors[str(int(pos))] = _interpolate(normalized, "#FFF176", "#E65100")

        return colors

    def _build_conservation_view(self) -> Dict[str, str]:
        """Conservation score: white → deep purple (0 → 1)."""
        colors: Dict[str, str] = {}
        if self.variants_df is None or self.variants_df.empty:
            return colors
        if 'conservation' not in self.variants_df.columns:
            return colors

        for pos, group in self.variants_df.groupby('protein_position'):
            cons_values = group['conservation'].dropna()
            if not cons_values.empty:
                score = float(cons_values.iloc[0])
                colors[str(int(pos))] = _interpolate(score, "#FFFFFF", "#4A148C")

        return colors

    def _build_alphamissense_view(self) -> Dict[str, str]:
        """AlphaMissense pathogenicity: green → red (0 → 1)."""
        colors: Dict[str, str] = {}
        if self.variants_df is None or self.variants_df.empty:
            return colors
        if 'am_score' not in self.variants_df.columns:
            return colors

        for pos, group in self.variants_df.groupby('protein_position'):
            am_values = group['am_score'].dropna()
            if not am_values.empty:
                score = float(am_values.max())
                colors[str(int(pos))] = _interpolate(score, "#388E3C", "#D32F2F")

        return colors

    def _build_cfr_view(self) -> Dict[str, str]:
        """Core functional residues — placeholder (empty until CFR data available)."""
        return {}

    def _build_aa_properties_view(self) -> Dict[str, str]:
        """Amino acid physicochemical property colors for all residues."""
        colors: Dict[str, str] = {}
        for seq, aa in self.residue_map.items():
            prop = _aa_property(aa)
            color = PROPERTY_COLORS.get(prop, "#FFFFFF")
            if color != "#FFFFFF":
                colors[str(seq)] = color
        return colors

    # ── Contact links ────────────────────────────────────────────────

    def _build_contact_links(self) -> List[Dict]:
        """Build contact link data for the most significant residue pairs.

        Returns list of {from_seq, to_seq, color, width, delta} dicts.
        Color: blue = active-favoring (positive delta), red = inactive-favoring.
        """
        dm = self.delta_matrix
        significant = dm[dm['abs_delta'] >= self.significance_threshold].copy()
        if significant.empty:
            return []

        # Only keep contacts where both residues are in helix segments
        # (TM1-7, H8) — loop/terminus residues are hidden by default
        helix = self.helix_residues
        significant = significant[
            significant['res1'].isin(helix) & significant['res2'].isin(helix)
        ]
        if significant.empty:
            return []

        significant = significant.nlargest(self.max_links, 'abs_delta')

        max_abs = significant['abs_delta'].max()
        min_abs = significant['abs_delta'].min()

        links = []
        for _, row in significant.iterrows():
            delta = float(row['delta_rrcs'])
            abs_d = float(row['abs_delta'])

            # Blue = active-favoring (positive), red = inactive-favoring (negative)
            color = "#2166AC" if delta > 0 else "#B2182B"

            # Width: 1–4 px proportional to magnitude
            if max_abs > min_abs:
                width = 1.0 + 3.0 * (abs_d - min_abs) / (max_abs - min_abs)
            else:
                width = 2.0

            links.append({
                'from_seq': int(row['res1']),
                'to_seq': int(row['res2']),
                'color': color,
                'width': round(width, 1),
                'delta': round(delta, 2),
            })

        return links

    # ── Public API ───────────────────────────────────────────────────

    def build(self, positions: Dict[int, Tuple[float, float]]) -> Dict:
        """Build all view data and return as a serializable dict.

        Args:
            positions: {seq_number: (cx, cy)} extracted from the base SVG.
        """
        views = {
            'delta_rrcs': {
                'colors': self._build_delta_rrcs_view(),
                'label': 'Delta RRCS',
                'description': 'Per-residue mean \u0394RRCS (active \u2212 inactive)',
                'legend_type': 'diverging',
                'legend_colors': ['#B2182B', '#FFFFFF', '#2166AC'],
                'legend_labels': ['Inactive-favoring', '0', 'Active-favoring'],
            },
            'active_rrcs': {
                'colors': self._build_active_rrcs_view(),
                'label': 'Active RRCS',
                'description': 'Total active-state contact strength per residue',
                'legend_type': 'sequential',
                'legend_colors': ['#FFFFFF', '#2166AC'],
                'legend_labels': ['Low', 'High'],
            },
            'inactive_rrcs': {
                'colors': self._build_inactive_rrcs_view(),
                'label': 'Inactive RRCS',
                'description': 'Total inactive-state contact strength per residue',
                'legend_type': 'sequential',
                'legend_colors': ['#FFFFFF', '#B2182B'],
                'legend_labels': ['Low', 'High'],
            },
            'variants': {
                'colors': self._build_variants_view(),
                'label': 'Variants',
                'description': 'Positions with known gnomAD variants (colored by variant count)',
                'legend_type': 'sequential',
                'legend_colors': ['#FFF176', '#E65100'],
                'legend_labels': ['1 variant', 'Many variants'],
            },
            'conservation': {
                'colors': self._build_conservation_view(),
                'label': 'Conservation',
                'description': 'Sequence conservation score (0\u20131)',
                'legend_type': 'sequential',
                'legend_colors': ['#FFFFFF', '#4A148C'],
                'legend_labels': ['Variable', 'Conserved'],
            },
            'alphamissense': {
                'colors': self._build_alphamissense_view(),
                'label': 'AlphaMissense',
                'description': 'AlphaMissense pathogenicity score (0\u20131)',
                'legend_type': 'sequential',
                'legend_colors': ['#388E3C', '#D32F2F'],
                'legend_labels': ['Benign', 'Pathogenic'],
            },
            'cfr': {
                'colors': self._build_cfr_view(),
                'label': 'Core Functional',
                'description': 'Core functional residues (data not yet available)',
                'legend_type': 'categorical',
                'legend_items': [
                    {'color': '#FF8F00', 'label': 'CFR'},
                ],
            },
            'aa_properties': {
                'colors': self._build_aa_properties_view(),
                'label': 'AA Properties',
                'description': 'Amino acid physicochemical properties',
                'legend_type': 'categorical',
                'legend_items': [
                    {'color': '#F0C05A', 'label': 'Hydrophobic'},
                    {'color': '#87CEEB', 'label': 'Polar'},
                    {'color': '#FF7F7F', 'label': 'Positive'},
                    {'color': '#5F9EA0', 'label': 'Negative'},
                    {'color': '#DDA0DD', 'label': 'Aromatic'},
                    {'color': '#98FB98', 'label': 'Special (C,G,P)'},
                ],
            },
        }

        positions_json = {str(k): [v[0], v[1]] for k, v in positions.items()}

        return {
            'views': views,
            'contact_links': self._build_contact_links(),
            'positions': positions_json,
        }

    def to_json(self, positions: Dict[int, Tuple[float, float]]) -> str:
        """Build all data and serialize to a JSON string."""
        return json.dumps(self.build(positions))
