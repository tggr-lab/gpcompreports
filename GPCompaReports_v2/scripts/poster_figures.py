#!/usr/bin/env python3
"""Generate standalone poster figures for the class A RRCS project.

Loads the full data store once, reuses the same computations the website runs
(so the numbers match the live statistics page), and writes six high-resolution
PNG figures restyled for print: large fonts, light background, the project
red/blue convention (red = active-favoring, blue = inactive-favoring), and
plain-language labels.

Run from the project root:
    python3 -m GPCompaReports_v2.scripts.poster_figures
or directly:
    python3 GPCompaReports_v2/scripts/poster_figures.py

Output: GPCompaReports_v2/output/poster/*.png
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2_ROOT = HERE.parent
PROJECT_ROOT = V2_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from GPCompaReports_v2.analysis.data_loader import GPCRDataStore
from GPCompaReports_v2.analysis.cfr_analysis import identify_cfrs, build_cfr_network
from GPCompaReports_v2.analysis.variant_correlation import (
    _get_cfr_position_map, pathogenicity_enrichment)
from GPCompaReports_v2.analysis.tm_domain_analysis import compute_domain_interactions

OUT_DIR = V2_ROOT / "output" / "poster"

# Project color convention: red = active-favoring, blue = inactive-favoring.
ACTIVE = "#B91C1C"     # red
INACTIVE = "#1D4ED8"   # blue
TEAL = "#0D9488"       # brand accent for single-series bars
FONT = "Arial, Helvetica, sans-serif"

# Pathogenicity tier colors (kept distinct from the active/inactive semantic so
# they are not read as conformational state).
TIER_COLORS = {
    "Likely pathogenic": "#E8820C",
    "Uncertain": "#F0C27A",
    "Likely benign": "#008080",
}


def short_gn(gn: str) -> str:
    """Format a generic number for a general audience: '3.50x50' -> '3.50'."""
    return str(gn).split("x")[0]


def fmt_p(p: float) -> str:
    """Readable p-value text, no raw scientific floats."""
    if p is None:
        return ""
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3f}"


def apply_poster_theme(fig, *, width=1100, height=800, show_legend=None):
    fig.update_layout(
        template="simple_white",
        font=dict(family=FONT, size=18, color="#1a1a1a"),
        title=dict(font=dict(size=26), x=0.5, xanchor="center"),
        width=width, height=height,
        margin=dict(l=120, r=60, t=110, b=90),
        paper_bgcolor="white", plot_bgcolor="white",
    )
    if show_legend is not None:
        fig.update_layout(showlegend=show_legend)
    fig.update_xaxes(title_font_size=20, tickfont_size=16)
    fig.update_yaxes(title_font_size=20, tickfont_size=16)
    return fig


def save(fig, name):
    path = OUT_DIR / f"{name}.png"
    fig.write_image(str(path), scale=3)
    print(f"  wrote {path}")


# --------------------------------------------------------------------------
# F1: Recurrent contacts between core functional residues
# --------------------------------------------------------------------------
def fig_cfr_contacts(cfr_network):
    top = cfr_network.head(20).copy()
    top["label"] = top.apply(
        lambda r: f"{short_gn(r['cfr_1'])} to {short_gn(r['cfr_2'])}", axis=1)
    top = top.sort_values("n_gpcrs", ascending=True)

    fig = go.Figure(go.Bar(
        x=top["n_gpcrs"], y=top["label"], orientation="h",
        marker_color=TEAL,
        text=top["n_gpcrs"], textposition="outside",
        cliponaxis=False,
    ))
    fig.update_layout(
        title="Recurrent contacts between core functional residues",
        xaxis_title="Number of receptors (of 283)",
        yaxis_title="Contact (residue pair)",
    )
    apply_poster_theme(fig, width=1100, height=900, show_legend=False)
    fig.update_layout(margin=dict(l=180, r=80, t=110, b=90))
    return fig


# --------------------------------------------------------------------------
# F2: AlphaMissense pathogenicity at core vs other positions
# --------------------------------------------------------------------------
def fig_pathogenicity(path_result):
    groups = ["Core functional residues", "Other positions"]
    tiers = [
        ("Likely pathogenic", "pathogenic"),
        ("Uncertain", "ambiguous"),
        ("Likely benign", "benign"),
    ]
    cfr, non = path_result["cfr_counts"], path_result["non_cfr_counts"]
    cfr_tot = path_result["cfr_total"] or 1
    non_tot = path_result["non_cfr_total"] or 1

    fig = go.Figure()
    for label, key in tiers:
        pcts = [cfr[key] / cfr_tot * 100, non[key] / non_tot * 100]
        fig.add_trace(go.Bar(
            name=label, x=groups, y=pcts,
            marker_color=TIER_COLORS[label],
            text=[f"{v:.0f}%" for v in pcts], textposition="outside",
            cliponaxis=False,
        ))

    p = path_result.get("stats", {}).get("p_value")
    contrast = (f"{path_result['cfr_pathogenic_pct']:.0f}% of variants at core "
                f"positions are likely pathogenic, versus "
                f"{path_result['non_cfr_pathogenic_pct']:.0f}% elsewhere")
    fig.add_annotation(
        x=0.5, y=1.13, xref="paper", yref="paper", showarrow=False,
        text=contrast + (f"  ({fmt_p(p)})" if p is not None else ""),
        font=dict(size=17, color="#444"),
    )

    fig.update_layout(
        title="Predicted variant effects at core functional residues",
        yaxis_title="Share of variants (%)",
        barmode="group", bargap=0.35, bargroupgap=0.1,
        legend=dict(orientation="h", y=-0.12, x=0.5, xanchor="center"),
        yaxis_range=[0, max(60, path_result["cfr_pathogenic_pct"] + 12)],
    )
    apply_poster_theme(fig, width=1000, height=820, show_legend=True)
    fig.update_layout(margin=dict(l=110, r=60, t=150, b=120))
    return fig


# --------------------------------------------------------------------------
# F3: Core functional residue hotspots (frequency vs magnitude)
# --------------------------------------------------------------------------
def fig_cfr_hotspots(cfr_table):
    # Keep recurrent positions only: a core functional residue should appear in
    # a meaningful share of receptors. This drops rare high-magnitude outliers
    # that would otherwise stretch the axis and squash the main cluster.
    df = cfr_table[cfr_table["frequency_pct"] >= 10].head(40).copy()
    df["pos"] = df["generic_number"].map(short_gn)
    highlight = {"3.50", "7.53", "6.48", "7.54", "8.50"}

    fig = px.scatter(
        df, x="mean_abs_delta", y="frequency_pct",
        size="cfr_score", color="segment",
        size_max=46,
        labels={
            "mean_abs_delta": "Average contact-score change",
            "frequency_pct": "Share of receptors affected (%)",
            "segment": "Helix or loop",
        },
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_traces(marker=dict(line=dict(width=1, color="white")), opacity=0.9)

    for _, r in df[df["pos"].isin(highlight)].iterrows():
        fig.add_annotation(
            x=r["mean_abs_delta"], y=r["frequency_pct"], text=r["pos"],
            showarrow=True, arrowhead=0, arrowwidth=1, arrowcolor="#888",
            ax=0, ay=-28, font=dict(size=16, color="#1a1a1a"),
        )

    fig.update_layout(
        title="Core functional residues across class A GPCRs",
        legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center"),
    )
    apply_poster_theme(fig, width=1100, height=850, show_legend=True)
    fig.update_layout(margin=dict(l=120, r=60, t=110, b=150))
    return fig


# --------------------------------------------------------------------------
# F4: Helix-pair rearrangement heatmap (7x7 TM interfaces)
# --------------------------------------------------------------------------
def fig_tm_heatmap(tm_matrix):
    z = tm_matrix.values
    fig = go.Figure(go.Heatmap(
        z=z, x=tm_matrix.columns.tolist(), y=tm_matrix.index.tolist(),
        colorscale=[[0, "#ffffff"], [1, TEAL]],
        text=np.round(z, 1), texttemplate="%{text}", textfont=dict(size=16),
        colorbar=dict(title=dict(text="Average<br>change", font=dict(size=16)),
                      tickfont=dict(size=15)),
        hovertemplate="%{y} and %{x}<br>Average change: %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title="Which helix pairs rearrange between states",
        yaxis_autorange="reversed",
        xaxis_title="", yaxis_title="",
    )
    apply_poster_theme(fig, width=900, height=820, show_legend=False)
    fig.update_layout(margin=dict(l=90, r=60, t=110, b=80))
    return fig


# --------------------------------------------------------------------------
# F5: Conformational change by ligand type (box plot)
# --------------------------------------------------------------------------
def fig_by_ligand(info_df):
    df = info_df[info_df["ligand_type"].astype(str).str.strip() != ""].copy()
    order = (df.groupby("ligand_type")["sum_abs_delta"].median()
               .sort_values(ascending=False).index.tolist())

    fig = px.box(
        df, x="ligand_type", y="sum_abs_delta", color="ligand_type",
        category_orders={"ligand_type": order},
        labels={
            "ligand_type": "Ligand type",
            "sum_abs_delta": "Total contact-score change per receptor",
        },
        color_discrete_sequence=px.colors.qualitative.Set2,
    )
    fig.update_layout(
        title="Conformational change by ligand type",
        xaxis_tickangle=-35,
    )
    apply_poster_theme(fig, width=1200, height=820, show_legend=False)
    fig.update_layout(margin=dict(l=120, r=60, t=110, b=190))
    return fig


# --------------------------------------------------------------------------
# F6: Distribution of contact-score changes across all receptors
# --------------------------------------------------------------------------
def fig_change_distribution(store):
    deltas = []
    for gid in store.gpcr_ids:
        df = store.delta_data.get(gid, pd.DataFrame())
        if not df.empty:
            deltas.extend(df["delta_rrcs"].tolist())
    arr = np.asarray(deltas, dtype=float)
    arr = arr[arr != 0]

    lo, hi = np.floor(arr.min()), np.ceil(arr.max())
    bins = dict(start=float(lo), end=float(hi), size=0.5)

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=arr[arr < 0], xbins=bins, marker_color=INACTIVE, opacity=0.85,
        name="Favors inactive state",
    ))
    fig.add_trace(go.Histogram(
        x=arr[arr > 0], xbins=bins, marker_color=ACTIVE, opacity=0.85,
        name="Favors active state",
    ))
    fig.update_layout(
        title="Distribution of contact-score changes",
        xaxis_title="Contact-score change (active vs inactive)",
        yaxis_title="Number of contacts",
        barmode="overlay", bargap=0,
        legend=dict(orientation="h", y=-0.14, x=0.5, xanchor="center"),
    )
    apply_poster_theme(fig, width=1150, height=780, show_legend=True)
    fig.update_layout(margin=dict(l=130, r=60, t=110, b=120))
    # The bulk of changes sit within a narrow window; crop the sparse tails so
    # the bars stay legible.
    fig.update_xaxes(range=[-12, 12])
    return fig


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("Loading data store (283 receptors)...")
    store = GPCRDataStore().load_all()

    print("Computing core functional residues...")
    cfr_table = identify_cfrs(store)
    cfr_network = build_cfr_network(store, cfr_table)

    print("Computing variant pathogenicity enrichment...")
    cfr_pos = _get_cfr_position_map(store, cfr_table)
    path_result = pathogenicity_enrichment(store, cfr_pos)

    print("Computing helix interface matrix...")
    tm_matrix = compute_domain_interactions(store)

    info_df = store.get_all_info_df()

    print("Rendering figures...")
    save(fig_cfr_contacts(cfr_network), "poster_cfr_contacts")
    save(fig_pathogenicity(path_result), "poster_pathogenicity")
    save(fig_cfr_hotspots(cfr_table), "poster_cfr_hotspots")
    save(fig_tm_heatmap(tm_matrix), "poster_tm_heatmap")
    save(fig_by_ligand(info_df), "poster_by_ligand")
    save(fig_change_distribution(store), "poster_change_distribution")
    print(f"\nDone. Six figures in {OUT_DIR}")


if __name__ == "__main__":
    main()
