"""Theme overrides for Plotly figures, keyed to v2 CSS tokens.

Figures are emitted once as Python/JSON (with the server defaults). The light
and dark layout overrides defined here are serialized into data-attrs on each
chart div, and applied by site_v2.js via Plotly.relayout() after the initial
newPlot(). This keeps analysis code theme-agnostic and avoids duplicating
figure data.

Only truly theme-dependent keys live here. Data-level encodings (RdBu colorscales,
brand teal/orange marks, semantic active/inactive colors) are theme-neutral
and stay untouched.

Dotted-key form is used because Plotly.relayout accepts it directly.
"""

FONT_FAMILY = (
    'system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", '
    'Arial, sans-serif'
)

_LIGHT = {
    'paper_bgcolor': '#ffffff',
    'plot_bgcolor':  '#ffffff',
    'font.color':    '#1f2328',
    'font.family':   FONT_FAMILY,
    'xaxis.gridcolor':      '#eaeef2',
    'xaxis.linecolor':      '#d1d9e0',
    'xaxis.tickcolor':      '#59636e',
    'xaxis.zerolinecolor':  '#d1d9e0',
    'xaxis.tickfont.color': '#59636e',
    'xaxis.title.font.color': '#1f2328',
    'yaxis.gridcolor':      '#eaeef2',
    'yaxis.linecolor':      '#d1d9e0',
    'yaxis.tickcolor':      '#59636e',
    'yaxis.zerolinecolor':  '#d1d9e0',
    'yaxis.tickfont.color': '#59636e',
    'yaxis.title.font.color': '#1f2328',
    'legend.font.color':    '#1f2328',
    'legend.bgcolor':       'rgba(255,255,255,0)',
    'legend.bordercolor':   'rgba(0,0,0,0)',
    'title.font.color':     '#1f2328',
}

_DARK = {
    'paper_bgcolor': '#0f1417',
    'plot_bgcolor':  '#0f1417',
    'font.color':    '#e2e8f0',
    'font.family':   FONT_FAMILY,
    'xaxis.gridcolor':      '#1f2937',
    'xaxis.linecolor':      '#334155',
    'xaxis.tickcolor':      '#94a3b8',
    'xaxis.zerolinecolor':  '#334155',
    'xaxis.tickfont.color': '#94a3b8',
    'xaxis.title.font.color': '#e2e8f0',
    'yaxis.gridcolor':      '#1f2937',
    'yaxis.linecolor':      '#334155',
    'yaxis.tickcolor':      '#94a3b8',
    'yaxis.zerolinecolor':  '#334155',
    'yaxis.tickfont.color': '#94a3b8',
    'yaxis.title.font.color': '#e2e8f0',
    'legend.font.color':    '#e2e8f0',
    'legend.bgcolor':       'rgba(15,20,23,0)',
    'legend.bordercolor':   'rgba(255,255,255,0)',
    'title.font.color':     '#e2e8f0',
}


def theme_overrides():
    """Return (light_layout, dark_layout) dicts. Each is a dotted-key mapping
    suitable for Plotly.relayout(). Copies so callers can safely mutate."""
    return dict(_LIGHT), dict(_DARK)
