"""Per-receptor derived numbers for the opt-in V3 analysis layer.

Pure functions over the frames the store already holds. Nothing here is
rendered unless the reader turns the analysis layer on, see
docs/superpowers/specs/2026-08-24-gpcompare-v3-design.md section 5.

GPCRdb numbering is read from `display_number`. The `generic_number` column
exists in every annotation CSV and is empty in all 283 of them.
"""

import statistics

SEGMENTS = [
    'N-term', 'TM1', 'ICL1', 'TM2', 'ECL1', 'TM3', 'ICL2', 'TM4',
    'ECL2', 'TM5', 'ICL3', 'TM6', 'ECL3', 'TM7', 'H8', 'C-term',
]

# AlphaFold-multistate is least reliable in termini and loops, and the single
# largest delta in a receptor usually lands there. Rows in these segments are
# marked low confidence and excluded from "largest structured change".
STRUCTURED = {'TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'H8'}

# A rank of 50 or better is the only bracket the CFR badge (and this count)
# claim to have verified; cfr_ranks holds ranks well past this on a full
# cfr_table, so anything worse than this cutoff does not count as a CFR.
CFR_RANK_CUTOFF = 50


def is_low_confidence(segment):
    """True for termini and loops, False for the seven helices and H8."""
    return segment not in STRUCTURED


def _seg(annot_map, position):
    entry = annot_map.get(int(position))
    if not entry:
        return 'unassigned'
    return entry.get('protein_segment') or 'unassigned'


def _num(annot_map, position):
    entry = annot_map.get(int(position))
    if not entry:
        return ''
    return entry.get('display_number') or entry.get('generic_number') or ''


def _aa(annot_map, position):
    entry = annot_map.get(int(position))
    if not entry:
        return ''
    return entry.get('amino_acid') or ''


def segment_profile(delta_df, annot_map, threshold):
    """Percent of above-threshold, annotated contact endpoints in each segment.

    Both ends of every above-threshold contact are counted, so a TM6-TM3
    contact contributes one to each. Returns a value for all 16 segments.

    Percentages are shares of endpoints that resolve to one of the 16 named
    segments only: an endpoint with no segment assignment is dropped from
    both the numerator and the denominator, not just the numerator. When at
    least one endpoint resolves, the 16 values sum to 100 within rounding;
    when none do, this returns all zeros. Because the base (annotated
    endpoints) differs from receptor to receptor, callers should report how
    much of the receptor that base actually covers via `segment_coverage`
    alongside these percentages, rather than presenting the bars alone.
    """
    profile = {s: 0.0 for s in SEGMENTS}
    if delta_df is None or delta_df.empty:
        return profile
    sig = delta_df[delta_df['abs_delta'] >= threshold]
    if sig.empty:
        return profile
    counts = {}
    accounted = 0
    for _, row in sig.iterrows():
        for col in ('res1', 'res2'):
            seg = _seg(annot_map, row[col])
            if seg not in profile:
                continue
            counts[seg] = counts.get(seg, 0) + 1
            accounted += 1
    if not accounted:
        return profile
    for seg in SEGMENTS:
        profile[seg] = round(100.0 * counts.get(seg, 0) / accounted, 1)
    return profile


def segment_coverage(delta_df, annot_map, threshold):
    """Percent of above-threshold contact endpoints that resolved to one of
    the 16 named segments, i.e. the base that `segment_profile` normalizes
    on. Returns 0.0 when there are no above-threshold contacts.
    """
    if delta_df is None or delta_df.empty:
        return 0.0
    sig = delta_df[delta_df['abs_delta'] >= threshold]
    if sig.empty:
        return 0.0
    total = 0
    accounted = 0
    for _, row in sig.iterrows():
        for col in ('res1', 'res2'):
            total += 1
            if _seg(annot_map, row[col]) in SEGMENTS:
                accounted += 1
    if not total:
        return 0.0
    return round(100.0 * accounted / total, 1)


def median_profile(profiles):
    """Per-segment median across many receptor profiles."""
    out = {}
    for seg in SEGMENTS:
        values = [p.get(seg, 0.0) for p in profiles]
        out[seg] = round(statistics.median(values), 1) if values else 0.0
    return out


def key_numbers(delta_df, annot_map, threshold, cfr_ranks, top_n=6):
    """Numbers for the analysis-layer strip.

    `cfr_ranks` maps a display_number to its cross-receptor CFR rank.
    `largest_structured` is the biggest |delta| whose two endpoints are both
    in a helix or H8, or None if the receptor has no such contact.
    `cfr_top_movers` counts top-mover ROWS (not endpoints), so it can never
    exceed `top_mover_count`: a row counts once if at least one of its two
    endpoints has a GPCRdb number ranked `CFR_RANK_CUTOFF` or better.
    """
    total = 0 if delta_df is None else len(delta_df)
    if delta_df is None or delta_df.empty:
        return {
            'total_contacts': 0, 'above_threshold': 0, 'above_threshold_pct': 0.0,
            'threshold': round(float(threshold), 2), 'top_segments': [],
            'largest_structured': None, 'cfr_top_movers': 0, 'top_mover_count': 0,
        }

    sig = delta_df[delta_df['abs_delta'] >= threshold]
    profile = segment_profile(delta_df, annot_map, threshold)
    ranked = sorted(profile.items(), key=lambda kv: -kv[1])
    top_segments = [(s, p) for s, p in ranked if p > 0][:3]

    largest = None
    for _, row in sig.sort_values('abs_delta', ascending=False).iterrows():
        s1, s2 = _seg(annot_map, row['res1']), _seg(annot_map, row['res2'])
        if s1 in STRUCTURED and s2 in STRUCTURED:
            largest = {
                'abs_delta': round(float(row['abs_delta']), 2),
                'delta': round(float(row['delta_rrcs']), 2),
                'seg1': s1, 'seg2': s2,
                'num1': _num(annot_map, row['res1']),
                'num2': _num(annot_map, row['res2']),
                'label1': '%s%d' % (_aa(annot_map, row['res1']), int(row['res1'])),
                'label2': '%s%d' % (_aa(annot_map, row['res2']), int(row['res2'])),
            }
            break

    movers = sig.sort_values('abs_delta', ascending=False).head(top_n)
    cfr_hits = 0
    for _, row in movers.iterrows():
        row_is_cfr = False
        for col in ('res1', 'res2'):
            num = _num(annot_map, row[col])
            rank = cfr_ranks.get(num) if num else None
            if rank is not None and rank <= CFR_RANK_CUTOFF:
                row_is_cfr = True
                break
        if row_is_cfr:
            cfr_hits += 1

    return {
        'total_contacts': total,
        'above_threshold': int(len(sig)),
        'above_threshold_pct': round(100.0 * len(sig) / total, 1) if total else 0.0,
        'threshold': round(float(threshold), 2),
        'top_segments': top_segments,
        'largest_structured': largest,
        'cfr_top_movers': cfr_hits,
        'top_mover_count': int(len(movers)),
    }
