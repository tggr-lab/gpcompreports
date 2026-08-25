import pandas as pd

from GPCompaReports_v2.analysis import receptor_profile as rp


def _delta(rows):
    df = pd.DataFrame(rows, columns=['res1', 'res2', 'delta_rrcs'])
    df['abs_delta'] = df['delta_rrcs'].abs()
    return df


ANNOT = {
    1: {'position': 1, 'amino_acid': 'A', 'protein_segment': 'TM6', 'display_number': '6.48x48'},
    2: {'position': 2, 'amino_acid': 'C', 'protein_segment': 'TM3', 'display_number': '3.50x50'},
    3: {'position': 3, 'amino_acid': 'D', 'protein_segment': 'N-term', 'display_number': ''},
    4: {'position': 4, 'amino_acid': 'E', 'protein_segment': 'ECL3', 'display_number': ''},
}


def test_segment_profile_counts_both_endpoints_of_each_contact():
    df = _delta([(1, 2, 5.0), (3, 4, 9.0)])
    prof = rp.segment_profile(df, ANNOT, threshold=1.0)
    assert prof['TM6'] == 25.0
    assert prof['TM3'] == 25.0
    assert prof['N-term'] == 25.0
    assert prof['ECL3'] == 25.0


def test_segment_profile_ignores_below_threshold_contacts():
    df = _delta([(1, 2, 5.0), (3, 4, 0.5)])
    prof = rp.segment_profile(df, ANNOT, threshold=1.0)
    assert prof['TM6'] == 50.0
    assert prof['N-term'] == 0.0


def test_segment_profile_of_empty_frame_is_all_zero():
    prof = rp.segment_profile(_delta([]), ANNOT, threshold=1.0)
    assert set(prof) == set(rp.SEGMENTS)
    assert sum(prof.values()) == 0.0


def test_segment_profile_renormalizes_on_resolved_endpoints_only():
    # Position 5 has no annotation entry, so its endpoint doesn't resolve.
    # Endpoints: (1=TM6, 2=TM3) from the first row, (1=TM6, 5=unassigned)
    # from the second. Only 3 of the 4 endpoints resolve: TM6 x2, TM3 x1.
    df = _delta([(1, 2, 5.0), (1, 5, 9.0)])
    prof = rp.segment_profile(df, ANNOT, threshold=1.0)
    assert round(sum(prof.values()), 1) == 100.0
    # Under the old total-endpoints denominator this would be 50.0/25.0 —
    # the unresolved endpoint must not dilute the bars.
    assert prof['TM6'] == round(100.0 * 2 / 3, 1)
    assert prof['TM3'] == round(100.0 * 1 / 3, 1)


def test_segment_coverage_is_full_when_every_endpoint_resolves():
    df = _delta([(1, 2, 5.0), (3, 4, 9.0)])
    assert rp.segment_coverage(df, ANNOT, threshold=1.0) == 100.0


def test_segment_coverage_is_partial_when_some_endpoints_do_not_resolve():
    # 4 endpoints total, 3 resolve (position 5 has no annotation entry).
    df = _delta([(1, 2, 5.0), (1, 5, 9.0)])
    assert rp.segment_coverage(df, ANNOT, threshold=1.0) == 75.0


def test_segment_coverage_is_zero_with_no_above_threshold_contacts():
    assert rp.segment_coverage(_delta([]), ANNOT, threshold=1.0) == 0.0
    # Contacts exist but none clear the threshold.
    df = _delta([(1, 2, 0.1)])
    assert rp.segment_coverage(df, ANNOT, threshold=1.0) == 0.0


def test_median_profile_is_per_segment():
    a = {s: 0.0 for s in rp.SEGMENTS}
    b = {s: 0.0 for s in rp.SEGMENTS}
    c = {s: 0.0 for s in rp.SEGMENTS}
    a['TM6'], b['TM6'], c['TM6'] = 10.0, 20.0, 30.0
    assert rp.median_profile([a, b, c])['TM6'] == 20.0


def test_low_confidence_covers_termini_and_loops_but_not_helices():
    assert rp.is_low_confidence('N-term')
    assert rp.is_low_confidence('C-term')
    assert rp.is_low_confidence('ECL3')
    assert rp.is_low_confidence('ICL1')
    assert not rp.is_low_confidence('TM6')
    assert not rp.is_low_confidence('H8')


def test_largest_structured_skips_the_terminus_artifact():
    # The 9.0 contact is N-term to ECL3 and must not win.
    df = _delta([(1, 2, 5.0), (3, 4, 9.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={})
    assert kn['largest_structured']['abs_delta'] == 5.0
    assert kn['largest_structured']['seg1'] == 'TM6'


def test_largest_structured_is_none_when_every_contact_is_low_confidence():
    df = _delta([(3, 4, 9.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={})
    assert kn['largest_structured'] is None


def test_key_numbers_counts_cfr_top_movers_by_row_not_by_endpoint():
    # Both endpoints of this one row are CFRs. The row must still count once,
    # not twice: numerator and denominator are both row counts.
    df = _delta([(1, 2, 5.0)])
    kn = rp.key_numbers(
        df, ANNOT, threshold=1.0,
        cfr_ranks={'3.50x50': 4, '6.48x48': 7},
    )
    assert kn['cfr_top_movers'] == 1
    assert kn['top_mover_count'] == 1
    assert kn['above_threshold'] == 1
    assert kn['total_contacts'] == 1


def test_key_numbers_cfr_top_movers_never_exceeds_top_mover_count():
    # Six rows, every endpoint of every row is a CFR: with the old
    # by-endpoint counting this returned 9 (distinct numbers) against 6 rows.
    annot = {
        1: {'position': 1, 'amino_acid': 'A', 'protein_segment': 'TM6', 'display_number': '6.48x48'},
        2: {'position': 2, 'amino_acid': 'C', 'protein_segment': 'TM3', 'display_number': '3.50x50'},
        3: {'position': 3, 'amino_acid': 'D', 'protein_segment': 'TM2', 'display_number': '2.50x50'},
        4: {'position': 4, 'amino_acid': 'E', 'protein_segment': 'TM7', 'display_number': '7.53x53'},
    }
    ranks = {'6.48x48': 1, '3.50x50': 2, '2.50x50': 3, '7.53x53': 4}
    rows = [(1, 2, 9.0), (1, 2, 8.0), (1, 2, 7.0), (3, 4, 6.0), (3, 4, 5.0), (3, 4, 4.0)]
    df = _delta(rows)
    kn = rp.key_numbers(df, annot, threshold=1.0, cfr_ranks=ranks, top_n=6)
    assert kn['top_mover_count'] == 6
    assert kn['cfr_top_movers'] == 6
    assert kn['cfr_top_movers'] <= kn['top_mover_count']


def test_key_numbers_cfr_top_movers_excludes_rank_worse_than_50():
    # A CFR ranked worse than the cutoff must not count.
    df = _delta([(1, 2, 5.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={'6.48x48': 51})
    assert kn['cfr_top_movers'] == 0
    assert kn['top_mover_count'] == 1


def test_key_numbers_cfr_top_movers_includes_rank_exactly_50():
    # Rank 50 is the boundary itself ("50 or better"): it must count. This
    # is the one rank value that distinguishes `<=` from `<` in the cutoff
    # check, so it must be pinned directly rather than only tested at 51.
    df = _delta([(1, 2, 5.0)])
    kn = rp.key_numbers(df, ANNOT, threshold=1.0, cfr_ranks={'6.48x48': 50})
    assert kn['cfr_top_movers'] == 1
    assert kn['top_mover_count'] == 1
