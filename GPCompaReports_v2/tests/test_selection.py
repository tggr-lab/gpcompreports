import pytest

from GPCompaReports_v2.website.site_generator import select_gpcr_ids

ALL = ['5ht2a_human', 'adrb2_human', 'cxcr4_human', 'par1_human', 'par2_human']


def test_no_filters_returns_everything():
    assert select_gpcr_ids(ALL) == ALL


def test_limit_takes_a_prefix():
    assert select_gpcr_ids(ALL, limit=2) == ['5ht2a_human', 'adrb2_human']


def test_only_selects_named_receptors_in_given_order():
    assert select_gpcr_ids(ALL, only=['par2_human', 'adrb2_human']) == [
        'par2_human', 'adrb2_human']


def test_only_beats_limit():
    assert select_gpcr_ids(ALL, limit=1, only=['par2_human']) == ['par2_human']


def test_unknown_id_raises_with_the_bad_name_in_the_message():
    with pytest.raises(ValueError, match='nosuch_human'):
        select_gpcr_ids(ALL, only=['nosuch_human'])
