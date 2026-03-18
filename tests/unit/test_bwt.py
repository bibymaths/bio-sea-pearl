from bio_sea_pearl.api import build_fm_index, search_fm_index


def test_bwt_search():
    seq = "banana"
    idx = build_fm_index(seq)
    # search for pattern "ana"
    positions = search_fm_index(idx, "ana")
    # In "banana", "ana" occurs starting at positions 1 and 3
    assert sorted(positions) == [1, 3]