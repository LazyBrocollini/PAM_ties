import pandas as pd

from src.nearby_pam import find_nearby
from src.ranking import rank_hits, numeric_indel


def make_hits(values):
    return pd.DataFrame([
        {"hit_id": f"h{i}", "gRNA Name": f"g{i}", "Ionis#": str(i), "source_sheet": "S", "source_nuclease": "Cas9", "strand": "+", "pam_anchor": 10+i, "135 Indel%": value, "spacer": "ACGT", "pam": "AGG"}
        for i, value in enumerate(values)
    ])


def test_highest_indel_and_ties_and_missing_values():
    hits = make_hits(["10", "20", "bad"])
    pairs, groups, links = find_nearby(hits)
    ranked, best = rank_hits(hits, pairs, groups, links, "135 Indel%")
    assert ranked.loc[1, "is_pair_or_group_winner"] == "Yes"
    assert numeric_indel("bad") is None
    assert any(best["winner_hit_id"] == "h1")

    tie_hits = make_hits(["20", "20"])
    pairs, groups, links = find_nearby(tie_hits)
    ranked, best = rank_hits(tie_hits, pairs, groups, links, "135 Indel%")
    assert (ranked["is_pair_or_group_winner"] == "Yes").sum() == 2
    assert "Tie" in " ".join(best["winner_reason"].astype(str))


def test_all_missing_indel_is_explicit():
    hits = make_hits(["", "not numeric"])
    pairs, groups, links = find_nearby(hits)
    ranked, best = rank_hits(hits, pairs, groups, links, "135 Indel%")
    assert all(ranked["is_pair_or_group_winner"] == "No")
    assert "No numeric Indel% available" in set(best["winner_reason"])
