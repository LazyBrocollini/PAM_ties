import pandas as pd

from src.guide_matching import match_guides
from src.nearby_pam import find_nearby


def guide(name, spacer, ion, nuclease="SpCas9", source_sheet="all_spacer_data.registered_JAG1", pam="NGG"):
    return {"gRNA Name": name, "Ionis#": ion, "spacer": spacer, "source_nuclease": nuclease,
            "source_sheet": source_sheet, "source_pam": pam, "Source Row": "2", "135 Indel%": "10"}


def test_forward_reverse_hits_and_coordinates_and_pam_anchors():
    # Forward protospacer at 1-4 with NGG at 5-7; reverse protospacer at 11-14 with NGG upstream at 8-10.
    reference = "ACGTAGGCCCATGCAT"
    guides = pd.DataFrame([guide("g1", "ACGT", "1"), guide("g2", "GCAT", "2")])
    hits, warnings = match_guides(reference, guides)
    fwd = hits[(hits["Ionis#"] == "1") & (hits["strand"] == "+")].iloc[0]
    rev = hits[(hits["Ionis#"] == "2") & (hits["strand"] == "-")].iloc[0]
    assert (fwd["spacer_start"], fwd["spacer_end"]) == (1, 4)
    assert (fwd["pam_start"], fwd["pam_end"], fwd["pam_anchor"]) == (5, 7, 5)
    assert (rev["spacer_start"], rev["spacer_end"]) == (11, 14)
    assert rev["pam_anchor"] == 10


def test_exact_one_nt_only_and_same_nuclease_strand_source_restrictions():
    hits = pd.DataFrame([
        {"hit_id": "a", "source_sheet": "S", "source_nuclease": "Cas9", "strand": "+", "pam_anchor": 10},
        {"hit_id": "b", "source_sheet": "S", "source_nuclease": "Cas9", "strand": "+", "pam_anchor": 11},
        {"hit_id": "c", "source_sheet": "S", "source_nuclease": "Cas9", "strand": "+", "pam_anchor": 13},
        {"hit_id": "d", "source_sheet": "S", "source_nuclease": "Cas12", "strand": "+", "pam_anchor": 11},
        {"hit_id": "e", "source_sheet": "S", "source_nuclease": "Cas9", "strand": "-", "pam_anchor": 11},
        {"hit_id": "f", "source_sheet": "T", "source_nuclease": "Cas9", "strand": "+", "pam_anchor": 12},
    ])
    pairs, groups, links = find_nearby(hits)
    assert {(p["hit_id_a"], p["hit_id_b"]) for _, p in pairs.iterrows()} == {("a", "b")}
    assert len(groups) == 1


def test_duplicate_hit_ids_include_coordinates():
    hits = pd.DataFrame([guide("g1", "ATGC", "2000471")])
    matched, _ = match_guides("ATGCAGCAT", hits)
    assert len(matched) == 2
    assert matched["hit_id"].nunique() == 2
