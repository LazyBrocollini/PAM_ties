from io import BytesIO

import pandas as pd

from src.io_utils import parse_workbook


def workbook_bytes():
    output = BytesIO()
    all_guides = pd.DataFrame([
        {"Source Sheet": "approved", "gRNA Name": "g1", "Ionis#": 2000471.0, "Spacer Sequence": "ACGT", "135 Indel%": 30, "45 Indel%": 20, "11 Indel%": "", "0 Indel%": 5, "PAM": "NGG", "Notes": "first"},
        {"Source Sheet": "approved", "gRNA Name": "g1-alt", "Ionis#": "2000471", "Spacer Sequence": "ACGT", "135 Indel%": 35, "45 Indel%": "", "11 Indel%": 4, "0 Indel%": "", "PAM": "NGG", "Notes": "duplicate"},
        {"Source Sheet": "approved", "gRNA Name": "fallback", "Ionis#": "other", "Spacer Sequence": "TTTT", "135 Indel%": 8, "45 Indel%": 2, "11 Indel%": 1, "0 Indel%": 0, "PAM": "NGG", "Notes": "fallback"},
    ])
    spacer = pd.DataFrame([
        {"Ion_ID": "2000471", "spacer": " ACGT ", "nuclease": "SpCas9", "strand": "+", "PAM": "NGG", "pam_site": "legacy-1", "cut_site": 3},
        {"Ion_ID": "not-in-all", "spacer": "TTTT", "nuclease": "SpCas9", "strand": "+", "PAM": "NGG", "pam_site": "legacy-2", "cut_site": 3},
    ])
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        all_guides.to_excel(writer, sheet_name="All Guides", index=False)
        spacer.to_excel(writer, sheet_name="all_spacer_data.registered_JAG1", index=False)
    return output.getvalue()


def test_ionis_join_duplicate_aggregation_and_spacer_fallback():
    info = parse_workbook(workbook_bytes())
    guides = info["guides"]
    assert len(guides) == 2
    joined = guides[guides["spacer"] == "ACGT"].iloc[0]
    assert "2000471" in joined["Ionis#"]
    assert "2" in joined["Source Row"] and "3" in joined["Source Row"]
    assert "g1" in joined["gRNA Name"] and "g1-alt" in joined["gRNA Name"]
    fallback = guides[guides["spacer"] == "TTTT"].iloc[0]
    assert fallback["join_method"] == "Spacer Sequence fallback"
    assert "all_spacer_data.registered_hJAG" in " ".join(info["warnings"])
