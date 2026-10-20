from __future__ import annotations

import re
from typing import Any

import pandas as pd

from .sequence_utils import find_all_exact, iupac_matches, reverse_complement

_IUPAC_CODES = set("ACGTNRYSWKMBDHV")


def _pattern_from_annotation(pam: object, nuclease: object) -> tuple[str, str]:
    text = str(pam or "").strip().upper().replace("-", "")
    candidate = "".join(ch for ch in text if ch in _IUPAC_CODES)
    if candidate and len(candidate) >= 2:
        return candidate, "workbook PAM"
    nuc = str(nuclease or "").lower()
    if any(term in nuc for term in ("cas12", "cpf1")):
        return "TTTV", "inferred Cas12/CPF1 PAM"
    if any(term in nuc for term in ("cas9", "spcas9", "s. pyogenes")):
        return "NGG", "inferred Cas9 PAM"
    return "", "PAM annotation unavailable"


def _pam_side(pattern: str, nuclease: object) -> str:
    nuc = str(nuclease or "").lower()
    if any(term in nuc for term in ("cas12", "cpf1")) or pattern.startswith("TTT"):
        return "5prime"
    return "3prime"


def _pam_for_hit(reference: str, spacer_start: int, spacer_end: int, orientation: str, pattern: str, nuclease: object) -> dict[str, Any]:
    length = len(pattern)
    side = _pam_side(pattern, nuclease)
    # Coordinates below are 1-based inclusive. The anchor is the PAM base adjacent to the protospacer.
    if side == "3prime":
        if orientation == "+":
            start, end = spacer_end + 1, spacer_end + length
            anchor = start
        else:
            start, end = spacer_start - length, spacer_start - 1
            anchor = end
    else:  # 5' PAM relative to the guide orientation
        if orientation == "+":
            start, end = spacer_start - length, spacer_start - 1
            anchor = end
        else:
            start, end = spacer_end + 1, spacer_end + length
            anchor = start
    if start < 1 or end > len(reference):
        return {"pam_start": None, "pam_end": None, "pam_anchor": None, "pam_sequence": "",
                "pam_valid": False, "pam_warning": f"PAM coordinates {start}-{end} fall outside the reference."}
    observed_reference = reference[start - 1 : end]
    observed = observed_reference if orientation == "+" else reverse_complement(observed_reference)
    valid = iupac_matches(observed, pattern)
    warning = None if valid else f"Observed PAM {observed} does not match {pattern} at {start}-{end}."
    return {"pam_start": start, "pam_end": end, "pam_anchor": anchor, "pam_sequence": observed,
            "pam_valid": valid, "pam_warning": warning}


def match_guides(reference: str, guides: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Find every forward/reverse hit for each whitelisted, joined guide."""
    hits: list[dict[str, Any]] = []
    warnings: list[str] = []
    if guides.empty:
        return pd.DataFrame(), warnings
    for guide_idx, guide in guides.iterrows():
        spacer = str(guide.get("spacer", "") or "").upper()
        if not spacer:
            warnings.append(f"Missing spacer for source row(s) {guide.get('Source Row', guide_idx)}.")
            continue
        nuclease = str(guide.get("source_nuclease", "") or "")
        pattern, pattern_source = _pattern_from_annotation(guide.get("source_pam") or guide.get("PAM", ""), nuclease)
        forward_starts = find_all_exact(reference, spacer)
        rc = reverse_complement(spacer)
        reverse_starts = find_all_exact(reference, rc)
        seen: set[tuple[int, str]] = set()
        for orientation, starts in (("+", forward_starts), ("-", reverse_starts)):
            for zero_start in starts:
                if (zero_start, orientation) in seen:
                    continue
                seen.add((zero_start, orientation))
                spacer_start = zero_start + 1
                spacer_end = spacer_start + len(spacer) - 1
                pam = _pam_for_hit(reference, spacer_start, spacer_end, orientation, pattern, nuclease) if pattern else {
                    "pam_start": None, "pam_end": None, "pam_anchor": None, "pam_sequence": "", "pam_valid": False,
                    "pam_warning": "No usable PAM pattern was available from the workbook or nuclease annotation."
                }
                if pam["pam_warning"]:
                    warnings.append(f"{guide.get('source_sheet', '')} row {guide.get('spacer_source_rows', '')}: {pam['pam_warning']}")
                hit = {k: v for k, v in guide.items()}
                hit.update({
                    "hit_index": len(hits), "strand": orientation, "spacer_start": spacer_start, "spacer_end": spacer_end,
                    "pam_start": pam["pam_start"], "pam_end": pam["pam_end"], "pam_anchor": pam["pam_anchor"],
                    "pam": pam["pam_sequence"] or guide.get("source_pam", guide.get("PAM", "")),
                    "pam_pattern": pattern, "pam_pattern_source": pattern_source,
                    "pam_valid": pam["pam_valid"], "pam_warning": pam["pam_warning"] or "",
                })
                source_sheet = str(hit.get("source_sheet", ""))
                ion = str(hit.get("Ionis#", hit.get("source_ion_id", "")))
                hit["hit_id"] = "|".join([source_sheet, ion, spacer, orientation, str(spacer_start), str(pam["pam_anchor"] or "NA")])
                hits.append(hit)
    result = pd.DataFrame(hits)
    if len(hits) > len(guides):
        warnings.append(f"Multiple reference hits were found for {len(hits) - len(guides)} additional guide occurrence(s); every exact hit is retained.")
    if not result.empty:
        result = result.drop_duplicates(subset=["hit_id"], keep="first").reset_index(drop=True)
    return result, sorted(set(warnings))
