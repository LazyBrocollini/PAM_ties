from __future__ import annotations

# Generated single-file deployment; rebuild with scripts/build_standalone.py


# --- src/sequence_utils.py ---

import re
from typing import Iterable

DNA_BASES = set("ACGTN")
_IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "N": set("ACGTN"), "V": set("ACG"), "B": set("CGT"), "D": set("AGT"),
    "H": set("ACT"), "K": set("GT"), "M": set("AC"), "R": set("AG"),
    "S": set("CG"), "W": set("AT"), "Y": set("CT"),
}


def normalize_sequence(value: object) -> str:
    """Normalize a DNA-like sequence, removing FASTA decoration and whitespace."""
    if value is None:
        return ""
    text = str(value).replace("\ufeff", "").strip()
    if not text:
        return ""
    lines = text.splitlines()
    if any(line.lstrip().startswith(">") for line in lines):
        lines = [line for line in lines if not line.lstrip().startswith(">")]
    text = "".join(lines)
    text = re.sub(r"\s+", "", text).upper()
    # Common formatting/gap characters are not sequence content.
    return text.replace("-", "").replace(".", "").replace("|", "")


def normalize_reference(value: object) -> tuple[str, str | None]:
    sequence = normalize_sequence(value)
    if not sequence:
        return "", "Reference sequence is empty after normalization."
    invalid = sorted(set(sequence) - DNA_BASES)
    if invalid:
        return "", f"Reference sequence contains invalid base(s): {', '.join(invalid)}. Allowed: A, C, G, T, N."
    return sequence, None


def reverse_complement(sequence: str) -> str:
    sequence = normalize_sequence(sequence)
    complement = str.maketrans("ACGTN", "TGCAN")
    return sequence.translate(complement)[::-1]


def find_all_exact(sequence: str, pattern: str) -> list[int]:
    """Return 0-based starts; N in either sequence acts as an ambiguity wildcard."""
    sequence = normalize_sequence(sequence)
    pattern = normalize_sequence(pattern)
    if not pattern or len(pattern) > len(sequence):
        return []
    hits: list[int] = []
    for start in range(len(sequence) - len(pattern) + 1):
        window = sequence[start : start + len(pattern)]
        if all(a == b or a == "N" or b == "N" for a, b in zip(window, pattern)):
            hits.append(start)
    return hits


def iupac_matches(observed: str, pattern: str) -> bool:
    observed = normalize_sequence(observed)
    pattern = normalize_sequence(pattern)
    if len(observed) != len(pattern) or not observed:
        return False
    return all(base in _IUPAC.get(code, {code}) for base, code in zip(observed, pattern))


def clean_identifier(value: object) -> str:
    """Normalize Excel IDs, including numeric cells rendered as ``2000471.0``."""
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    text = str(value).strip()
    if re.fullmatch(r"[+-]?\d+\.0+", text):
        return text.split(".", 1)[0]
    return text


def display_values(values: Iterable[object]) -> str:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = clean_identifier(value)
        if text and text not in seen:
            result.append(text)
            seen.add(text)
    return "; ".join(result)

# --- src/io_utils.py ---

from collections import defaultdict
from io import BytesIO
from pathlib import Path
from typing import Any

import pandas as pd


ALL_GUIDES_SHEET = "All Guides"
SPACER_SHEETS = ("all_spacer_data.registered_JAG1", "all_spacer_data.registered_hJAG")
INDel_FIELDS = ("135 Indel%", "45 Indel%", "11 Indel%", "0 Indel%")


def _norm_col(value: object) -> str:
    return "".join(ch.lower() for ch in str(value).strip() if ch.isalnum())


def _find_column(columns: list[object], candidates: list[str]) -> str | None:
    lookup = {_norm_col(col): str(col) for col in columns}
    for candidate in candidates:
        if _norm_col(candidate) in lookup:
            return lookup[_norm_col(candidate)]
    return None


def _value(row: pd.Series, column: str | None) -> object:
    return row[column] if column and column in row.index else ""


def _nonempty(values: list[object]) -> list[str]:
    return [clean_identifier(v) for v in values if clean_identifier(v)]


def _aggregate(rows: list[pd.Series], column: str | None, normalize_seq: bool = False) -> str:
    vals = []
    for row in rows:
        raw = _value(row, column)
        val = normalize_sequence(raw) if normalize_seq else clean_identifier(raw)
        if val:
            vals.append(val)
    return display_values(vals)


def _read_excel(source: bytes | str | Path) -> tuple[pd.ExcelFile, str]:
    if isinstance(source, (str, Path)):
        return pd.ExcelFile(source), Path(source).name
    return pd.ExcelFile(BytesIO(source)), "uploaded_workbook.xlsx"


def parse_workbook(source: bytes | str | Path) -> dict[str, Any]:
    """Parse, normalize, and join workbook data while preserving provenance."""
    xls, workbook_name = _read_excel(source)
    sheet_names = list(xls.sheet_names)
    warnings: list[str] = []
    checks: list[dict[str, Any]] = []
    if ALL_GUIDES_SHEET not in sheet_names:
        warnings.append(f"Missing required sheet: {ALL_GUIDES_SHEET}")
        return {"workbook_name": workbook_name, "sheet_names": sheet_names, "guides": pd.DataFrame(),
                "warnings": warnings, "checks": checks, "unmatched_all_guides": pd.DataFrame(),
                "missing_spacers": pd.DataFrame(), "duplicate_identifiers": pd.DataFrame(),
                "available_nucleases": [], "available_strands": [], "available_sources": []}

    all_df = pd.read_excel(xls, sheet_name=ALL_GUIDES_SHEET, dtype=object)
    all_columns = list(all_df.columns)
    ion_col = _find_column(all_columns, ["Ionis#", "Ion_ID", "Ion ID", "Ionis ID"])
    spacer_col = _find_column(all_columns, ["Spacer Sequence", "spacer", "Spacer"])
    if not ion_col and not spacer_col:
        warnings.append("All Guides is missing both an Ionis#/Ion_ID column and a Spacer Sequence column.")
        return {"workbook_name": workbook_name, "sheet_names": sheet_names, "guides": pd.DataFrame(),
                "warnings": warnings, "checks": checks, "unmatched_all_guides": all_df,
                "missing_spacers": pd.DataFrame(), "duplicate_identifiers": pd.DataFrame(),
                "available_nucleases": [], "available_strands": [], "available_sources": []}

    all_rows: list[dict[str, Any]] = []
    for idx, row in all_df.iterrows():
        record = {str(col): row[col] for col in all_columns}
        record.update({
            "_all_row": int(idx) + 2,
            "_ion_key": clean_identifier(_value(row, ion_col)),
            "_spacer_key": normalize_sequence(_value(row, spacer_col)),
        })
        all_rows.append(record)

    spacer_frames: list[pd.DataFrame] = []
    for sheet in SPACER_SHEETS:
        if sheet not in sheet_names:
            warnings.append(f"Optional spacer-data sheet not found: {sheet}")
            continue
        frame = pd.read_excel(xls, sheet_name=sheet, dtype=object)
        frame["_source_sheet"] = sheet
        frame["_spacer_source_row"] = frame.index + 2
        spacer_frames.append(frame)
    if not spacer_frames:
        warnings.append("Neither expected spacer-data sheet was available; no joined guides can be analyzed.")
        return {"workbook_name": workbook_name, "sheet_names": sheet_names, "guides": pd.DataFrame(),
                "warnings": warnings, "checks": checks, "unmatched_all_guides": all_df,
                "missing_spacers": pd.DataFrame(), "duplicate_identifiers": pd.DataFrame(),
                "available_nucleases": [], "available_strands": [], "available_sources": []}

    spacer_df = pd.concat(spacer_frames, ignore_index=True, sort=False)
    sp_ion_col = _find_column(list(spacer_df.columns), ["Ion_ID", "Ionis#", "Ion ID", "Ionis ID"])
    sp_seq_col = _find_column(list(spacer_df.columns), ["spacer", "Spacer Sequence", "Spacer"])
    nuclease_col = _find_column(list(spacer_df.columns), ["nuclease", "enzyme", "Nuclease/Enzyme", "Nuclease"])
    strand_col = _find_column(list(spacer_df.columns), ["strand", "guide strand", "orientation"])
    pam_col = _find_column(list(spacer_df.columns), ["PAM", "pam"])
    pam_site_col = _find_column(list(spacer_df.columns), ["pam_site", "PAM site", "pam site"])
    cut_site_col = _find_column(list(spacer_df.columns), ["cut_site", "cut site"])
    if not sp_seq_col:
        warnings.append("Spacer-data sheets are missing a spacer column.")

    # Deduplicate spacer rows by the requested normalized tuple, preserving row provenance.
    dedup: dict[tuple[str, str, str, str, str], list[pd.Series]] = defaultdict(list)
    for _, row in spacer_df.iterrows():
        key = (
            normalize_sequence(_value(row, sp_seq_col)),
            clean_identifier(_value(row, nuclease_col)),
            clean_identifier(_value(row, strand_col)),
            clean_identifier(_value(row, pam_site_col)),
            clean_identifier(_value(row, cut_site_col)),
        )
        dedup[key].append(row)

    joined: list[dict[str, Any]] = []
    matched_all_rows: set[int] = set()
    missing_spacers: list[dict[str, Any]] = []
    duplicate_identifier_rows: list[dict[str, Any]] = []
    all_lookup_ion: dict[str, list[dict[str, Any]]] = defaultdict(list)
    all_lookup_spacer: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for rec in all_rows:
        if rec["_ion_key"]:
            all_lookup_ion[rec["_ion_key"]].append(rec)
        if rec["_spacer_key"]:
            all_lookup_spacer[rec["_spacer_key"]].append(rec)

    for key, rows in dedup.items():
        row = rows[0]
        spacer = key[0]
        if not spacer:
            missing_spacers.append({"source_sheet": row["_source_sheet"], "source_row": int(row["_spacer_source_row"]), "reason": "Missing spacer sequence"})
            continue
        ion_key = clean_identifier(_value(row, sp_ion_col))
        candidates = all_lookup_ion.get(ion_key, []) if ion_key else []
        join_method = "Ionis#"
        if not candidates:
            candidates = all_lookup_spacer.get(spacer, [])
            join_method = "Spacer Sequence fallback"
        if not candidates:
            continue
        for candidate in candidates:
            matched_all_rows.add(candidate["_all_row"])
        all_source_rows = sorted({c["_all_row"] for c in candidates})
        all_ionis = _aggregate([pd.Series(c) for c in candidates], ion_col)
        all_spacers = _aggregate([pd.Series(c) for c in candidates], spacer_col, normalize_seq=True)
        source_fields = {
            "Source Sheet": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["Source Sheet", "source sheet"])),
            "Source Row": "; ".join(str(x) for x in all_source_rows),
            "source_row": "; ".join(str(x) for x in all_source_rows),
            "gRNA Name": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["gRNA Name", "gRNA_Name", "guide name"])),
            "Ionis#": all_ionis or ion_key,
            "Description": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["Description"])),
            "135 Indel%": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["135 Indel%", "135 Indel"])),
            "45 Indel%": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["45 Indel%", "45 Indel"])),
            "11 Indel%": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["11 Indel%", "11 Indel"])),
            "0 Indel%": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["0 Indel%", "0 Indel"])),
            "Spacer Sequence": all_spacers or spacer,
            "PAM": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["PAM", "pam"])),
            "Notes": _aggregate([pd.Series(c) for c in candidates], _find_column(all_columns, ["Notes", "Note"])),
        }
        base = {str(k): v for k, v in row.items() if not str(k).startswith("_")}
        base.update({
            "source_sheet": str(row["_source_sheet"]),
            "spacer_source_rows": "; ".join(str(int(r["_spacer_source_row"])) for r in rows),
            "spacer": spacer,
            "source_ion_id": ion_key,
            "source_spacer": spacer,
            "source_nuclease": clean_identifier(_value(row, nuclease_col)),
            "nuclease": clean_identifier(_value(row, nuclease_col)),
            "source_strand": clean_identifier(_value(row, strand_col)),
            "source_pam": clean_identifier(_value(row, pam_col)),
            "source_pam_site": clean_identifier(_value(row, pam_site_col)),
            "source_cut_site": clean_identifier(_value(row, cut_site_col)),
            "join_method": join_method,
        })
        # Keep all spacer-data values available without collisions with the normalized result fields.
        for col in spacer_df.columns:
            if col not in base and not str(col).startswith("_"):
                base[f"spacer_data__{col}"] = display_values([r[col] for r in rows])
        base.update(source_fields)
        joined.append(base)

    unmatched = all_df.iloc[[i - 2 for i in range(2, len(all_df) + 2) if i not in matched_all_rows]].copy()
    if not unmatched.empty:
        warnings.append(f"{len(unmatched)} All Guides row(s) did not join to spacer data.")
    for column, label in ((nuclease_col, "nuclease/enzyme"), (strand_col, "strand"), (pam_col, "PAM")):
        if not column:
            warnings.append(f"Spacer-data sheets are missing an optional annotation column: {label}.")
    duplicate_ids = all_df[all_df.apply(lambda r: clean_identifier(_value(r, ion_col)), axis=1).duplicated(keep=False)] if ion_col else pd.DataFrame()
    if not duplicate_ids.empty:
        duplicate_identifier_rows = duplicate_ids.to_dict("records")
        warnings.append(f"{len(duplicate_ids)} All Guides row(s) share duplicate Ionis# identifiers; values were aggregated.")

    guides = pd.DataFrame(joined)
    if guides.empty:
        warnings.append("No whitelisted guides could be joined from All Guides to the spacer-data sheets.")
    for field in INDel_FIELDS:
        if field not in guides.columns:
            guides[field] = ""
        def _numeric_present(value: object) -> bool:
            return any(token.strip().replace("%", "") for token in str(value or "").split(";") if token.strip().replace("%", "").strip()) and any(
                _is_float(token.strip().replace("%", "")) for token in str(value or "").split(";")
            )
        def _is_float(token: str) -> bool:
            try:
                float(token)
                return True
            except (TypeError, ValueError):
                return False
        missing_count = int(sum(not _numeric_present(value) for value in guides[field])) if not guides.empty else 0
        if missing_count:
            warnings.append(f"{missing_count} joined guide row(s) have no numeric value for {field}.")
    for col in ["source_nuclease", "source_strand", "source_pam"]:
        if col not in guides.columns:
            guides[col] = ""
    checks.extend([
        {"check": "Workbook sheets", "status": "warning" if warnings else "pass", "detail": "; ".join(warnings) if warnings else "Required and optional sheets loaded."},
        {"check": "Missing spacer sequences", "status": "warning" if missing_spacers else "pass", "detail": f"{len(missing_spacers)} spacer-data row(s) missing a spacer."},
        {"check": "Unmatched All Guides rows", "status": "warning" if not unmatched.empty else "pass", "detail": f"{len(unmatched)} row(s) did not join."},
        {"check": "Duplicate identifiers", "status": "warning" if duplicate_identifier_rows else "pass", "detail": f"{len(duplicate_identifier_rows)} duplicate-ID row(s); distinct values preserved."},
    ])
    available = lambda col: sorted([x for x in guides.get(col, pd.Series(dtype=object)).dropna().astype(str).unique().tolist() if x and x != "nan"])
    return {"workbook_name": workbook_name, "sheet_names": sheet_names, "guides": guides,
            "warnings": warnings, "checks": checks, "unmatched_all_guides": unmatched,
            "missing_spacers": pd.DataFrame(missing_spacers), "duplicate_identifiers": pd.DataFrame(duplicate_identifier_rows),
            "available_nucleases": available("source_nuclease"), "available_strands": available("source_strand"),
            "available_sources": sorted(guides.get("source_sheet", pd.Series(dtype=object)).dropna().astype(str).unique().tolist())}

# --- src/guide_matching.py ---

import re
from typing import Any

import pandas as pd


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

# --- src/nearby_pam.py ---

from collections import defaultdict
from typing import Any

import pandas as pd


def _norm(value: object) -> str:
    return str(value or "").strip().lower()


def find_nearby(hits: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[str]]]:
    """Return pair rows, connected-component group rows, and hit-to-entity maps."""
    if hits.empty:
        return pd.DataFrame(), pd.DataFrame(), {"hit_to_pairs": {}, "hit_to_groups": {}}
    pairs: list[dict[str, Any]] = []
    parent = {str(idx): str(idx) for idx in hits.index}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(len(hits)):
        a = hits.iloc[i]
        for j in range(i + 1, len(hits)):
            b = hits.iloc[j]
            if a.get("pam_anchor") is None or b.get("pam_anchor") is None:
                continue
            if pd.isna(a.get("pam_anchor")) or pd.isna(b.get("pam_anchor")):
                continue
            if (_norm(a.get("source_sheet")) != _norm(b.get("source_sheet")) or
                _norm(a.get("source_nuclease")) != _norm(b.get("source_nuclease")) or
                _norm(a.get("strand")) != _norm(b.get("strand"))):
                continue
            distance = abs(int(a["pam_anchor"]) - int(b["pam_anchor"]))
            if distance != 1:
                continue
            union(str(a.name), str(b.name))
            pair_id = f"P{len(pairs) + 1:04d}"
            pairs.append({
                "pair_id": pair_id,
                "hit_id_a": a["hit_id"], "hit_id_b": b["hit_id"],
                "source_sheet": a.get("source_sheet", ""), "nuclease": a.get("source_nuclease", ""), "strand": a.get("strand", ""),
                "pam_anchor_a": int(a["pam_anchor"]), "pam_anchor_b": int(b["pam_anchor"]),
                "pam_distance_nt": distance,
                "index_a": int(a.name), "index_b": int(b.name),
            })

    root_to_indices: dict[str, list[int]] = defaultdict(list)
    for idx in hits.index:
        root_to_indices[find(str(idx))].append(int(idx))
    components = sorted(root_to_indices.values(), key=lambda x: min(x))
    group_for_index: dict[int, str] = {}
    groups: list[dict[str, Any]] = []
    for n, indices in enumerate(components, 1):
        if not any(any(p["index_a"] == idx or p["index_b"] == idx for p in pairs) for idx in indices):
            continue
        group_id = f"G{n:04d}"
        for idx in indices:
            group_for_index[idx] = group_id
        group_hits = hits.loc[indices]
        groups.append({
            "group_id": group_id, "hit_indices": indices, "hit_ids": "; ".join(group_hits["hit_id"].astype(str)),
            "source_sheet": group_hits.iloc[0].get("source_sheet", ""), "nuclease": group_hits.iloc[0].get("source_nuclease", ""),
            "strand": group_hits.iloc[0].get("strand", ""), "pam_anchor_min": int(group_hits["pam_anchor"].min()),
            "pam_anchor_max": int(group_hits["pam_anchor"].max()), "group_size": len(indices),
        })
    hit_to_pairs = defaultdict(list)
    for p in pairs:
        hit_to_pairs[p["hit_id_a"]].append(p["pair_id"])
        hit_to_pairs[p["hit_id_b"]].append(p["pair_id"])
    hit_to_groups = {hits.loc[idx, "hit_id"]: group_id for idx, group_id in group_for_index.items()}
    pair_df = pd.DataFrame(pairs)
    group_df = pd.DataFrame(groups)
    return pair_df, group_df, {"hit_to_pairs": dict(hit_to_pairs), "hit_to_groups": hit_to_groups}

# --- src/ranking.py ---

from typing import Any

import pandas as pd


def numeric_indel(value: object) -> float | None:
    """Return the highest numeric value in a preserved semicolon-delimited field."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    values: list[float] = []
    for token in str(value).replace("%", "").split(";"):
        token = token.strip()
        if not token:
            continue
        try:
            parsed = float(token)
        except (TypeError, ValueError):
            continue
        if pd.notna(parsed):
            values.append(parsed)
    return max(values) if values else None


def rank_hits(hits: pd.DataFrame, pairs: pd.DataFrame, groups: pd.DataFrame, links: dict[str, Any], metric: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if hits.empty:
        return hits.copy(), pd.DataFrame()
    out = hits.copy()
    out["selected_indel_metric"] = metric
    out["selected_indel_value"] = out.get(metric, pd.Series(index=out.index, dtype=object)).map(numeric_indel)
    out["group_id"] = out["hit_id"].map(links.get("hit_to_groups", {})).fillna("")
    pair_map = links.get("hit_to_pairs", {})
    out["pair_id"] = out["hit_id"].map(lambda x: "; ".join(pair_map.get(x, [])))
    out["pam_distance_nt"] = out["hit_id"].map(lambda x: 1 if pair_map.get(x, []) else "")
    out["is_pair_or_group_winner"] = "No"
    out["winner_reason"] = ""
    best_rows: list[dict[str, Any]] = []

    def choose(indices: list[int], entity_type: str, entity_id: str) -> None:
        vals = [(idx, out.loc[idx, "selected_indel_value"]) for idx in indices]
        numeric = [(idx, val) for idx, val in vals if val is not None]
        if not numeric:
            winner_indices: list[int] = []
            reason = "No numeric Indel% available"
            winner_text = "No numeric Indel% available"
            value_text = ""
        else:
            max_value = max(val for _, val in numeric)
            winner_indices = [idx for idx, val in numeric if val == max_value]
            reason = "Tie for highest selected Indel%" if len(winner_indices) > 1 else "Highest selected Indel%"
            winner_text = "; ".join(str(out.loc[idx, "hit_id"]) for idx in winner_indices)
            value_text = max_value
        for idx in winner_indices:
            out.at[idx, "is_pair_or_group_winner"] = "Yes"
            prior = out.at[idx, "winner_reason"]
            out.at[idx, "winner_reason"] = "; ".join(x for x in [prior, f"{entity_type} winner" + (" (tie)" if len(winner_indices) > 1 else "")] if x)
        for idx in winner_indices or [None]:
            row = {"entity_type": entity_type, "entity_id": entity_id, "winner_hit_id": out.loc[idx, "hit_id"] if idx is not None else "",
                   "winner_reason": reason, "selected_indel_metric": metric, "selected_indel_value": value_text}
            if idx is not None:
                for field in ("gRNA Name", "Ionis#", "source_nuclease", "strand", "spacer", "pam", "spacer_start", "spacer_end", "pam_start", "pam_end", "pam_anchor", "source_row", "Source Row"):
                    if field in out.columns:
                        row[field] = out.loc[idx, field]
            best_rows.append(row)

    for _, pair in pairs.iterrows():
        choose([int(pair["index_a"]), int(pair["index_b"])], "pair", str(pair["pair_id"]))
    for _, group in groups.iterrows():
        choose([int(x) for x in group["hit_indices"]], "group", str(group["group_id"]))
    return out, pd.DataFrame(best_rows)

# --- src/visualization.py ---

from html import escape
from typing import Any

import pandas as pd


def sequence_html(reference: str, hits: pd.DataFrame, start: int = 1, end: int | None = None) -> str:
    end = end or len(reference)
    start = max(1, start)
    end = min(len(reference), end)
    colors = {"+": "#13795b", "-": "#7c3aed"}
    spans: list[str] = []
    for pos in range(start, end + 1):
        base = escape(reference[pos - 1])
        labels: list[str] = []
        css = ""
        if not hits.empty:
            for _, hit in hits.iterrows():
                in_spacer = hit.get("spacer_start") is not None and int(hit["spacer_start"]) <= pos <= int(hit["spacer_end"])
                pam_start, pam_end = hit.get("pam_start"), hit.get("pam_end")
                in_pam = pam_start is not None and pd.notna(pam_start) and pam_end is not None and pd.notna(pam_end) and int(pam_start) <= pos <= int(pam_end)
                if in_spacer or in_pam:
                    strand = str(hit.get("strand", "+"))
                    css = f"background:{colors.get(strand, '#555')};color:white;"
                    labels.append(f"{hit.get('gRNA Name', hit.get('Ionis#', 'guide'))} {strand} {'PAM' if in_pam else 'spacer'}")
        title = escape(" | ".join(labels)) if labels else ""
        spans.append(f'<span style="{css}" title="{title}">{base}</span>')
    return '<div style="font-family:monospace;line-height:1.8;word-break:break-all;border:1px solid #d1d5db;padding:12px;border-radius:8px;">' + "".join(spans) + "</div>"


def pair_table(hits: pd.DataFrame, pairs: pd.DataFrame, metric: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, pair in pairs.iterrows():
        a = hits.loc[int(pair["index_a"])]
        b = hits.loc[int(pair["index_b"])]
        av, bv = a.get("selected_indel_value"), b.get("selected_indel_value")
        if av is None and bv is None:
            winner = "No numeric Indel% available"
        elif av is None:
            winner = str(b.get("gRNA Name", b.get("Ionis#", b["hit_id"])))
        elif bv is None or av > bv:
            winner = str(a.get("gRNA Name", a.get("Ionis#", a["hit_id"])))
        elif bv > av:
            winner = str(b.get("gRNA Name", b.get("Ionis#", b["hit_id"])))
        else:
            winner = f"Tie: {a.get('gRNA Name', a['hit_id'])}; {b.get('gRNA Name', b['hit_id'])}"
        rows.append({"pair_id": pair["pair_id"], "guide_a": a.get("gRNA Name", a.get("Ionis#", "")), "guide_b": b.get("gRNA Name", b.get("Ionis#", "")),
                     "Ionis# A": a.get("Ionis#", ""), "Ionis# B": b.get("Ionis#", ""), "spacer A": a.get("spacer", ""), "spacer B": b.get("spacer", ""),
                     "strand": pair.get("strand", ""), "nuclease": pair.get("nuclease", ""), "PAM anchor A": pair["pam_anchor_a"], "PAM anchor B": pair["pam_anchor_b"],
                     "exact PAM distance (nt)": pair["pam_distance_nt"], f"{metric} A": av, f"{metric} B": bv, "Winner": winner})
    return pd.DataFrame(rows)

# --- Streamlit interface ---

from io import BytesIO
from pathlib import Path
import sys

import pandas as pd
import streamlit as st

# Resolve sibling modules independently of the deployment working directory.
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


st.set_page_config(page_title="Nearby PAM Guide Selector", page_icon="🧬", layout="wide")
DEFAULT_NAME = "Sum_samples_guides_hJAG1_dKO_consolidated.xlsx"


def default_path() -> Path | None:
    candidates = [Path(__file__).parent / "fixtures" / DEFAULT_NAME, Path.cwd() / DEFAULT_NAME,
                  Path("/home/user/output") / DEFAULT_NAME, Path("/mnt/data") / DEFAULT_NAME]
    return next((p for p in candidates if p.exists()), None)


@st.cache_data(show_spinner=False)
def cached_parse_workbook(content: bytes, filename: str):
    return parse_workbook(content)


@st.cache_data(show_spinner=False)
def cached_match(reference: str, guides: pd.DataFrame):
    return match_guides(reference, guides)


def excel_bytes(tables: dict[str, pd.DataFrame]) -> bytes:
    output = BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        for name, table in tables.items():
            table.to_excel(writer, sheet_name=name[:31], index=False)
    return output.getvalue()


def warning_list(workbook_info, match_warnings, reference_error):
    warnings = list(workbook_info.get("warnings", [])) + list(match_warnings)
    if reference_error:
        warnings.append(reference_error)
    return sorted(set(warnings))


def main():
    st.title("Nearby PAM Guide Selector")
    st.caption("Find whitelisted guides in a supplied reference sequence, compare exact 1-nt PAM-anchor neighbors, and rank each pair or connected group by the selected Indel% metric.")
    with st.sidebar:
        st.header("Inputs")
        uploaded_workbook = st.file_uploader("Guide workbook (.xlsx)", type=["xlsx"], help=f"Expected default filename: {DEFAULT_NAME}")
        fixture = default_path()
        if uploaded_workbook is not None:
            workbook_content, workbook_name = uploaded_workbook.getvalue(), uploaded_workbook.name
        elif fixture:
            workbook_content, workbook_name = fixture.read_bytes(), fixture.name
            st.info(f"Using default fixture: {fixture}")
        else:
            workbook_content, workbook_name = None, ""
            st.warning(f"Upload a workbook. The app will also auto-load {DEFAULT_NAME} if placed in the project fixtures folder.")
        sequence_text = st.text_area("Raw gDNA sequence", height=140, placeholder="Paste A/C/G/T/N sequence or FASTA text")
        sequence_file = st.file_uploader("Reference FASTA or plain-text file", type=["fa", "fasta", "txt"], help="If both inputs are supplied, this uploaded file takes precedence.")
        metric = st.selectbox("Indel% metric", list(INDel_FIELDS), index=0)
        workbook_info = cached_parse_workbook(workbook_content, workbook_name) if workbook_content else {"guides": pd.DataFrame(), "warnings": [], "available_nucleases": [], "available_strands": [], "available_sources": [], "checks": []}
        nuclease_filter = st.multiselect("Nuclease/enzyme", workbook_info.get("available_nucleases", []))
        strand_filter = st.multiselect("Source strand annotation", workbook_info.get("available_strands", []))
        source_filter = st.multiselect("Source dataset/sheet", workbook_info.get("available_sources", []))
        analyze = st.button("Analyze nearby PAM guides", type="primary", use_container_width=True)

    if not analyze:
        st.info("Provide a workbook and reference sequence, then select Analyze nearby PAM guides.")
        if workbook_info.get("guides") is not None and not workbook_info["guides"].empty:
            st.metric("Whitelisted guides loaded", len(workbook_info["guides"]))
        return
    if not workbook_content:
        st.error("No workbook was provided.")
        return
    if sequence_file is not None:
        raw_reference = sequence_file.getvalue().decode("utf-8", errors="replace")
        st.info("The uploaded FASTA/plain-text file takes precedence over the text area.")
    else:
        raw_reference = sequence_text
    reference, reference_error = normalize_reference(raw_reference)
    if reference_error:
        st.error(reference_error)
        return

    guides = workbook_info["guides"].copy()
    if nuclease_filter:
        guides = guides[guides["source_nuclease"].isin(nuclease_filter)]
    if strand_filter:
        guides = guides[guides["source_strand"].isin(strand_filter)]
    if source_filter:
        guides = guides[guides["source_sheet"].isin(source_filter)]
    hits, match_warnings = cached_match(reference, guides)
    pairs, groups, links = find_nearby(hits)
    hits, best = rank_hits(hits, pairs, groups, links, metric)
    warnings = warning_list(workbook_info, match_warnings, None)

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Reference length", f"{len(reference):,} nt")
    c2.metric("Whitelisted guides loaded", len(guides))
    c3.metric("Guide hits", len(hits))
    c4.metric("Exact 1-nt pairs", len(pairs))
    c5.metric("Connected groups", len(groups))
    st.caption(f"Workbook: {workbook_name} · sheets used: {', '.join(workbook_info.get('sheet_names', []))} · PAM comparison uses exact abs(anchor A − anchor B) = 1 nt; 2-nt differences are excluded.")
    if warnings:
        with st.expander(f"Warnings ({len(warnings)})", expanded=False):
            for warning in warnings:
                st.warning(warning)

    tab_best, tab_pairs, tab_seq, tab_all, tab_checks = st.tabs(["Best guides", "Nearby pairs", "Sequence view", "All matched guides", "Data checks"])
    with tab_best:
        st.subheader("Highest selected Indel% within each pair or connected group")
        st.dataframe(best, use_container_width=True, hide_index=True)
    with tab_pairs:
        st.subheader("Qualifying nearby pairs")
        st.dataframe(pair_table(hits, pairs, metric), use_container_width=True, hide_index=True)
    with tab_seq:
        st.subheader("Reference sequence view")
        st.markdown('<span style="background:#13795b;color:white;padding:3px 8px;border-radius:4px;">Forward (+)</span> <span style="background:#7c3aed;color:white;padding:3px 8px;border-radius:4px;">Reverse (−)</span> &nbsp; highlighted regions include matched spacers and validated/attempted PAMs.', unsafe_allow_html=True)
        if hits.empty:
            st.info("No guide hits were found in the submitted reference.")
        else:
            window = st.slider("Coordinate window", 1, len(reference), (1, min(len(reference), 500)))
            st.markdown(sequence_html(reference, hits, window[0], window[1]), unsafe_allow_html=True)
            st.dataframe(hits[[c for c in ["hit_id", "gRNA Name", "Ionis#", "strand", "spacer_start", "spacer_end", "pam_start", "pam_end", "pam_anchor", "pam", "pam_valid"] if c in hits.columns]], use_container_width=True, hide_index=True)
    with tab_all:
        st.subheader("Every matched whitelisted guide hit")
        st.dataframe(hits, use_container_width=True, hide_index=True)
    with tab_checks:
        st.subheader("Data checks and provenance")
        check_rows = list(workbook_info.get("checks", []))
        check_rows += [{"check": "PAM validation", "status": "warning" if any(hits.get("pam_valid", pd.Series(dtype=bool)) == False) else "pass", "detail": f"{int((~hits.get('pam_valid', pd.Series(dtype=bool))).sum()) if not hits.empty and 'pam_valid' in hits else 0} hit(s) lacked a validated PAM."}]
        check_rows += [{"check": "Reference sequence", "status": "pass", "detail": f"Normalized and validated as {len(reference):,} nt."},
                       {"check": "Nearby-PAM rule", "status": "pass", "detail": "Only same source sheet, nuclease, strand, and exact 1-nt PAM-anchor differences are included."}]
        st.dataframe(pd.DataFrame(check_rows), use_container_width=True, hide_index=True)
        if not workbook_info.get("unmatched_all_guides", pd.DataFrame()).empty:
            st.write("Unmatched All Guides rows")
            st.dataframe(workbook_info["unmatched_all_guides"], use_container_width=True, hide_index=True)
        if not workbook_info.get("missing_spacers", pd.DataFrame()).empty:
            st.write("Missing spacer sequences")
            st.dataframe(workbook_info["missing_spacers"], use_container_width=True, hide_index=True)

    export_tables = {"Matched Guides": hits, "Nearby Pairs": pair_table(hits, pairs, metric), "Best Guides": best, "Data Checks": pd.DataFrame(check_rows)}
    st.download_button("Download matched-guide CSV", hits.to_csv(index=False).encode("utf-8"), file_name="nearby_pam_matched_guides.csv", mime="text/csv")
    st.download_button("Download Excel results", excel_bytes(export_tables), file_name="nearby_pam_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


if __name__ == "__main__":
    main()
