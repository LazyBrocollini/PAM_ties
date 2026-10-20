from __future__ import annotations

from collections import defaultdict
from io import BytesIO
from pathlib import Path
from typing import Any

import pandas as pd

from .sequence_utils import clean_identifier, display_values, normalize_sequence

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
