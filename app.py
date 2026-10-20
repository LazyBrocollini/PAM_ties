from __future__ import annotations

from io import BytesIO
from pathlib import Path
import sys

import pandas as pd
import streamlit as st

# Resolve sibling modules independently of the deployment working directory.
PROJECT_ROOT = Path(__file__).resolve().parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.guide_matching import match_guides
from src.io_utils import INDel_FIELDS, parse_workbook
from src.nearby_pam import find_nearby
from src.ranking import rank_hits
from src.sequence_utils import normalize_reference
from src.visualization import pair_table, sequence_html

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
