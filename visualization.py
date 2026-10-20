from __future__ import annotations

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
