from __future__ import annotations

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
