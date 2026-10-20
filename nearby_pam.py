from __future__ import annotations

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
