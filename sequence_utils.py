from __future__ import annotations

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
