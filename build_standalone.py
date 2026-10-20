"""Generate a single-file deployment without any local package imports."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULES = ("sequence_utils", "io_utils", "guide_matching", "nearby_pam", "ranking", "visualization")

def strip_internal_imports(text):
    return "\n".join(line for line in text.splitlines()
                     if not line.startswith(("from __future__ import", "from .", "from src.")))


def main():
    chunks = ["from __future__ import annotations\n", "# Generated single-file deployment; rebuild with scripts/build_standalone.py\n"]
    for name in MODULES:
        chunks.append(f"\n# --- src/{name}.py ---\n" + strip_internal_imports((ROOT / "src" / f"{name}.py").read_text()))
    chunks.append("\n# --- Streamlit interface ---\n" + strip_internal_imports((ROOT / "app.py").read_text()))
    target = ROOT / "standalone" / "app.py"
    target.parent.mkdir(exist_ok=True)
    target.write_text("\n".join(chunks) + "\n")
    print(target)

if __name__ == "__main__":
    main()
