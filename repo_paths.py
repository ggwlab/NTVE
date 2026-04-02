from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parent
FIGURE2_LOADED_DIR = ROOT / "Figure2" / "figure2_loaded_data"


def find_repo_root(start: Path | None = None) -> Path:
    if start is None:
        return ROOT

    here = start.resolve()
    for candidate in (here.parent, here.parent.parent, here.parent.parent.parent):
        if (candidate / "resources").exists() and (candidate / "repo_paths.py").exists():
            return candidate

    raise FileNotFoundError(f"Could not locate repo root from {start}")
