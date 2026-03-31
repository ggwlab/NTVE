"""
Figure 7c full Python pipeline.

This script:
1. recreates the spline-model input CSVs from quantfiles
2. runs the DESeq2 cubic-spline model in Docker
3. extracts fitted trajectories in Docker
4. generates the heatmap plots in Python
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

def find_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here.parent, here.parent.parent):
        if (candidate / "resources").exists():
            return candidate
    return here.parent.parent


ROOT = find_root()
DOCKER_DIR = ROOT / "resources" / "docker_cubic_splines"
INPUT_DIR = ROOT / "Figure7" / "cardio_deseq2_cubicsplines_input"
OUTPUT_DIR = ROOT / "Figure7" / "cardio_deseq2_cubicsplines_output"
SCRIPTS_DIR = Path(__file__).parent
IMAGE = "docker_cubic_splines-deseq2-cubicsplines"


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("=== Step 1: Prepare spline-model inputs ===")
    run([sys.executable, str(SCRIPTS_DIR / "fig7c_prepare_inputs.py")])

    print("\n=== Step 2: Build Docker image ===")
    run(["docker", "compose", "build"], cwd=DOCKER_DIR)

    print("\n=== Step 3: Run DESeq2 cubic-spline analysis ===")
    run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{INPUT_DIR}:/workspace/input:ro",
            "-v",
            f"{OUTPUT_DIR}:/workspace/output:rw",
            IMAGE,
            "Rscript",
            "/workspace/input/run_deseq2_analysis.R",
        ]
    )

    print("\n=== Step 4: Extract fitted trajectories ===")
    run(
        [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{OUTPUT_DIR}:/workspace/output:rw",
            "-v",
            f"{SCRIPTS_DIR}:/workspace/scripts:ro",
            IMAGE,
            "Rscript",
            "/workspace/scripts/fig7c_extract_curves.R",
        ]
    )

    print("\n=== Step 5: Plot heatmaps ===")
    run([sys.executable, str(SCRIPTS_DIR / "fig7c.py")])
    print("\nDone.")


if __name__ == "__main__":
    main()
