"""
Figure 7c full Python pipeline.

This script:
1. recreates the spline-model input CSVs from quantfiles
2. runs the DESeq2 cubic-spline model in Docker
3. extracts fitted trajectories in Docker
4. generates the heatmap plots in Python
"""

from __future__ import annotations

import argparse
import shutil
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
IMAGE = "ntve-deseq2-cubicsplines:2025-11-07"


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=cwd)


def docker_image_exists(image: str) -> bool:
    result = subprocess.run(
        ["docker", "image", "inspect", image],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Figure 7c spline pipeline.")
    parser.add_argument(
        "--build",
        action="store_true",
        help="Force a Docker rebuild even if the tagged image already exists.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip the Docker build step and require the tagged image to already exist.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Only regenerate the input CSVs and analysis R script.",
    )
    parser.add_argument(
        "--docker-only",
        action="store_true",
        help="Run the Docker DESeq2 and curve-extraction steps, but skip plotting.",
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Only regenerate the Python heatmap from existing trajectory output.",
    )
    args = parser.parse_args()
    if args.build and args.skip_build:
        parser.error("--build and --skip-build cannot be used together")
    if args.prepare_only and (args.docker_only or args.plot_only):
        parser.error("--prepare-only cannot be combined with --docker-only or --plot-only")
    if args.docker_only and args.plot_only:
        parser.error("--docker-only and --plot-only cannot be combined")
    return args


def main() -> None:
    args = parse_args()
    if shutil.which("docker") is None and not (args.prepare_only or args.plot_only):
        raise SystemExit("`docker` is not available in PATH.")

    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    should_prepare = not args.plot_only
    should_run_docker = not (args.prepare_only or args.plot_only)
    should_plot = not (args.prepare_only or args.docker_only)

    if should_prepare:
        print("=== Step 1: Prepare spline-model inputs ===", flush=True)
        run([sys.executable, str(SCRIPTS_DIR / "fig7c_prepare_inputs.py")])
        if args.prepare_only:
            print("\nDone.", flush=True)
            return

    if should_run_docker:
        force_build = args.build
        skip_build = args.skip_build
        image_present = docker_image_exists(IMAGE)

        if force_build:
            print("\n=== Step 2: Build Docker image ===", flush=True)
            run(["docker", "compose", "build", "--progress=plain"], cwd=DOCKER_DIR)
        else:
            print(
                f"\n=== Step 2: Build Docker image ===\nSkipping build; using existing image `{IMAGE}`.",
                flush=True,
            )
            if not image_present:
                raise SystemExit(
                    f"Docker image `{IMAGE}` does not exist. Build it explicitly with `--build`."
                )

        print("\n=== Step 3: Run DESeq2 cubic-spline analysis ===", flush=True)
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

        print("\n=== Step 4: Extract fitted trajectories ===", flush=True)
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

    if should_plot:
        print("\n=== Step 5: Plot heatmaps ===", flush=True)
        run([sys.executable, str(SCRIPTS_DIR / "fig7c.py")])
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
