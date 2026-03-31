"""
Figure 7c upstream input generation.

Recreates the notebook-generated intermediate files:
- Figure7/cardio_deseq2_cubicsplines_input/count_matrix.csv
- Figure7/cardio_deseq2_cubicsplines_input/sample_metadata.csv

Source notebook:
- Figure7/DESeq2_Cardio_Cubic_Splines.ipynb
"""

from pathlib import Path

from fig7_shared import ROOT, build_cardio_count_matrix

OUT_DIR = ROOT / "Figure7" / "cardio_deseq2_cubicsplines_input"
OUT_DIR.mkdir(exist_ok=True)


def main() -> None:
    print("Preparing Figure 7c DESeq2 spline inputs...")
    count_matrix, sample_metadata, _ = build_cardio_count_matrix(include_day10=False)

    count_path = OUT_DIR / "count_matrix.csv"
    metadata_path = OUT_DIR / "sample_metadata.csv"

    count_matrix.to_csv(count_path)
    sample_metadata.to_csv(metadata_path, index=False)

    print(f"Saved: {count_path}")
    print(f"  Genes: {count_matrix.shape[0]}, Samples: {count_matrix.shape[1]}")
    print(f"Saved: {metadata_path}")
    print(f"  Timepoints: {sample_metadata['Time'].nunique()}")
    print(f"  Replicates per timepoint: {sorted(sample_metadata.groupby('Time').size().unique())}")


if __name__ == "__main__":
    main()
