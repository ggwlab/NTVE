#!/usr/bin/env bash
# Stage a minimal explicit refactoring bundle into /tmp, run it there, and
# sync generated outputs back into refactoring_roadmap.
#
# Run from repo root:
#   bash refactoring_roadmap/run_stage_tmp_complete.sh
#
# This script intentionally does NOT copy the whole repository. It writes an
# explicit manifest into the staged workspace so copied inputs are auditable.

set -uo pipefail

CONDA_ENV="ntve"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
STAGE_ROOT="${TMPDIR:-/tmp}/ntve_refactor_stage_${STAMP}"
STAGED_REPO="${STAGE_ROOT}/repo"
STAGED_LOG_DIR="${STAGE_ROOT}/logs"
MPLCONFIGDIR_STAGE="${STAGE_ROOT}/mplconfig"
PYTHONPATH_STAGE="${STAGED_REPO}:${PYTHONPATH:-}"

PRESERVE_ROOT="${REPO_ROOT}/staged_runs/${STAMP}"
RUN_LOG="${PRESERVE_ROOT}/run.log"
MANIFEST_FILE="${PRESERVE_ROOT}/manifest.txt"

mkdir -p "${STAGED_REPO}" "${STAGED_LOG_DIR}" "${MPLCONFIGDIR_STAGE}" "${PRESERVE_ROOT}"

ok()  { echo "[OK]  $1"; }
err() { echo "[FAIL] $1"; }
note() { echo "[INFO] $1"; }

SCRIPT_FILES=(
  "1e.py"
  "1f.py"
  "1g.py"
  "1j.py"
  "2a.py"
  "2b.py"
  "2c.py"
  "3a.py"
  "3bc.py"
  "3d.py"
  "4cde.py"
  "figure4_prepare_deseq2.py"
  "5b.py"
  "6c.py"
  "figure2_prepare_loaded_data.py"
  "fig7_2vs2.py"
  "fig7_shared.py"
  "fig7c.py"
  "fig7c_extract_curves.R"
  "fig7c_pipeline.py"
  "fig7c_prepare_inputs.py"
  "run_deseq2_analysis.R"
  "rediscovery_analysis_lib.py"
  "refactoring_roadmap_figure2_shared.py"
  "suppl10.py"
  "suppl11a.py"
  "suppl11b.py"
  "suppl11b_2a.py"
  "suppl12a.py"
  "suppl12b.py"
  "suppl13a.py"
  "suppl13b.py"
  "suppl14.py"
  "suppl15a.py"
  "suppl15bc.py"
  "suppl16ab.py"
  "suppl19.py"
  "suppl20a.py"
  "suppl21a.py"
  "suppl21bc.py"
  "suppl21bc_rankings.py"
  "suppl5.py"
  "suppl5d.py"
  "suppl5e.py"
  "suppl6_data.py"
  "suppl6ab.py"
  "suppl8ab.py"
  "ntvetools"
)

REFERENCE_INPUTS=(
  "resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
  "resources/dominant_isoform_genes.csv.gz"
  "resources/harmonized_harmonized_gene_counts_rv_stranded.txt"
  "resources/comparison_harmonized_gene_counts.txt.gz"
  "resources/Sample_barcode.csv"
  "resources/figure4_Project_1716_lims_simplified.csv"
  "resources/Project_1716_lims_simplified_Spam2_deleted.csv"
  "resources/Project_1716_lims_simplified_Spam2_without_missing.csv"
  "resources/c2.all.v2025.1.Hs.json"
  "resources/c5.go.v2025.1.Hs.json"
  "resources/docker_cubic_splines"
  "resources/docker_DeSeq_with_Shrinkage"
  "resources/docker_fgsea"
)

STAGED_INPUT_MAPPINGS=(
  "resources/figure4_Project_1716_lims_simplified.csv|Figure4/Project_1716_lims_simplified.csv"
  "resources/suppl14|Suppl14"
  "resources/suppl20/blast_samples|Suppl20/blast_samples"
  "resources/suppl20/blast_results|Suppl20/blast_results"
  "resources/sup5_comparison_data|Suppl5/comparison_data"
  "resources/sup5_comparison_data|Suppl10/comparison_data"
  "resources/cardio_data.pkl|Suppl21/cardio_data.pkl"
  "resources/fig3/Dox_vs_IFNg_lysate_base_mean_gt_0.csv|Figure3/Dox_vs_IFNg_lysate_base_mean_gt_0.csv"
  "resources/fig3/Dox_vs_IFNg_SN_base_mean_gt_0.csv|Figure3/Dox_vs_IFNg_SN_base_mean_gt_0.csv"
  "resources/fig3/greaterAbs_Dox_vs_IFNg_lysate_base_mean_gt_0.csv|Figure3/greaterAbs_Dox_vs_IFNg_lysate_base_mean_gt_0.csv"
  "resources/fig3/lessAbs_Dox_vs_IFNg_lysate_base_mean_gt_0.csv|Figure3/lessAbs_Dox_vs_IFNg_lysate_base_mean_gt_0.csv"
  "resources/fig3/greaterAbs_Dox_vs_IFNg_SN_base_mean_gt_0.csv|Figure3/greaterAbs_Dox_vs_IFNg_SN_base_mean_gt_0.csv"
  "resources/figure6_three_lineage_data.pkl|Figure6/pipeline/three_lineage_data.pkl"
  "resources/gsea_ifn-gamma/lysate/gsea_report_for_na_pos_1727164776729.tsv|Figure3/IFN-gamma_lysate.GseaPreranked.1727164776729/gsea_report_for_na_pos_1727164776729.tsv"
  "resources/gsea_ifn-gamma/ntve/gsea_report_for_na_pos_1727164795980.tsv|Figure3/IFN-gamma_SN.GseaPreranked.1727164795980/gsea_report_for_na_pos_1727164795980.tsv"
)

PRIMARY_DATA_INPUTS=(
  "resources/salmon_harmonized"
  "resources/quantfiles_filtered_pipeline"
)

PRECOMPUTED_UPSTREAM_INPUTS=(
  "resources/cardio_data.pkl"
  "resources/coverage_arrays"
  "resources/data_dict.joblib"
  "resources/fig3"
  "resources/figure6_three_lineage_data.pkl"
  "resources/gsea_ifn-gamma"
  "resources/suppl14"
  "resources/suppl20"
  "resources/sup5_comparison_data"
  "resources/transcript_coverage.db"
)

RUN_SCRIPTS=(
  "refactoring_roadmap/1e.py"
  "refactoring_roadmap/1f.py"
  "refactoring_roadmap/1g.py"
  "refactoring_roadmap/1j.py"
  "refactoring_roadmap/figure2_prepare_loaded_data.py"
  "refactoring_roadmap/2a.py"
  "refactoring_roadmap/2b.py"
  "refactoring_roadmap/2c.py"
  "refactoring_roadmap/3a.py"
  "refactoring_roadmap/3bc.py"
  "refactoring_roadmap/3d.py"
  "refactoring_roadmap/figure4_prepare_deseq2.py"
  "refactoring_roadmap/4cde.py"
  "refactoring_roadmap/5b.py"
  "refactoring_roadmap/6c.py"
  "refactoring_roadmap/suppl5.py"
  "refactoring_roadmap/suppl5d.py"
  "refactoring_roadmap/suppl5e.py"
  "refactoring_roadmap/suppl6_data.py"
  "refactoring_roadmap/suppl6ab.py"
  "refactoring_roadmap/suppl8ab.py"
  "refactoring_roadmap/suppl10.py"
  "refactoring_roadmap/suppl11a.py"
  "refactoring_roadmap/suppl11b.py"
  "refactoring_roadmap/suppl11b_2a.py"
  "refactoring_roadmap/suppl12a.py"
  "refactoring_roadmap/suppl12b.py"
  "refactoring_roadmap/suppl13a.py"
  "refactoring_roadmap/suppl13b.py"
  "refactoring_roadmap/suppl14.py"
  "refactoring_roadmap/suppl15a.py"
  "refactoring_roadmap/suppl15bc.py"
  "refactoring_roadmap/suppl16ab.py"
  "refactoring_roadmap/fig7c_pipeline.py"
  "refactoring_roadmap/fig7_2vs2.py"
  "refactoring_roadmap/suppl19.py"
  "refactoring_roadmap/suppl20a.py"
  "refactoring_roadmap/suppl21a.py"
  "refactoring_roadmap/suppl21bc_rankings.py"
  "refactoring_roadmap/suppl21bc.py"
)

write_manifest() {
  {
    echo "Staged refactoring run manifest"
    echo "Timestamp: ${STAMP}"
    echo "Repo root: ${REPO_ROOT}"
    echo "Stage root: ${STAGE_ROOT}"
    echo
    echo "Scripts copied:"
    for item in "${SCRIPT_FILES[@]}"; do
      if [ "${item}" = "ntvetools" ] || [ "${item}" = "rediscovery_analysis_lib.py" ]; then
        echo "  - ${item} -> ${item}"
      else
        echo "  - ${item} -> refactoring_roadmap/${item}"
      fi
    done
    echo
    echo "Reference inputs copied:"
    echo "  These are reference assets, metadata, helper code, or Docker build contexts."
    for item in "${REFERENCE_INPUTS[@]}"; do
      echo "  - ${item}"
    done
    echo
    echo "Primary data inputs copied:"
    echo "  These are the main quantification/count tables used as source input."
    for item in "${PRIMARY_DATA_INPUTS[@]}"; do
      echo "  - ${item}"
    done
    echo
    echo "Precomputed upstream inputs copied:"
    echo "  These are NOT raw data. They are cached or notebook-derived inputs that"
    echo "  some roadmap scripts currently depend on and do not yet regenerate."
    for item in "${PRECOMPUTED_UPSTREAM_INPUTS[@]}"; do
      echo "  - ${item}"
    done
    for item in "${STAGED_INPUT_MAPPINGS[@]}"; do
      echo "  - ${item%%|*} -> ${item#*|}"
    done
    echo
    echo "Outputs synced back:"
    echo "  - refactoring_roadmap/*_plots"
    echo "  - refactoring_roadmap/*_output"
    echo "  - refactoring_roadmap/*_csv"
    echo "  - refactoring_roadmap/figure2_loaded_data"
    echo "  - staged_runs/${STAMP}/Figure7/*"
    echo "  - staged_runs/${STAMP}/run.log"
    echo "  - staged_runs/${STAMP}/manifest.txt"
  } > "${MANIFEST_FILE}"
}

require_path() {
  local path="$1"
  if [ ! -e "${REPO_ROOT}/${path}" ]; then
    err "Missing required input: ${path}"
    exit 1
  fi
}

copy_file_to_stage() {
  local src_rel="$1"
  local dst_rel="$2"
  require_path "${src_rel}"
  if [ -d "${REPO_ROOT}/${src_rel}" ]; then
    mkdir -p "${STAGED_REPO}/${dst_rel}"
    rsync -a "${REPO_ROOT}/${src_rel}/" "${STAGED_REPO}/${dst_rel}/"
  else
    mkdir -p "${STAGED_REPO}/$(dirname "${dst_rel}")"
    rsync -a "${REPO_ROOT}/${src_rel}" "${STAGED_REPO}/${dst_rel}"
  fi
}

copy_into_stage() {
  note "Copying explicit manifest into ${STAGED_REPO}"
  mkdir -p "${STAGED_REPO}/refactoring_roadmap"

  local script
  for script in "${SCRIPT_FILES[@]}"; do
    if [ "${script}" = "ntvetools" ] || [ "${script}" = "rediscovery_analysis_lib.py" ]; then
      copy_file_to_stage "${script}" "${script}"
    else
      copy_file_to_stage "${script}" "refactoring_roadmap/${script}"
    fi
  done

  local path
  for path in "${REFERENCE_INPUTS[@]}" "${PRIMARY_DATA_INPUTS[@]}" "${PRECOMPUTED_UPSTREAM_INPUTS[@]}"; do
    copy_file_to_stage "${path}" "${path}"
  done

  local mapping src_rel dst_rel
  for mapping in "${STAGED_INPUT_MAPPINGS[@]}"; do
    src_rel="${mapping%%|*}"
    dst_rel="${mapping#*|}"
    copy_file_to_stage "${src_rel}" "${dst_rel}"
  done
}

run_py() {
  local script="$1"
  echo "" | tee -a "${RUN_LOG}"
  echo ">>> ${script}" | tee -a "${RUN_LOG}"
  if (
    cd "${STAGED_REPO}" && \
    MPLCONFIGDIR="${MPLCONFIGDIR_STAGE}" \
    PYTHONPATH="${PYTHONPATH_STAGE}" \
    conda run --no-capture-output -n "${CONDA_ENV}" python "${script}"
  ) 2>&1 | tee -a "${RUN_LOG}"; then
    ok "${script}" | tee -a "${RUN_LOG}"
    return 0
  else
    err "${script}" | tee -a "${RUN_LOG}"
    return 1
  fi
}

sync_back_outputs() {
  note "Syncing staged refactoring outputs back into repo"
  while IFS= read -r dir; do
    base="$(basename "${dir}")"
    mkdir -p "${REPO_ROOT}/${base}"
    rsync -a "${dir}/" "${REPO_ROOT}/${base}/"
  done < <(find "${STAGED_REPO}/refactoring_roadmap" -maxdepth 1 -type d \( -name '*_plots' -o -name '*_output' -o -name '*_csv' \) | sort)

  if [ -d "${STAGED_REPO}/refactoring_roadmap/figure2_loaded_data" ]; then
    mkdir -p "${REPO_ROOT}/figure2_loaded_data"
    rsync -a \
      "${STAGED_REPO}/refactoring_roadmap/figure2_loaded_data/" \
      "${REPO_ROOT}/figure2_loaded_data/"
  fi

  note "Preserving staged Figure7 pipeline outputs under staged_runs/${STAMP}"
  for dir in \
    "Figure7/cardio_deseq2_cubicsplines_input" \
    "Figure7/cardio_deseq2_cubicsplines_output" \
    "Figure7/deseq2_2vs2_causal_output"
  do
    if [ -e "${STAGED_REPO}/${dir}" ]; then
      mkdir -p "${PRESERVE_ROOT}/$(dirname "${dir}")"
      rsync -a "${STAGED_REPO}/${dir}" "${PRESERVE_ROOT}/$(dirname "${dir}")/"
    fi
  done
}

main() {
  local failures=0

  write_manifest
  copy_into_stage

  {
    echo "======================================================"
    echo " Running staged refactoring bundle in /tmp"
    echo "======================================================"
    echo "Stage root: ${STAGE_ROOT}"
    echo "Outputs preserved under: ${PRESERVE_ROOT}"
    echo
  } | tee -a "${RUN_LOG}"

  for script in "${RUN_SCRIPTS[@]}"; do
    if ! run_py "${script}"; then
      failures=1
    fi
  done

  sync_back_outputs

  {
    echo
    echo "======================================================"
    if [ "${failures}" -eq 0 ]; then
      echo " Staged run finished without script-level failures."
    else
      echo " Staged run finished with one or more failures."
    fi
    echo " Stage root: ${STAGE_ROOT}"
    echo " Preserved logs/manifest: ${PRESERVE_ROOT}"
    echo "======================================================"
  } | tee -a "${RUN_LOG}"

  return "${failures}"
}

main "$@"
