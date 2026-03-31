import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_and_clean_featurecounts(counts_file, suffixes_to_remove=None, dataset_name="Dataset", verbose=False):
    if verbose:
        print(f"Loading FeatureCounts Gene Expression Data: {dataset_name}")
        print("=" * 60)

    try:
        raw_counts_df = pd.read_csv(counts_file, sep="\t", skiprows=1)
        if verbose:
            print(f"Successfully loaded data from {counts_file}")
    except Exception as e:
        if verbose:
            print(f"Error loading file: {e}")
        raise e

    gene_info_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    count_cols = [col for col in raw_counts_df.columns if col not in gene_info_cols]

    clean_sample_names = []
    for col in count_cols:
        sample_id = col.split("/")[-1]
        if not suffixes_to_remove:
            suffixes_to_remove = ["_untrimmed_harmonizedAligned.sortedByCoord.out.bam"]
        for suffix in suffixes_to_remove:
            if sample_id.endswith(suffix):
                sample_id = sample_id.replace(suffix, "")
                break
        clean_sample_names.append(sample_id)

    if len(clean_sample_names) != len(set(clean_sample_names)):
        clean_sample_names_fixed = []
        name_counts = {}
        for name in clean_sample_names:
            if name in name_counts:
                name_counts[name] += 1
                clean_sample_names_fixed.append(f"{name}_dup{name_counts[name]}")
            else:
                name_counts[name] = 0
                clean_sample_names_fixed.append(name)
        clean_sample_names = clean_sample_names_fixed

    total_counts_df = raw_counts_df.copy()
    total_counts_df.columns = gene_info_cols + clean_sample_names
    total_counts_df = total_counts_df.set_index("Geneid")

    count_data_only = total_counts_df[clean_sample_names]
    total_mapped_reads = count_data_only.sum(axis=0)
    ensg_mask = total_counts_df.index.str.startswith("ENSG", na=False)

    return {
        "total_counts": total_counts_df,
        "sample_names": clean_sample_names,
        "gene_info_cols": gene_info_cols,
        "summary_stats": {
            "total_genes": len(total_counts_df),
            "ensg_genes": int(ensg_mask.sum()),
            "non_ensg_genes": int(len(total_counts_df) - ensg_mask.sum()),
            "total_samples": len(clean_sample_names),
            "min_count": count_data_only.values.min(),
            "max_count": count_data_only.values.max(),
            "total_mapped_reads_per_sample": total_mapped_reads.to_dict(),
        },
    }


def calculate_rpm_normalization(counts_data, ensg_only=True, verbose=False):
    total_counts_df = counts_data["total_counts"]
    sample_names = counts_data["sample_names"]
    count_data_only = total_counts_df[sample_names]
    total_mapped_reads = count_data_only.sum(axis=0)

    rpm_df = total_counts_df.copy()
    for sample in sample_names:
        rpm_df[sample] = (count_data_only[sample] / total_mapped_reads[sample]) * 1e6

    result = {
        "rpm": rpm_df,
        "summary_stats": {"total_mapped_reads_per_sample": total_mapped_reads.to_dict()},
    }

    if ensg_only:
        ensg_mask = rpm_df.index.str.startswith("ENSG", na=False)
        ensg_genes = rpm_df[ensg_mask]
        rpm_renormalized_df = ensg_genes.copy()
        ensg_rpm_totals = ensg_genes[sample_names].sum(axis=0)
        for sample in sample_names:
            rpm_renormalized_df[sample] = (ensg_genes[sample] / ensg_rpm_totals[sample]) * 1e6
        result["rpm_ensg_only"] = rpm_renormalized_df
        result["ensg_genes_count"] = len(ensg_genes)
        result["summary_stats"]["ensg_rpm_totals_per_sample"] = ensg_rpm_totals.to_dict()
        result["summary_stats"]["total_ensg_genes"] = len(ensg_genes)
        result["summary_stats"]["total_non_ensg_genes"] = len(rpm_df) - len(ensg_genes)

    return result


def load_and_preprocess_featurecounts(
    counts_file,
    sample_mapping_dict=None,
    suffixes_to_remove=None,
    dataset_name="Dataset",
    ensg_only=True,
    verbose=False,
):
    counts_data = load_and_clean_featurecounts(
        counts_file=counts_file,
        suffixes_to_remove=suffixes_to_remove,
        dataset_name=dataset_name,
        verbose=verbose,
    )
    rpm_data = calculate_rpm_normalization(
        counts_data=counts_data,
        ensg_only=ensg_only,
        verbose=verbose,
    )

    result = {
        "total_counts": counts_data["total_counts"],
        "rpm": rpm_data["rpm"],
        "sample_names": counts_data["sample_names"],
        "gene_info_cols": counts_data["gene_info_cols"],
        "summary_stats": {**counts_data["summary_stats"], **rpm_data["summary_stats"]},
    }
    if ensg_only and "rpm_ensg_only" in rpm_data:
        result["rpm_ensg_only"] = rpm_data["rpm_ensg_only"]
    return result
