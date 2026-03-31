"""
Rediscovery Analysis Library
Quantile-based rediscovery analysis for gene expression data
"""

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from itertools import combinations


def calculate_rediscovery_curve(source_samples, target_samples, source_data,
                                direction='bottom', threshold=1.0, n_quantiles=50,
                                detection_threshold=0):
    """
    Calculate rediscovery rate across quantiles using gene ranking.

    Parameters:
    -----------
    source_samples : list
        List of sample column names to use for gene ranking
    target_samples : list
        List of sample column names to check for gene detection
    source_data : pd.DataFrame
        Gene expression dataframe with sample columns (no _TPM suffix needed)
    direction : str
        'bottom' (low expression) or 'top' (high expression)
    threshold : float
        Minimum mean TPM to include genes
    n_quantiles : int
        Number of quantile points to calculate
    detection_threshold : float
        TPM threshold for detection (genes > this value count as detected)

    Returns:
    --------
    tuple : (quantiles, rediscovery_rates, auc, gene_count)
        - quantiles: array of quantile percentages (0-100)
        - rediscovery_rates: array of rediscovery rates (0-100)
        - auc: area under curve
        - gene_count: total genes in analysis
    """

    # Get column names - check if they exist as-is or need _TPM suffix
    def get_valid_cols(samples, df):
        """Get column names that exist in dataframe"""
        valid = []
        for s in samples:
            if s in df.columns:
                valid.append(s)
            elif f'{s}_TPM' in df.columns:
                valid.append(f'{s}_TPM')
        return valid

    source_tpm_cols = get_valid_cols(source_samples, source_data)
    target_tpm_cols = get_valid_cols(target_samples, source_data)

    if len(source_tpm_cols) == 0:
        raise ValueError(f"No columns found for source samples: {source_samples}. Available: {source_data.columns.tolist()[:10]}")

    if len(target_tpm_cols) == 0:
        raise ValueError(f"No columns found for target samples: {target_samples}")

    # Calculate mean expression in source samples
    mean_expr = source_data[source_tpm_cols].mean(axis=1)

    # Filter genes by threshold
    mask = mean_expr > threshold
    filtered_data = source_data.loc[mask, source_tpm_cols + target_tpm_cols].copy()
    filtered_mean_expr = mean_expr[mask]

    gene_count = len(filtered_data)

    if gene_count == 0:
        return np.array([0, 100]), np.array([0, 0]), 0, 0

    # Sort by expression (ascending)
    sort_idx = np.argsort(filtered_mean_expr.values)
    sorted_indices = filtered_data.index[sort_idx]
    sorted_data = source_data.loc[sorted_indices, source_tpm_cols + target_tpm_cols].copy()

    # Calculate quantiles – use equal-width 2% bins starting at 2%
    # so every bin contains roughly the same number of genes
    step = 100 / n_quantiles
    quantile_points = np.arange(step, 100 + step / 2, step)  # e.g. [2, 4, ..., 100]
    rediscovery_rates = []

    for quantile in quantile_points:
        n_genes = max(1, int(np.round(gene_count * quantile / 100)))

        if direction == 'bottom':
            # Bottom quantile: genes with LOWEST expression
            genes_in_quantile = sorted_data.iloc[:n_genes]
        elif direction == 'top':
            # Top quantile: genes with HIGHEST expression
            genes_in_quantile = sorted_data.iloc[-n_genes:]
        else:
            raise ValueError(f"direction must be 'bottom' or 'top', got {direction}")

        # Check detection in ALL target replicates
        detected_in_all = (genes_in_quantile[target_tpm_cols] > detection_threshold).all(axis=1).sum()

        rediscovery_rate = 100 * (detected_in_all / len(genes_in_quantile))
        rediscovery_rates.append(rediscovery_rate)

    rediscovery_rates = np.array(rediscovery_rates)
    auc = trapezoid(rediscovery_rates, quantile_points)

    return quantile_points, rediscovery_rates, auc, gene_count


def reorganize_samples_by_number(samples, mapping):
    """
    Reorganize flat sample list into nested dict by sample number.

    Parameters:
    -----------
    samples : list
        List of sample column names
    mapping : dict
        Sample name to metadata mapping dict with keys: Sample, Replicate, Type, Experiment

    Returns:
    --------
    dict : {sample_num: [rep1_name, rep2_name, ...]}
        Samples organized by number with sorted replicates
    """
    organized = {}
    for sample_name in samples:
        if sample_name in mapping:
            sample_num = mapping[sample_name].get('Sample')
            if sample_num not in organized:
                organized[sample_num] = []
            organized[sample_num].append(sample_name)

    # Sort replicates within each sample
    for sample_num in organized:
        organized[sample_num] = sorted(organized[sample_num],
                                      key=lambda s: mapping[s].get('Replicate', 0))

    return organized


def run_per_sample_analysis(organized_samples, gene_df, sample_mapping,
                            source_type='L', target_type='SN'):
    """
    Run per-sample rediscovery analysis with local ranking.

    Parameters:
    -----------
    organized_samples : dict
        {sample_type: {sample_num: [rep1, rep2, ...]}}
    gene_df : pd.DataFrame
        Gene expression dataframe
    sample_mapping : dict
        Sample name to metadata mapping
    source_type : str
        Source sample type ('L', 'SN', or 'Lysate')
    target_type : str
        Target sample type ('L', 'SN', or 'Lysate')

    Returns:
    --------
    dict : {sample_num: {analysis_key: {quantiles, rates, auc, gene_count}}}
    """

    analysis_types = [
        (source_type, target_type, 'bottom'),
        (source_type, target_type, 'top'),
        (target_type, source_type, 'bottom'),
        (target_type, source_type, 'top'),
    ]

    results = {}

    for sample_num in sorted(organized_samples[source_type].keys()):
        source_reps = organized_samples[source_type][sample_num]
        target_reps = organized_samples[target_type][sample_num]

        results[sample_num] = {}

        for src, tgt, direction in analysis_types:
            analysis_key = f"{src}→{tgt}_{direction}"

            if src == source_type:
                reps_for_ranking = source_reps
                reps_for_testing = target_reps
            else:
                reps_for_ranking = target_reps
                reps_for_testing = source_reps

            quantiles, rates, auc, gene_count = calculate_rediscovery_curve(
                reps_for_ranking, reps_for_testing, gene_df,
                direction=direction, threshold=1.0, n_quantiles=50, detection_threshold=0
            )

            results[sample_num][analysis_key] = {
                'quantiles': quantiles,
                'rates': rates,
                'auc': auc,
                'gene_count': gene_count
            }

    return results


def run_global_ranking_analysis(all_samples, test_samples, gene_df, sample_mapping,
                                source_type='L', target_type='SN', organized_samples=None):
    """
    Run global ranking rediscovery analysis.

    Parameters:
    -----------
    all_samples : dict
        {sample_type: [all sample names]}
    test_samples : dict
        {sample_type: {sample_num: [rep1, rep2, ...]}} - subset to test
    gene_df : pd.DataFrame
        Gene expression dataframe
    sample_mapping : dict
        Sample name to metadata mapping
    source_type : str
        Source sample type for ranking
    target_type : str
        Target sample type for ranking
    organized_samples : dict, optional
        Pre-organized samples dict

    Returns:
    --------
    dict : {sample_num: {analysis_key: {quantiles, rates, auc, gene_count}}}
    """

    analysis_types = [
        (source_type, target_type, 'bottom'),
        (source_type, target_type, 'top'),
        (target_type, source_type, 'bottom'),
        (target_type, source_type, 'top'),
    ]

    results = {}

    for sample_num in sorted(test_samples[source_type].keys()):
        source_test_reps = test_samples[source_type][sample_num]
        target_test_reps = test_samples[target_type][sample_num]

        results[sample_num] = {}

        for src, tgt, direction in analysis_types:
            analysis_key = f"{src}→{tgt}_{direction}"

            # Use ALL samples for ranking, but test subset for detection
            if src == source_type:
                reps_for_ranking = all_samples[source_type]
                reps_for_testing = target_test_reps
            else:
                reps_for_ranking = all_samples[target_type]
                reps_for_testing = source_test_reps

            quantiles, rates, auc, gene_count = calculate_rediscovery_curve(
                reps_for_ranking, reps_for_testing, gene_df,
                direction=direction, threshold=1.0, n_quantiles=50, detection_threshold=0
            )

            results[sample_num][analysis_key] = {
                'quantiles': quantiles,
                'rates': rates,
                'auc': auc,
                'gene_count': gene_count
            }

    return results


def run_replicate_combination_analysis(all_samples, test_samples, gene_df, sample_mapping,
                                       source_type='L', target_type='SN'):
    """
    Run analysis with all possible 2-replicate combinations.

    Parameters:
    -----------
    all_samples : dict
        {sample_type: [all sample names]} - used for global ranking
    test_samples : dict
        {sample_type: {sample_num: [rep1, rep2, ...]}} - samples to test
    gene_df : pd.DataFrame
        Gene expression dataframe
    sample_mapping : dict
        Sample name to metadata mapping
    source_type : str
        Source sample type
    target_type : str
        Target sample type

    Returns:
    --------
    dict : {sample_num: {analysis_key: {quantiles, rates_mean, rates_std, auc_mean, auc_std, gene_count}}}
    """

    analysis_types = [
        (source_type, target_type, 'bottom'),
        (source_type, target_type, 'top'),
        (target_type, source_type, 'bottom'),
        (target_type, source_type, 'top'),
    ]

    results = {}

    for sample_num in sorted(test_samples[source_type].keys()):
        source_reps_all = test_samples[source_type][sample_num]
        target_reps_all = test_samples[target_type][sample_num]

        source_combinations = list(combinations(range(len(source_reps_all)), 2))
        target_combinations = list(combinations(range(len(target_reps_all)), 2))

        results[sample_num] = {}

        for src, tgt, direction in analysis_types:
            analysis_key = f"{src}→{tgt}_{direction}"

            rates_list = []
            auc_list = []
            quantiles_first = None
            gene_count_first = None

            if src == source_type:
                src_combinations = source_combinations
                tgt_combinations = target_combinations
                src_reps_all_local = source_reps_all
                tgt_reps_all_local = target_reps_all
                ranking_samples = all_samples[source_type]
            else:
                src_combinations = target_combinations
                tgt_combinations = source_combinations
                src_reps_all_local = target_reps_all
                tgt_reps_all_local = source_reps_all
                ranking_samples = all_samples[target_type]

            for src_idx_tuple, tgt_idx_tuple in zip(src_combinations, tgt_combinations):
                src_reps = [src_reps_all_local[i] for i in src_idx_tuple]
                tgt_reps = [tgt_reps_all_local[i] for i in tgt_idx_tuple]

                quantiles, rates, auc, gene_count = calculate_rediscovery_curve(
                    ranking_samples, tgt_reps, gene_df,
                    direction=direction, threshold=1.0, n_quantiles=50, detection_threshold=0
                )

                rates_list.append(rates)
                auc_list.append(auc)

                if quantiles_first is None:
                    quantiles_first = quantiles
                    gene_count_first = gene_count

            rates_array = np.array(rates_list)
            rates_mean = np.mean(rates_array, axis=0)
            rates_std = np.std(rates_array, axis=0)
            auc_mean = np.mean(auc_list)
            auc_std = np.std(auc_list)

            results[sample_num][analysis_key] = {
                'quantiles': quantiles_first,
                'rates_mean': rates_mean,
                'rates_std': rates_std,
                'auc_mean': auc_mean,
                'auc_std': auc_std,
                'gene_count': gene_count_first
            }

    return results


def create_auc_comparison_table(results):
    """
    Create AUC comparison dataframe from results.

    Parameters:
    -----------
    results : dict
        {sample_num: {analysis_key: {..., 'auc': value}}} or {..., 'auc_mean': value}

    Returns:
    --------
    pd.DataFrame : AUC comparison matrix
    """
    auc_matrix = {}

    for sample_num in sorted(results.keys()):
        auc_matrix[sample_num] = {}
        for analysis_key in results[sample_num].keys():
            if 'auc_mean' in results[sample_num][analysis_key]:
                auc = results[sample_num][analysis_key]['auc_mean']
            else:
                auc = results[sample_num][analysis_key]['auc']
            auc_matrix[sample_num][analysis_key] = auc

    return pd.DataFrame(auc_matrix).T.round(2)


def print_auc_rankings(comparison_df, analysis_name="Analysis"):
    """
    Print AUC rankings for each analysis method.

    Parameters:
    -----------
    comparison_df : pd.DataFrame
        AUC comparison matrix
    analysis_name : str
        Name of the analysis for printing
    """
    print(f"\n{'='*100}")
    print(f"RANKINGS BY METHOD ({analysis_name})")
    print(f"{'='*100}")

    for analysis_key in comparison_df.columns:
        sorted_samples = comparison_df[analysis_key].sort_values(ascending=False)
        print(f"\n{analysis_key}:")
        print("  Rank | Sample | AUC")
        print("  -----|--------|------")
        for rank, (sample, auc) in enumerate(sorted_samples.items(), 1):
            print(f"   {rank}   |   {sample}   | {auc:.2f}")
