def calculate_enrichment_factors(
    gene_id_is_mt,
    gene_id_to_biotype,
    count_df,
    sn_list,
    lysate_list,
    pseudocount=1,
    verbose=False,
):
    import numpy as np

    ensg_df = count_df["rpm_ensg_only"]
    sn_mean = ensg_df[sn_list].mean(axis=1) + pseudocount
    lysate_mean = ensg_df[lysate_list].mean(axis=1) + pseudocount
    all_enrichment_factors = sn_mean / lysate_mean

    protein_coding_genes = [
        gene_id for gene_id in all_enrichment_factors.index if gene_id_to_biotype.get(gene_id) == "protein_coding"
    ]
    protein_coding_enrichment_factors = all_enrichment_factors[protein_coding_genes]

    mt_genes = [gene_id for gene_id in all_enrichment_factors.index if gene_id_is_mt.get(gene_id, False) is True]
    mt_enrichment_factors = all_enrichment_factors[mt_genes]

    if verbose:
        print(f"Total genes: {len(all_enrichment_factors)}")
        print(f"Protein coding genes: {len(protein_coding_enrichment_factors)}")
        print(f"Mitochondrial genes: {len(mt_enrichment_factors)}")
        print(f"Mean enrichment factor (all genes): {all_enrichment_factors.mean():.3f}")
        print(f"Median enrichment factor (all genes): {all_enrichment_factors.median():.3f}")
        if np.isinf(all_enrichment_factors).any():
            print("Warning: Found infinite values in enrichment factors")
        if all_enrichment_factors.isna().any():
            print("Warning: Found NaN values in enrichment factors")

    return {
        "all_enrichment_factors": all_enrichment_factors,
        "protein_coding_enrichment_factors": protein_coding_enrichment_factors,
        "mt_enrichment_factors": mt_enrichment_factors,
    }
