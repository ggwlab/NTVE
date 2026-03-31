import hashlib
import os
import pickle
from pathlib import Path

import pandas as pd


def get_file_hash(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def load_gtf_df(path_to_gtf=None, cache_dir="/tmp/gtf_cache", use_cache=True):
    if path_to_gtf is None:
        env_path = os.environ.get("NTVE_GTF_PATH")
        if env_path:
            path_to_gtf = env_path
        else:
            default_repo_path = (
                Path(__file__).resolve().parents[1]
                / "resources"
                / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"
            )
            path_to_gtf = str(default_repo_path)

    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)

    base_name = Path(path_to_gtf).stem.replace(".csv", "")
    cache_file = cache_path / f"{base_name}_processed.pkl"
    hash_file = cache_path / f"{base_name}_hash.txt"

    if use_cache and cache_file.exists() and hash_file.exists():
        try:
            current_hash = get_file_hash(path_to_gtf)
            with open(hash_file, "r") as f:
                cached_hash = f.read().strip()
            if current_hash == cached_hash:
                print(f"Loading cached GTF data from {cache_file}")
                with open(cache_file, "rb") as f:
                    cached_data = pickle.load(f)
                print(
                    f"Loaded cached GTF data: {len(cached_data['gene_id_to_name'])} gene names, "
                    f"{len(cached_data['gene_id_to_biotype'])} biotypes, "
                    f"{sum(cached_data['gene_id_is_mt'].values())} mitochondrial genes"
                )
                return cached_data
            print("Source file has changed, rebuilding cache...")
        except Exception as e:
            print(f"Error loading cache: {e}, rebuilding...")

    try:
        print(f"Loading GTF data from {path_to_gtf}")
        gtf_df = pd.read_csv(path_to_gtf)
        df_unique = gtf_df.drop_duplicates("gene_id").set_index("gene_id")

        gene_id_to_name = {
            gene_id: gene_name if pd.notna(gene_name) else gene_id
            for gene_id, gene_name in df_unique["gene_name"].items()
        }
        gene_id_to_biotype = {
            gene_id: biotype if pd.notna(biotype) else gene_id
            for gene_id, biotype in df_unique["gene_biotype"].items()
        }
        gene_id_is_mt = {
            gene_id: bool(pd.notna(gene_name) and str(gene_name).startswith("MT-"))
            for gene_id, gene_name in df_unique["gene_name"].items()
        }

        result = {
            "gtf_df": gtf_df,
            "gene_id_to_name": gene_id_to_name,
            "gene_id_to_biotype": gene_id_to_biotype,
            "gene_id_is_mt": gene_id_is_mt,
        }

        if use_cache:
            try:
                print(f"Caching processed data to {cache_file}")
                with open(cache_file, "wb") as f:
                    pickle.dump(result, f)
                current_hash = get_file_hash(path_to_gtf)
                with open(hash_file, "w") as f:
                    f.write(current_hash)
                print("Data cached successfully")
            except Exception as e:
                print(f"Warning: Could not cache data: {e}")

        mt_genes_count = sum(gene_id_is_mt.values())
        print(
            f"Loaded GTF data: {len(gene_id_to_name)} gene names, "
            f"{len(gene_id_to_biotype)} biotypes, "
            f"{mt_genes_count} mitochondrial genes"
        )
    except Exception as e:
        print(f"Error: Could not load GTF file: {e}")
        return {
            "gtf_df": None,
            "gene_id_to_name": {},
            "gene_id_to_biotype": {},
            "gene_id_is_mt": {},
        }

    return result
