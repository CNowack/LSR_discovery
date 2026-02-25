import os
import pandas as pd

def main():
    lsr_clusters_file = snakemake.input.lsr_clusters
    target_clusters_file = snakemake.input.target_clusters
    out_file = snakemake.output.specificity

    # 1. Load the MMseqs2 cluster TSVs
    # Format: [cluster_representative, member_id]
    try:
        df_target = pd.read_csv(target_clusters_file, sep='\t', names=['target_cluster', 'lsr_id'])
    except pd.errors.EmptyDataError:
        df_target = pd.DataFrame(columns=['target_cluster', 'lsr_id'])

    try:
        df_lsr = pd.read_csv(lsr_clusters_file, sep='\t', names=['lsr_cluster', 'lsr_id'])
    except pd.errors.EmptyDataError:
        df_lsr = pd.DataFrame(columns=['lsr_cluster', 'lsr_id'])

    # Failsafe: if clustering failed or produced empty results
    if df_target.empty or df_lsr.empty:
        with open(out_file, 'w') as f:
            f.write("lsr_id\tspecificity\n")
        return

    # 2. Merge target gene clusters with the LSR clusters
    df_merged = pd.merge(df_lsr, df_target, on='lsr_id', how='left')

    # 3. Calculate Specificity
    specificity_map = {}
    
    # Group by the LSR family
    for lsr_clust, group in df_merged.groupby('lsr_cluster'):
        # Find all unique target gene families this LSR family integrates into
        unique_targets = group['target_cluster'].dropna().unique()
        
        if len(unique_targets) == 1:
            # All members hit the exact same gene family
            specificity_map[lsr_clust] = "Site-Specific"
        elif len(unique_targets) > 1:
            # Members hit multiple different, unrelated genes
            specificity_map[lsr_clust] = "Multi-Targeting"
        else:
            # No target gene was able to be extracted/clustered
            specificity_map[lsr_clust] = "Unknown"

    # 4. Apply the prediction back to individual LSR IDs
    df_merged['specificity'] = df_merged['lsr_cluster'].map(specificity_map)
    df_merged['specificity'] = df_merged['specificity'].fillna("Unknown")

    # 5. Format and write output
    df_out = df_merged[['lsr_id', 'specificity']].drop_duplicates()
    df_out.to_csv(out_file, sep='\t', index=False)

if __name__ == "__main__":
    main()