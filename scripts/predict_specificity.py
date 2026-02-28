import os
import pandas as pd

def main():
    lsr_clusters_file = snakemake.input.lsr_clusters
    target_clusters_file = snakemake.input.target_clusters
    
    out_file = snakemake.output.specificity
    site_specific_file = snakemake.output.site_specific
    multi_targeting_file = snakemake.output.multi_targeting

    # 1. Load the MMseqs2 cluster TSVs
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
        # Create empty files to satisfy Snakemake
        for f in [out_file, site_specific_file, multi_targeting_file]:
            with open(f, 'w') as file:
                file.write("lsr_id\tspecificity\n")
        return

    # 2. Merge target gene clusters with the LSR clusters
    df_merged = pd.merge(df_lsr, df_target, on='lsr_id', how='left')

    # 3. Calculate Specificity
    specificity_map = {}
    
    for lsr_clust, group in df_merged.groupby('lsr_cluster'):
        unique_targets = group['target_cluster'].dropna().unique()
        
        if len(unique_targets) == 1:
            specificity_map[lsr_clust] = "Site-Specific"
        elif len(unique_targets) > 1:
            specificity_map[lsr_clust] = "Multi-Targeting"
        else:
            specificity_map[lsr_clust] = "Unknown"

    # 4. Apply the prediction back to individual LSR IDs
    df_merged['specificity'] = df_merged['lsr_cluster'].map(specificity_map)
    df_merged['specificity'] = df_merged['specificity'].fillna("Unknown")

    # 5. Format and write the main output
    df_out = df_merged[['lsr_id', 'specificity']].drop_duplicates()
    df_out.to_csv(out_file, sep='\t', index=False)

    # 6. Filter and write the sub-category outputs
    df_site_specific = df_out[df_out['specificity'] == 'Site-Specific']
    df_site_specific.to_csv(site_specific_file, sep='\t', index=False)

    df_multi = df_out[df_out['specificity'] == 'Multi-Targeting']
    df_multi.to_csv(multi_targeting_file, sep='\t', index=False)

if __name__ == "__main__":
    main()