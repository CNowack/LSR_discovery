import os
import pandas as pd

def main():
    filtered_file = snakemake.input.filtered
    spec_file = snakemake.input.specificity
    clust_file = snakemake.input.clusters_50
    gt_file = snakemake.input.genome_targeting
    out_db = snakemake.output.database

    # Define the final schema for the database
    final_cols = [
        "lsr_id", "source_genome", "genome_label", "contig", "lsr_start", "lsr_end", 
        "distance_to_att", "mge_size", "cluster_50", "specificity",
        "attA_id", "human_chr", "human_start", "human_end", "pident"
    ]

    # 1. Load Filtered LSR Candidates
    try:
        df_lsr = pd.read_csv(filtered_file, sep='\t')
    except pd.errors.EmptyDataError:
        df_lsr = pd.DataFrame()

    # Failsafe: If no candidates survived filtering, write an empty DB and exit
    if df_lsr.empty:
        pd.DataFrame(columns=final_cols).to_csv(out_db, sep='\t', index=False)
        return

    # Standardize the genome column for merging
    df_lsr = df_lsr.rename(columns={"genome": "genome_label"})

    # Map the labels back to accessions for the validation script
    test_csv = "resources/test_genomes.csv"
    label_to_acc = {}
    if os.path.exists(test_csv):
        try:
            df_test = pd.read_csv(test_csv)
            if 'label' in df_test.columns and 'accession' in df_test.columns:
                label_to_acc = dict(zip(df_test['label'], df_test['accession']))
        except Exception:
            pass
    
    df_lsr['source_genome'] = df_lsr['genome_label'].map(lambda x: label_to_acc.get(x, x))

    # 2. Merge Specificity Predictions
    try:
        df_spec = pd.read_csv(spec_file, sep='\t')
        if not df_spec.empty and 'lsr_id' in df_spec.columns:
            df_lsr = pd.merge(df_lsr, df_spec[['lsr_id', 'specificity']], on='lsr_id', how='left')
        else:
            df_lsr['specificity'] = "Unknown"
    except (pd.errors.EmptyDataError, FileNotFoundError):
        df_lsr['specificity'] = "Unknown"

    # 3. Merge MMseqs2 Clusters
    try:
        # MMseqs2 TSV format is typically headerless: [cluster_rep, member]
        df_clust = pd.read_csv(clust_file, sep='\t', header=None, names=['cluster_50', 'lsr_id'])
        if not df_clust.empty:
            df_lsr = pd.merge(df_lsr, df_clust, on='lsr_id', how='left')
        else:
            df_lsr['cluster_50'] = df_lsr['lsr_id']
    except (pd.errors.EmptyDataError, FileNotFoundError):
        df_lsr['cluster_50'] = df_lsr['lsr_id']

    # 4. Merge Genome Targeting (BLAST) Hits
    try:
        df_gt = pd.read_csv(gt_file, sep='\t')
        if not df_gt.empty:
            df_gt = df_gt.rename(columns={"genome": "genome_label"})
            
            # Left join on genome_label. If a genome has multiple human targets,
            # this will safely duplicate the LSR row to capture each targeting event.
            merge_cols = ['genome_label', 'attA_id', 'human_chr', 'human_start', 'human_end', 'pident']
            existing_gt_cols = [c for c in merge_cols if c in df_gt.columns]
            
            df_lsr = pd.merge(df_lsr, df_gt[existing_gt_cols], on='genome_label', how='left')
        else:
            for c in ['attA_id', 'human_chr', 'human_start', 'human_end', 'pident']:
                df_lsr[c] = None
    except (pd.errors.EmptyDataError, FileNotFoundError):
        for c in ['attA_id', 'human_chr', 'human_start', 'human_end', 'pident']:
            df_lsr[c] = None

    # Ensure schema enforcement and handle missing values
    for col in final_cols:
        if col not in df_lsr.columns:
            df_lsr[col] = None

    df_lsr = df_lsr[final_cols]
    
    # Write the compiled database
    df_lsr.to_csv(out_db, sep='\t', index=False)

if __name__ == "__main__":
    main()