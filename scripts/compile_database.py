import pandas as pd
import os
import re

def main():
    mge_metadata_path = snakemake.input[0]
    specificity_path = snakemake.input[1]
    clusters_path = snakemake.input[2]
    targeting_path = snakemake.input[3]
    output_path = snakemake.output[0]

    def safe_read(path, is_clust=False):
        if not (os.path.exists(path) and os.path.getsize(path) > 0):
            return pd.DataFrame()
        
        df = pd.read_csv(path, sep='\t')
        
        if is_clust and not df.empty:
            first_col = str(df.columns[0])
            if re.search(r'[_.]', first_col) and not any(h in first_col.lower() for h in ['id', 'cluster', 'member']):
                df = pd.read_csv(path, sep='\t', header=None, names=['cluster', 'lsr_id'])
        return df

    df_merged = safe_read(mge_metadata_path)
    df_spec = safe_read(specificity_path)
    df_clust = safe_read(clusters_path, is_clust=True)
    df_target = safe_read(targeting_path)

    if df_merged.empty:
        pd.DataFrame().to_csv(output_path, sep='\t', index=False)
        return

    # Merge Clustering Results
    if not df_clust.empty:
        potential_id_cols = ['lsr_id', 'member', 'id', 'query', 'representative', 'accession', 'sequence_id']
        id_col = next((c for c in potential_id_cols if c in df_clust.columns), None)
        if id_col is None:
            id_col = df_clust.columns[0]
        
        if id_col in df_clust.columns:
            if df_clust[id_col].astype(str).str.contains('_').any() and 'lsr_id' not in df_clust.columns:
                df_clust['lsr_id'] = df_clust[id_col].apply(lambda x: str(x).rsplit('_', 1)[0])
                merge_key = 'lsr_id'
            else:
                merge_key = id_col
            df_merged = pd.merge(df_merged, df_clust, left_on='lsr_id', right_on=merge_key, how='left', suffixes=('', '_clust'))

    # Merge Specificity & Targeting
    if not df_spec.empty:
        spec_key = next((c for c in ['lsr_id', 'id', 'query'] if c in df_spec.columns), df_spec.columns[0])
        df_merged = pd.merge(df_merged, df_spec, left_on='lsr_id', right_on=spec_key, how='left', suffixes=('', '_spec'))

    if not df_target.empty:
        target_key = next((c for c in ['lsr_id', 'id', 'query'] if c in df_target.columns), df_target.columns[0])
        df_merged = pd.merge(df_merged, df_target, left_on='lsr_id', right_on=target_key, how='left', suffixes=('', '_target'))

    # Cleanup and Output
    cols_priority = [
        'source_genome', 'lsr_id', 'pair_id', 'cluster', 'cluster_id',
        'predicted_specificity', 'mge_size', 'core_homology', 'lsr_dna_seq', 
        'attL_seq', 'attR_seq'
    ]
    
    if 'cluster' not in df_merged.columns:
        for c in ['cluster_id', 'rep', 'representative']:
            if c in df_merged.columns:
                df_merged = df_merged.rename(columns={c: 'cluster'})
                break

    existing_cols = [c for c in cols_priority if c in df_merged.columns]
    other_cols = [c for c in df_merged.columns if c not in existing_cols]
    
    final_cols = []
    seen = set()
    for c in (existing_cols + other_cols):
        if c not in seen and not any(c.endswith(s) for s in ['_clust', '_spec', '_target']):
            final_cols.append(c)
            seen.add(c)

    df_merged[final_cols].to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    main()