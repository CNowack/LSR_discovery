import os
import pandas as pd

def main():
    # 1. Access inputs/outputs/params
    mge_files = snakemake.input.mge_data
    hmmer_files = snakemake.input.hmmer_data
    att_files = snakemake.input.att_data
    
    out_tsv = snakemake.output.filtered
    out_fasta = snakemake.output.fasta

    MAX_MGE_SIZE = snakemake.params.max_mge_kb * 1000

    filtered_candidates = []

    # 2. Map files to their respective genomes
    genome_map = {}
    for f in mge_files:
        gid = f.split('/')[-2]
        genome_map.setdefault(gid, {})['mge'] = f
    
    for f in hmmer_files:
        gid = os.path.basename(f).replace('.tblout', '')
        genome_map.setdefault(gid, {})['hmmer'] = f

    for f in att_files:
        gid = os.path.basename(f).replace('.tsv', '')
        genome_map.setdefault(gid, {})['att'] = f

    # 3. Process each genome
    for genome_id, files in genome_map.items():
        mge_path = files.get('mge')
        hmmer_path = files.get('hmmer')
        att_path = files.get('att')
        
        if not all([mge_path, hmmer_path, att_path]) or not all(os.path.exists(p) for p in [mge_path, hmmer_path, att_path]):
            continue
            
        try:
            # Load MGE and ATT files
            df_mge = pd.read_csv(mge_path, sep='\t')
            df_att = pd.read_csv(att_path, sep='\t')
            
            # Load HMMER file using regex whitespace separator
            df_hmmer = pd.read_csv(hmmer_path, sep=r'\s+', comment='#', header=None)
        except Exception:
            continue

        if df_hmmer.empty or df_mge.empty or df_att.empty:
            continue

        # 4. Filter and Merge Data
        for _, mge_row in df_mge.iterrows():
            # Extract MGE size safely
            mge_size = mge_row.get('inferred_seq_length', 0)
            if pd.isna(mge_size) or mge_size > MAX_MGE_SIZE:
                continue
            
            pair_id = mge_row.get('pair_id')
            
            # Fetch corresponding attachment sites
            att_match = df_att[df_att['pair_id'] == pair_id]
            attL_seq = att_match['attL_seq'].values[0] if not att_match.empty else "NOT_FOUND"
            attR_seq = att_match['attR_seq'].values[0] if not att_match.empty else "NOT_FOUND"

            # Cross-reference with HMMER hits
            for _, hmm_hit in df_hmmer.iterrows():
                protein_id = hmm_hit[0] # The PROKKA or NCBI protein ID
                
                filtered_candidates.append({
                    'source_genome': genome_id,
                    'lsr_id': protein_id,
                    'pair_id': pair_id,
                    'mge_size': mge_size,
                    'attL_seq': attL_seq,
                    'attR_seq': attR_seq
                })

    # 5. Write outputs
    df_out = pd.DataFrame(filtered_candidates)
    
    if not df_out.empty:
        df_out.to_csv(out_tsv, sep='\t', index=False)
        with open(out_fasta, 'w') as f:
            for _, row in df_out.iterrows():
                f.write(f">{row['lsr_id']}\nSEQUENCE_PLACEHOLDER\n")
    else:
        # Failsafe for empty output
        pd.DataFrame(columns=['source_genome', 'lsr_id', 'pair_id', 'mge_size', 'attL_seq', 'attR_seq']).to_csv(out_tsv, sep='\t', index=False)
        open(out_fasta, 'w').close()

if __name__ == "__main__":
    main()