import os
import pandas as pd

def main():
    # 1. Access named inputs/outputs/params from the Snakemake object
    # These are lists because of the 'expand' function in the Snakefile
    mge_files = snakemake.input.mge_data
    hmmer_files = snakemake.input.hmmer_data
    att_files = snakemake.input.att_data
    
    out_tsv = snakemake.output.filtered
    out_faa = snakemake.output.fasta

    # 2. Get thresholds from Snakemake params (passed from config.yaml)
    MIN_LSR_LEN = snakemake.params.min_len
    MAX_LSR_LEN = snakemake.params.max_len
    MAX_DIST_TO_ATT = snakemake.params.max_dist
    MAX_MGE_SIZE = snakemake.params.max_mge_kb * 1000

    filtered_candidates = []
    
    print("\n--- STARTING QC DEBUG LOG ---")

    # 3. Create a map to pair the expanded file lists by Genome ID
    # We use the filename to ensure we are comparing the right data for the right genome
    genome_map = {}
    for f in mge_files:
        gid = f.split('/')[-2] # Extract {genome} from .../{genome}/mge_boundaries.tsv
        genome_map.setdefault(gid, {})['mge'] = f
    
    for f in hmmer_files:
        gid = os.path.basename(f).replace('.tblout', '')
        genome_map.setdefault(gid, {})['hmmer'] = f

    # 4. Process each genome
    for genome_id, files in genome_map.items():
        mge_path = files.get('mge')
        hmmer_path = files.get('hmmer')
        
        if not mge_path or not hmmer_path or not os.path.exists(mge_path) or not os.path.exists(hmmer_path):
            continue
            
        try:
            df_mge = pd.read_csv(mge_path, sep='\t')
            # HMMER tblout parsing (ignoring headers)
            df_hmmer = pd.read_csv(hmmer_path, sep='\t', comment='#', header=None, delim_whitespace=True)
        except Exception:
            continue

        if df_hmmer.empty or df_mge.empty:
            continue

        print(f"Evaluating: {genome_id}")

        for _, mge_row in df_mge.iterrows():
            mge_start, mge_end = mge_row['start'], mge_row['end']
            mge_size = abs(mge_end - mge_start)

            if mge_size > MAX_MGE_SIZE:
                print(f"  [REJECTED] MGE size {mge_size} > {MAX_MGE_SIZE}")
                continue

            for _, hmm_hit in df_hmmer.iterrows():
                protein_id = hmm_hit[0]
                # Note: In the full pipeline, coordinates are pulled from the GFF
                # For this debug logic, we verify the HMMER hit exists
                
                # Mock length check for demonstration (replace with actual logic if available)
                prot_len = 500 
                
                if not (MIN_LSR_LEN <= prot_len <= MAX_LSR_LEN):
                    print(f"  [REJECTED] {protein_id}: Length {prot_len} aa")
                    continue

                print(f"  [PASSED] {protein_id}")
                filtered_candidates.append({
                    'source_genome': genome_id,
                    'lsr_id': protein_id,
                    'mge_size': mge_size
                })

    # 5. Write outputs
    df_out = pd.DataFrame(filtered_candidates)
    df_out.to_csv(out_tsv, sep='\t', index=False)
    
    with open(out_faa, 'w') as f:
        if not df_out.empty:
            for _, row in df_out.iterrows():
                f.write(f">{row['lsr_id']}\nSEQUENCE_PLACEHOLDER\n")

    print("--- QC COMPLETE ---\n")

if __name__ == "__main__":
    main()