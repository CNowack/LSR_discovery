import os
import pandas as pd

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def load_fasta(fasta_path):
    seqs = {}
    if not os.path.exists(fasta_path):
        return seqs
    with open(fasta_path, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            seqs[header] = "".join(seq)
    return seqs

def get_core_homology(attL, attR, min_bp, max_bp):
    if attL == "NOT_FOUND" or attR == "NOT_FOUND":
        return None
    end_L = str(attL)[-10:]
    start_R = str(attR)[:10]
    for length in range(max_bp, min_bp - 1, -1):
        for i in range(len(end_L) - length + 1):
            kmer = end_L[i:i+length]
            if kmer in start_R:
                return kmer
    return None

def parse_prodigal_faa(faa_path):
    """Parses both the coordinates from the header AND the amino acid sequence."""
    data = {}
    if not os.path.exists(faa_path):
        return data
    with open(faa_path, 'r') as f:
        prot_id = None
        coords = {}
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save the previous record before starting a new one
                if prot_id:
                    data[prot_id] = {'coords': coords, 'seq': "".join(seq)}
                
                parts = line.split(' # ')
                if len(parts) >= 4:
                    prot_id = parts[0][1:].split()[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    strand = int(parts[3])
                    contig = prot_id.rsplit('_', 1)[0]
                    coords = {'contig': contig, 'start': start, 'end': end, 'strand': strand}
                else:
                    prot_id = None
                seq = []
            else:
                seq.append(line)
        # Catch the final record in the file
        if prot_id:
            data[prot_id] = {'coords': coords, 'seq': "".join(seq)}
    return data

def main():
    mge_files = snakemake.input.mge_data
    hmmer_files = snakemake.input.hmmer_data
    att_files = snakemake.input.att_data
    
    out_tsv = snakemake.output.filtered
    out_fasta = snakemake.output.fasta

    MAX_MGE_SIZE = snakemake.params.get("max_mge_size_kb", 500) * 1000
    MIN_MGE_SIZE = snakemake.params.get("min_mge_size_bp", 30)
    MAX_LSR_DIST = snakemake.params.get("max_lsr_distance_bp", 10000)
    MIN_CORE = snakemake.params.get("min_core_bp", 2)
    MAX_CORE = snakemake.params.get("max_core_bp", 4)

    filtered_candidates = []

    genome_map = {}
    for f in mge_files:
        parts = f.split('/')
        batch = parts[-3]
        gid = parts[-2]
        genome_map.setdefault(gid, {})['mge'] = f
        genome_map[gid]['batch'] = batch
    
    for f in hmmer_files:
        gid = os.path.basename(f).replace('.tblout', '')
        genome_map.setdefault(gid, {})['hmmer'] = f
        genome_map[gid]['faa'] = f.replace('hmmer', 'annotations').replace('.tblout', '.faa')

    for f in att_files:
        gid = os.path.basename(f).replace('.tsv', '')
        genome_map.setdefault(gid, {})['att'] = f

    for genome_id, files in genome_map.items():
        mge_path = files.get('mge')
        hmmer_path = files.get('hmmer')
        att_path = files.get('att')
        faa_path = files.get('faa')
        batch = files.get('batch')
        
        if not all([mge_path, hmmer_path, att_path]) or not all(os.path.exists(p) for p in [mge_path, hmmer_path, att_path]):
            continue
            
        try:
            df_mge = pd.read_csv(mge_path, sep='\t')
            df_att = pd.read_csv(att_path, sep='\t')
            df_hmmer = pd.read_csv(hmmer_path, sep=r'\s+', comment='#', header=None)
        except Exception:
            continue

        if df_hmmer.empty or df_mge.empty or df_att.empty:
            continue

        if 'method' in df_mge.columns:
            df_mge = df_mge.sort_values(by='method')
        df_mge = df_mge.drop_duplicates(subset=['pair_id'], keep='first')

        # Extract both coordinates AND amino acid sequences
        prot_data = parse_prodigal_faa(faa_path)
        
        genome_fasta_path = f"data/genomes/{batch}/{genome_id}.fna"
        genome_seqs = load_fasta(genome_fasta_path)

        for _, mge_row in df_mge.iterrows():
            mge_size = mge_row.get('inferred_seq_length', 0)
            
            if pd.isna(mge_size) or not (MIN_MGE_SIZE <= mge_size <= MAX_MGE_SIZE):
                continue
            
            pair_id = mge_row.get('pair_id')
            mge_loc = str(mge_row.get('loc'))
            
            try:
                mge_contig = mge_loc.split(":")[0]
                mge_start = int(mge_loc.split(":")[1].split("-")[0])
                mge_end = int(mge_loc.split(":")[1].split("-")[1])
            except (IndexError, ValueError):
                continue
            
            att_match = df_att[df_att['pair_id'] == pair_id]
            attL_seq = att_match['attL_seq'].values[0] if not att_match.empty else "NOT_FOUND"
            attR_seq = att_match['attR_seq'].values[0] if not att_match.empty else "NOT_FOUND"

            core_seq = get_core_homology(attL_seq, attR_seq, MIN_CORE, MAX_CORE)
            if not core_seq:
                continue 

            valid_lsr_found = False
            associated_lsr = None

            for _, hmm_hit in df_hmmer.iterrows():
                protein_id = hmm_hit[0]
                if protein_id not in prot_data:
                    continue
                
                p_data = prot_data[protein_id]['coords']
                if p_data['contig'] != mge_contig:
                    continue
                
                if p_data['end'] >= mge_start and p_data['start'] <= mge_end:
                    dist = 0
                else:
                    dist = min(abs(mge_start - p_data['end']), abs(p_data['start'] - mge_end))
                
                if dist <= MAX_LSR_DIST:
                    valid_lsr_found = True
                    associated_lsr = protein_id
                    break 

            if valid_lsr_found:
                # 1. Grab Amino Acid sequence (from the .faa file parsing)
                lsr_aa_seq = prot_data[associated_lsr]['seq']

                # 2. Grab DNA sequence (from the whole genome .fna slicing)
                lsr_dna_seq = "NOT_FOUND"
                p_data = prot_data[associated_lsr]['coords']
                
                matched_contig = p_data['contig']
                if matched_contig not in genome_seqs:
                    matched_contig = next((k for k in genome_seqs.keys() if k.startswith(p_data['contig'])), None)
                    
                if matched_contig and matched_contig in genome_seqs:
                    seq_slice = genome_seqs[matched_contig][p_data['start']-1 : p_data['end']]
                    if p_data['strand'] == -1:
                        seq_slice = reverse_complement(seq_slice)
                    lsr_dna_seq = seq_slice

                filtered_candidates.append({
                    'source_genome': genome_id,
                    'lsr_id': associated_lsr,
                    'pair_id': pair_id,
                    'mge_size': mge_size,
                    'core_homology': core_seq,
                    'lsr_dna_seq': lsr_dna_seq,
                    'lsr_aa_seq': lsr_aa_seq, # Added amino acid sequence to dictionary
                    'attL_seq': attL_seq,
                    'attR_seq': attR_seq
                })

    df_out = pd.DataFrame(filtered_candidates)
    
    if not df_out.empty:
        df_out = df_out.drop_duplicates()
        
        # Save the TSV
        tsv_cols = [c for c in df_out.columns if c != 'lsr_aa_seq']
        df_out[tsv_cols].to_csv(out_tsv, sep='\t', index=False)
        
        # Save the FASTA
        with open(out_fasta, 'w') as f:
            for _, row in df_out.iterrows():
                # Write the actual amino acid sequence here
                f.write(f">{row['lsr_id']}_{row['pair_id']}\n{row['lsr_aa_seq']}\n")
    else:
        pd.DataFrame(columns=['source_genome', 'lsr_id', 'pair_id', 'mge_size', 'core_homology', 'lsr_dna_seq', 'attL_seq', 'attR_seq']).to_csv(out_tsv, sep='\t', index=False)
        open(out_fasta, 'w').close()

if __name__ == "__main__":
    main()