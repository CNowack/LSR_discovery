import os
import pandas as pd

def load_fasta(fasta_path):
    """Simple FASTA parser to load the genome into a dictionary."""
    seqs = {}
    with open(fasta_path, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                # Store only the primary ID (everything before the first space)
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header:
            seqs[header] = "".join(seq)
    return seqs

def main():
    # 1. Access inputs and outputs from Snakemake
    boundary_file = snakemake.input[0]  # mge_boundaries.tsv
    genome_file = snakemake.input[1]    # .fna assembly
    out_file = snakemake.output.att_sites
    
    # Define how many base pairs to extract on either side of the insertion
    flank_size = 50 

    # 2. FAILSAFE: Initialize the output file immediately
    # This prevents Snakemake MissingOutputExceptions if no MGEs were found
    with open(out_file, 'w') as f:
        f.write("sample\tpair_id\tcontig\tattL_start\tattL_end\tattL_seq\tattR_start\tattR_end\tattR_seq\n")

    # 3. Check if boundary file exists and contains data
    if not os.path.exists(boundary_file):
        return

    try:
        df = pd.read_csv(boundary_file, sep='\t')
        if df.empty or 'loc' not in df.columns:
            return
    except Exception:
        return

    # 4. Load the full genome sequence into memory
    genome_seqs = load_fasta(genome_file)

    results = []
    
    # 5. Parse coordinates and extract flanking DNA
    for _, row in df.iterrows():
        loc = str(row['loc'])
        
        # Ensure the string matches the expected CONTIG:START-END format
        if ":" not in loc or "-" not in loc:
            continue

        # Split the string to isolate the variables
        contig = loc.split(":")[0]
        coords = loc.split(":")[1]
        
        try:
            start = int(coords.split("-")[0])
            end = int(coords.split("-")[1])
        except ValueError:
            continue

        # Handle potential mismatch in NCBI FASTA headers vs MGEfinder contig tracking
        if contig not in genome_seqs:
            matched_contig = next((k for k in genome_seqs.keys() if k.startswith(contig)), None)
            if not matched_contig:
                continue
            contig = matched_contig

        seq = genome_seqs[contig]

        # Extract attL (DNA immediately upstream of the START coordinate)
        attL_start = max(0, start - flank_size)
        attL_seq = seq[attL_start:start]

        # Extract attR (DNA immediately downstream of the END coordinate)
        attR_end = min(len(seq), end + flank_size)
        attR_seq = seq[end:attR_end]

        results.append({
            'sample': row['sample'],
            'pair_id': row['pair_id'],
            'contig': contig,
            'attL_start': attL_start,
            'attL_end': start,
            'attL_seq': attL_seq,
            'attR_start': end,
            'attR_end': attR_end,
            'attR_seq': attR_seq
        })

    # 6. Write the successfully extracted attachment sites to the final TSV
    if results:
        res_df = pd.DataFrame(results)
        res_df.to_csv(out_file, sep='\t', index=False)

if __name__ == "__main__":
    main()