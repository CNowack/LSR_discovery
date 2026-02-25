import pandas as pd
import os

def main():
    # Snakemake automatically passes inputs and outputs to Python scripts
    input_files = snakemake.input
    output_fasta = snakemake.output[0]
    
    with open(output_fasta, 'w') as fasta_out:
        for tsv_file in input_files:
            # Check if file is empty to avoid Pandas errors on genomes with no LSRs
            if os.path.getsize(tsv_file) == 0:
                continue
                
            try:
                df = pd.read_csv(tsv_file, sep='\t')
                
                # Iterate through each predicted attachment site
                for index, row in df.iterrows():
                    # Extract identifying info for the FASTA header
                    genome_id = row.get('genome', 'unknown_genome')
                    site_type = row.get('site_type', 'att') # e.g., attB, attP
                    sequence = row.get('sequence', '')
                    
                    # Create a unique header so we can trace BLAST hits back to the source
                    if pd.notna(sequence) and sequence.strip():
                        header = f">{genome_id}_{site_type}_{index}"
                        fasta_out.write(f"{header}\n{sequence}\n")
                        
            except pd.errors.EmptyDataError:
                pass # Gracefully skip if the TSV only has headers

if __name__ == "__main__":
    main()