import os

def main():
    input_files = snakemake.input
    output_file = snakemake.output[0]

    valid_genomes = []

    for tblout in input_files:
        has_hit = False
        
        # Check if the file exists and is not completely empty
        if os.path.exists(tblout) and os.path.getsize(tblout) > 0:
            with open(tblout, 'r') as f:
                for line in f:
                    # HMMER tblout files use '#' for headers and footers.
                    # Any line not starting with '#' is a detected hit.
                    if not line.startswith('#') and line.strip():
                        has_hit = True
                        break # We only need 1 hit to keep the genome in the pipeline
                        
        if has_hit:
            # Extract the genome name and batch ID from the file path
            # Example path: data/hmmer/test_batch/Streptomyces_phage_phiC31.tblout
            genome_name = os.path.basename(tblout).replace('.tblout', '')
            batch_id = os.path.basename(os.path.dirname(tblout))
            
            valid_genomes.append(f"{batch_id}\t{genome_name}\n")

    # Write the compiled list of LSR-containing genomes
    with open(output_file, 'w') as f:
        f.write("batch_id\tgenome\n")
        for record in valid_genomes:
            f.write(record)

if __name__ == "__main__":
    main()