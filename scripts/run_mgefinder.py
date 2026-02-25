import os
import subprocess
import glob
import shutil

def main():
    lsr_genome = snakemake.input.lsr_genome
    related_file = snakemake.input.related_genomes
    out_boundaries = snakemake.output.boundaries
    out_sequences = snakemake.output.sequences
    outdir = os.path.dirname(out_boundaries)
    
    # Initialize empty outputs to prevent Snakemake MissingOutputExceptions
    # in the event that no MGEs are found in any comparison
    with open(out_boundaries, 'w') as f:
        f.write("genome\tinsertion_id\tlength\tsequence\n") # Dummy header
        
    with open(out_sequences, 'w') as f:
        pass

    # Load related genomes
    with open(related_file, 'r') as f:
        related_genomes = [line.strip() for line in f.readlines() if line.strip()]

    if not related_genomes:
        return

    header_written = False

    # MGEfinder wholegenome performs 1-vs-1 comparisons. 
    # We loop through all related genomes (references lacking the LSR) 
    # and compare them to the query (genome containing the LSR candidate).
    for i, ref_genome in enumerate(related_genomes):
        temp_outdir = os.path.join(outdir, f"pair_{i}")
        
        cmd = [
            "mgefinder", "wholegenome",
            ref_genome, lsr_genome,
            "--outdir", temp_outdir
        ]
        
        try:
            # Execute MGEfinder for this pair
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Search for output TSV and FASTA files dynamically
            found_tsvs = glob.glob(os.path.join(temp_outdir, "**", "*.tsv"), recursive=True)
            found_fas = glob.glob(os.path.join(temp_outdir, "**", "*.fa"), recursive=True) + \
                        glob.glob(os.path.join(temp_outdir, "**", "*.fasta"), recursive=True)
            
            # Aggregate TSV data
            for tsv in found_tsvs:
                with open(tsv, 'r') as infile:
                    lines = infile.readlines()
                    if len(lines) <= 1: 
                        continue # Skip if empty or header-only
                        
                    # Overwrite dummy header on first valid find, append thereafter
                    mode = 'w' if not header_written else 'a'
                    with open(out_boundaries, mode) as outfile:
                        if not header_written:
                            outfile.writelines(lines)
                            header_written = True
                        else:
                            outfile.writelines(lines[1:]) # Skip header for subsequent files
                            
            # Aggregate FASTA data
            for fa in found_fas:
                with open(fa, 'r') as infile:
                    content = infile.read()
                    if content.strip():
                        with open(out_sequences, 'a') as outfile:
                            outfile.write(content)
                            if not content.endswith('\n'):
                                outfile.write('\n')
                                
        except subprocess.CalledProcessError:
            # If a pair fails (e.g., genomes are too divergent for BWA to align),
            # catch the error and continue to the next related genome
            pass
        finally:
            # Clean up the temporary directory to manage SCC storage constraints
            if os.path.exists(temp_outdir):
                shutil.rmtree(temp_outdir)

if __name__ == "__main__":
    main()