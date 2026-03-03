import os
import shutil
import subprocess
import pandas as pd

def run_cmd(cmd):
    """Utility to run shell commands with error checking."""
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def parse_fasta(fasta_path):
    """Simple native FASTA parser to remove Biopython dependency."""
    with open(fasta_path, 'r') as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header:
            yield header, "".join(seq_parts)

def make_reads(genome, out_fastq, readlength, depth):
    """Simulates single-end reads by sliding a window across the genome."""
    with open(out_fastq, 'w') as out:
        step = int(readlength / depth)
        recnum = 0
        for header, seq in parse_fasta(genome):
            for i in range(0, len(seq) - readlength, step):
                recnum += 1
                outread = seq[i:i + readlength]
                out.write('@SEQ' + str(recnum) + '\n')
                out.write(str(outread) + '\n+\n')
                out.write('F' * len(outread) + '\n')
    return out_fastq

def wholegenome_task(reference, query, output_prefix, query_id):
    """
    Executes the wholegenome discovery workflow via CLI to bypass 
    Snakemake's Python 3.13 version forcing.
    Simulation -> Alignment -> Find -> Pair -> InferSeq
    """
    # Parameters from author's snippet
    read_length = 150
    read_depth = 10
    min_alignment_quality = 20
    max_direct_repeat_length = 20
    large_insertion_cutoff = 30
    min_alignment_inner_length = max_direct_repeat_length + 1

    # 1. Indexing (Check if index exists to save time)
    if not os.path.exists(reference + ".bwt"):
        run_cmd(f"bwa index {reference}")
    
    if not os.path.exists(query + ".1.bt2"):
        run_cmd(f"bowtie2-build -o 0 -q {query} {query}")

    # 2. Simulate Reads
    query_reads_path = output_prefix + '.query.tmp.fq'
    make_reads(query, query_reads_path, read_length, read_depth)

    # 3. Alignment
    out_sam = output_prefix + '.query.reference.tmp.sam'
    run_cmd(f"bwa mem -t 1 {reference} {query_reads_path} > {out_sam}")

    # 4. Process BAM
    out_bam = output_prefix + '.query.reference.tmp.bam'
    run_cmd(f"samtools sort -o {out_bam} {out_sam}")
    run_cmd(f"samtools index {out_bam}")

    # 5. MGEfinder discovery steps via CLI mappings
    find_file = output_prefix + '.find.tsv'
    run_cmd(f"mgefinder find {out_bam} "
            f"--min_softclip_length 8 "
            f"--min_softclip_count 1 "
            f"--min_alignment_quality {min_alignment_quality} "
            f"--min_alignment_inner_length {min_alignment_inner_length} "
            f"--min_distance_to_mate {max_direct_repeat_length + 2} "
            f"--min_softclip_ratio 0.01 "
            f"--max_indel_ratio 0.0 "
            f"--large_insertion_cutoff {large_insertion_cutoff} "
            f"--min_count_consensus 1 "
            f"--sample_id {query_id} "
            f"--output_file {find_file}")

    pair_file = output_prefix + '.pair.tsv'
    run_cmd(f"mgefinder pair {find_file} {out_bam} {reference} "
            f"--max_direct_repeat_length {max_direct_repeat_length} "
            f"--min_alignment_quality {min_alignment_quality} "
            f"--min_alignment_inner_length {min_alignment_inner_length} "
            f"--max_junction_spanning_prop 0.01 "
            f"--large_insertion_cutoff {large_insertion_cutoff} "
            f"--output_file {pair_file}")

    inferseq_file = output_prefix + '.inferseq.tsv'
    run_cmd(f"mgefinder inferseq-assembly {pair_file} {out_bam} {query} {reference} "
            f"--min_perc_identity 0.95 "
            f"--max_internal_softclip_prop 0.01 "
            f"--max_inferseq_size 500000 "
            f"--min_inferseq_size 30 "
            f"--no-keep-intermediate "
            f"--output_file {inferseq_file}")

    # Cleanup intermediate large files
    for f in [query_reads_path, out_sam, out_bam, out_bam + '.bai']:
        if os.path.exists(f): 
            os.remove(f)
    
    return inferseq_file

def main():
    lsr_genome = snakemake.input.lsr_genome
    related_file = snakemake.input.related_genomes
    out_boundaries = snakemake.output.boundaries
    out_sequences = snakemake.output.sequences
    outdir = os.path.dirname(out_boundaries)
    
    # Parse genome ID (wildcard)
    query_id = snakemake.wildcards.genome

    # Load related references
    with open(related_file, 'r') as f:
        related_genomes = [line.strip() for line in f.readlines() if line.strip()]

    # To satisfy Snakemake, we MUST create the output files even if they are empty
    final_results = []

    if related_genomes:
        for i, ref_genome in enumerate(related_genomes):
            prefix = os.path.join(outdir, f"pair_{i}")
            try:
                # Execute the wholegenome logic via CLI
                res_file = wholegenome_task(ref_genome, lsr_genome, prefix, query_id)
                
                if os.path.exists(res_file):
                    df = pd.read_csv(res_file, sep='\t')
                    if not df.empty:
                        final_results.append(df)
            except Exception as e:
                print(f"Error processing pair {i} with reference {ref_genome}: {e}")
                continue

    # Aggregate and write results
    if final_results:
        combined_df = pd.concat(final_results, ignore_index=True)
        # Write boundaries TSV
        combined_df.to_csv(out_boundaries, sep='\t', index=False)
        
        # Write sequences FASTA
        with open(out_sequences, 'w') as f:
            for _, row in combined_df.iterrows():
                if pd.notna(row['inferred_seq']):
                    f.write(f">{row['sample']}_{row['pair_id']}\n{row['inferred_seq']}\n")
    else:
        # Create empty files with headers if no candidates found
        pd.DataFrame(columns=['sample', 'pair_id', 'method', 'loc', 'inferred_seq_length', 'inferred_seq']).to_csv(out_boundaries, sep='\t', index=False)
        open(out_sequences, 'w').close()

if __name__ == "__main__":
    main()