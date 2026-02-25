import os
import subprocess
import pandas as pd
import glob
import shutil

def parse_fasta_for_id(fasta_path, target_id):
    """Extracts the sequence for a specific protein ID from a FASTA file."""
    if not os.path.exists(fasta_path):
        return ""
    with open(fasta_path, 'r') as f:
        capture = False
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if capture:
                    return "".join(seq)
                header_id = line[1:].split()[0]
                if header_id == target_id:
                    capture = True
            elif capture:
                seq.append(line)
        if capture:
            return "".join(seq)
    return ""

def main():
    filtered_lsrs_file = snakemake.input[0]
    out_clusters = snakemake.output[0]
    threads = str(snakemake.threads)

    # 1. Load Filtered LSR Candidates
    try:
        df_lsr = pd.read_csv(filtered_lsrs_file, sep='\t')
    except pd.errors.EmptyDataError:
        df_lsr = pd.DataFrame()

    # Failsafe: If no candidates exist, write an empty output and exit
    if df_lsr.empty:
        with open(out_clusters, 'w') as f:
            f.write("target_cluster\tlsr_id\n")
        return

    tmp_fasta = "tmp_target_genes.faa"
    tmp_out_prefix = "tmp_target_mmseqs"
    
    records_written = 0

    # 2. Identify and extract the target gene for each LSR
    with open(tmp_fasta, 'w') as f_out:
        for _, row in df_lsr.iterrows():
            genome = row.get('genome', '')
            contig = row.get('contig', '')
            lsr_id = row.get('lsr_id', '')
            
            # Dynamically locate the genome's batch directory
            gff_files = glob.glob(f"data/annotations/*/{genome}.gff")
            if not gff_files:
                continue
                
            gff_path = gff_files[0]
            batch_id = os.path.basename(os.path.dirname(gff_path))
            
            att_path = f"data/att_sites/{batch_id}/{genome}.tsv"
            faa_path = f"data/annotations/{batch_id}/{genome}.faa"
            
            if not os.path.exists(att_path):
                continue

            # Read the attachment site data to find the insertion position
            try:
                df_att = pd.read_csv(att_path, sep='\t')
                # Filter to the same contig
                df_att = df_att[df_att['contig'] == contig]
                if df_att.empty:
                    continue
                # Assuming the first matching MGE on this contig is the correct one
                insertion_pos = df_att.iloc[0].get('insertion_pos', df_att.iloc[0].get('start', 0))
            except Exception:
                continue

            # Parse GFF to find the CDS closest to (or overlapping) the insertion position
            target_gene_id = None
            min_dist = float('inf')
            
            with open(gff_path, 'r') as f_gff:
                for line in f_gff:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9 or parts[2] != 'CDS' or parts[0] != contig:
                        continue
                        
                    start, end = int(parts[3]), int(parts[4])
                    
                    # Calculate distance to insertion (0 if overlapping)
                    if start <= insertion_pos <= end:
                        dist = 0
                    else:
                        dist = min(abs(start - insertion_pos), abs(end - insertion_pos))
                        
                    if dist < min_dist:
                        min_dist = dist
                        # Extract the Prodigal protein ID
                        for attr in parts[8].split(';'):
                            if attr.startswith('ID='):
                                target_gene_id = attr.split('=')[1]
                                break

            # If a valid target gene was found within a reasonable distance (e.g., 2000bp)
            if target_gene_id and min_dist <= 2000:
                seq = parse_fasta_for_id(faa_path, target_gene_id)
                if seq:
                    # We name the FASTA header with the LSR_ID. 
                    # This ensures the MMseqs2 cluster output directly maps LSR_IDs to Target_Clusters.
                    f_out.write(f">{lsr_id}\n{seq}\n")
                    records_written += 1

    # 3. Cluster target genes using MMseqs2
    if records_written > 0:
        cmd = [
            "mmseqs", "easy-cluster",
            tmp_fasta,
            tmp_out_prefix,
            "tmp_mmseqs_work",
            "--min-seq-id", "0.50",
            "-c", "0.8",
            "--threads", threads
        ]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
            # Move the cluster TSV to the final output destination
            cluster_tsv = f"{tmp_out_prefix}_cluster.tsv"
            if os.path.exists(cluster_tsv):
                shutil.move(cluster_tsv, out_clusters)
            else:
                # Failsafe if mmseqs ran but produced no output
                with open(out_clusters, 'w') as f:
                    pass
                    
        except subprocess.CalledProcessError:
            # Failsafe if mmseqs crashes
            with open(out_clusters, 'w') as f:
                pass
    else:
        # Failsafe if no valid target genes were found
        with open(out_clusters, 'w') as f:
            pass

    # 4. Cleanup temporary MMseqs2 files
    for tmp_file in glob.glob(f"{tmp_out_prefix}*") + glob.glob("tmp_mmseqs_work*") + [tmp_fasta]:
        if os.path.isdir(tmp_file):
            shutil.rmtree(tmp_file)
        elif os.path.exists(tmp_file):
            os.remove(tmp_file)

if __name__ == "__main__":
    main()