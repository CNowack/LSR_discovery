import os
import sys
import pandas as pd

def parse_fasta(fasta_path):
    """Simple FASTA parser yielding (header, sequence)."""
    if not os.path.exists(fasta_path):
        return
    with open(fasta_path, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line[1:].split()[0] # Get ID only
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

def parse_gff(gff_path, valid_ids):
    """Parse GFF to get contig, start, end for specific protein IDs."""
    coords = {}
    if not os.path.exists(gff_path):
        return coords
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            
            # Prodigal outputs ID=... in the attributes column
            for attr in attributes.split(';'):
                if attr.startswith("ID="):
                    prot_id = attr.split("=")[1]
                    # Sometimes Prodigal appends numbers differently in GFF vs FAA, 
                    # but usually they match the FASTA header exactly.
                    if prot_id in valid_ids:
                        coords[prot_id] = {'contig': contig, 'start': start, 'end': end}
    return coords

def parse_hmmer(tblout_path):
    """Parse HMMER tblout to get a list of protein IDs that matched the HMM."""
    hits = set()
    if not os.path.exists(tblout_path):
        return hits
    with open(tblout_path, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.split()
                if len(parts) > 0:
                    hits.add(parts[0]) # Target name (protein ID)
    return hits

def main():
    # 1. Load Snakemake I/O and parameters
    att_site_files = snakemake.input.att_sites
    hmmer_files = snakemake.input.hmmer_results
    
    out_filtered = snakemake.output.filtered
    out_fasta = snakemake.output.fasta
    
    # Paper thresholds
    min_len = int(snakemake.params.get("min_len", 400))
    max_len = int(snakemake.params.get("max_len", 650))
    max_att_center = int(snakemake.params.get("max_att_center", 20))
    max_mge_kb = float(snakemake.params.get("max_mge_kb", 200.0))
    max_dist = int(snakemake.params.get("max_dist", 500))
    
    max_mge_bp = max_mge_kb * 1000

    valid_candidates = []
    fasta_records_to_write = []

    # 2. Process each genome's HMMER results
    hmmer_dict = {}
    for h_file in hmmer_files:
        genome_name = os.path.basename(h_file).replace('.tblout', '')
        hmmer_dict[genome_name] = h_file

    # 3. Iterate through att_site files and correlate with HMMER/Annotations
    for att_file in att_site_files:
        if os.path.getsize(att_file) == 0:
            continue
            
        genome_name = os.path.basename(att_file).replace('.tsv', '')
        
        try:
            df_att = pd.read_csv(att_file, sep='\t')
        except pd.errors.EmptyDataError:
            continue

        if df_att.empty:
            continue

        # Get the corresponding annotation files
        # (Assuming standard folder structures defined in your Snakefile)
        batch_dir = os.path.basename(os.path.dirname(att_file))
        faa_path = os.path.join("data", "annotations", batch_dir, f"{genome_name}.faa")
        gff_path = os.path.join("data", "annotations", batch_dir, f"{genome_name}.gff")
        tblout_path = hmmer_dict.get(genome_name)

        if not tblout_path or not os.path.exists(faa_path) or not os.path.exists(gff_path):
            continue

        # Parse HMMER to get candidate IDs
        hmmer_hits = parse_hmmer(tblout_path)
        if not hmmer_hits:
            continue

        # Filter proteins by sequence criteria
        valid_proteins = {}
        for header, seq in parse_fasta(faa_path):
            if header in hmmer_hits:
                seq_len = len(seq)
                ambig_frac = seq.upper().count('X') / seq_len if seq_len > 0 else 1.0
                
                # Paper Filter: Length 400-650, <5% ambiguous AAs
                if min_len <= seq_len <= max_len and ambig_frac < 0.05:
                    valid_proteins[header] = seq

        if not valid_proteins:
            continue

        # Get genomic coordinates for valid proteins
        prot_coords = parse_gff(gff_path, valid_proteins.keys())

        # 4. Filter MGEs and check Distance
        for _, mge_row in df_att.iterrows():
            mge_size = mge_row.get('mge_size', 0)
            att_center_len = mge_row.get('att_center_len', 0)
            mge_contig = mge_row.get('contig', '')
            insertion_pos = mge_row.get('insertion_pos', mge_row.get('start', 0))
            att_seq = str(mge_row.get('sequence', '')).upper()
            
            # Paper Filter: MGE < 200kb, att center <= 20bp
            if mge_size > max_mge_bp or att_center_len > max_att_center:
                continue
                
            # Paper Filter: < 5% ambiguous nucleotides in att site
            if len(att_seq) > 0 and (att_seq.count('N') / len(att_seq)) >= 0.05:
                continue

            # Check distance for all valid LSRs on the same contig
            for prot_id, coords in prot_coords.items():
                if coords['contig'] == mge_contig:
                    # Calculate distance from LSR boundaries to MGE insertion site
                    dist_to_start = abs(coords['start'] - insertion_pos)
                    dist_to_end = abs(coords['end'] - insertion_pos)
                    min_dist = min(dist_to_start, dist_to_end)
                    
                    # Paper Filter: LSR within 500nt of predicted att site
                    if min_dist <= max_dist:
                        valid_candidates.append({
                            'genome': genome_name,
                            'lsr_id': prot_id,
                            'contig': coords['contig'],
                            'lsr_start': coords['start'],
                            'lsr_end': coords['end'],
                            'distance_to_att': min_dist,
                            'mge_size': mge_size
                        })
                        
                        fasta_records_to_write.append((prot_id, valid_proteins[prot_id]))

    # 5. Deduplicate and Write Outputs
    df_filtered = pd.DataFrame(valid_candidates).drop_duplicates()
    
    # Write TSV Failsafe
    if not df_filtered.empty:
        df_filtered.to_csv(out_filtered, sep='\t', index=False)
    else:
        with open(out_filtered, 'w') as f:
            f.write("genome\tlsr_id\tcontig\tlsr_start\tlsr_end\tdistance_to_att\tmge_size\n")

    # Write FASTA Failsafe
    with open(out_fasta, 'w') as f_out:
        written_ids = set()
        for prot_id, seq in fasta_records_to_write:
            if prot_id not in written_ids: # Avoid writing the same LSR twice if it matched multiple MGE rows
                f_out.write(f">{prot_id}\n{seq}\n")
                written_ids.add(prot_id)

if __name__ == "__main__":
    main()