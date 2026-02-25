import pandas as pd
import argparse
import os

def parse_fasta(fasta_file):
    """Parses a FASTA file into a dictionary of {header: sequence}."""
    sequences = {}
    if not os.path.exists(fasta_file):
        return sequences
        
    with open(fasta_file, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq)
                header = line[1:] # Remove '>'
                seq = []
            else:
                seq.append(line)
        if header:
            sequences[header] = "".join(seq)
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Annotate genome-targeting LSR candidates.")
    parser.add_argument("--blast", required=True, help="BLAST TSV output file")
    parser.add_argument("--att-sites", required=True, help="FASTA file of predicted att sites")
    parser.add_argument("--out", required=True, help="Output TSV file")
    args = parser.parse_args()

    # 1. Parse the FASTA sequences
    fasta_seqs = parse_fasta(args.att_sites)
    
    # 2. Group the att sites to find pairs (attB and attP from the same integration event)
    # The header format defined in aggregate_fasta.py is: {genome_id}_{site_type}_{index}
    grouped_by_base = {}
    for header in fasta_seqs.keys():
        parts = header.rsplit('_', 2)
        if len(parts) == 3:
            genome_id, site_type, idx = parts
            base_key = f"{genome_id}_{idx}"
            if base_key not in grouped_by_base:
                grouped_by_base[base_key] = {}
            grouped_by_base[base_key][site_type] = header

    # Create a mapping of each site to its counterpart (Donor vs Acceptor)
    pairs = {}
    for base_key, sites in grouped_by_base.items():
        site_types = list(sites.keys())
        if len(site_types) == 2:
            t1, t2 = site_types[0], site_types[1]
            pairs[sites[t1]] = sites[t2]
            pairs[sites[t2]] = sites[t1]
        elif len(site_types) == 1:
            pairs[sites[site_types[0]]] = None

    # 3. Read the BLAST results and handle empty outcomes
    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    
    # Failsafe: if BLAST found no hits, write the empty header and exit
    if os.path.getsize(args.blast) == 0:
        with open(args.out, 'w') as f:
            f.write("genome\tinsertion_id\tattA_id\tattA_seq\tattD_id\tattD_seq\thuman_chr\thuman_start\thuman_end\tpident\tevalue\tbitscore\n")
        return

    df_blast = pd.read_csv(args.blast, sep='\t', names=columns)
    
    # Sort hits to ensure the top alignment (lowest e-value, highest bitscore) appears first
    df_blast = df_blast.sort_values(by=["evalue", "bitscore"], ascending=[True, False])

    # 4. Annotate Acceptor (attA) and Donor (attD) designations
    results = []
    for _, row in df_blast.iterrows():
        qseqid = row["qseqid"]
        
        parts = qseqid.rsplit('_', 2)
        if len(parts) == 3:
            genome_id, _, idx = parts
        else:
            genome_id = qseqid
            idx = "0"
            
        attA_id = qseqid
        attA_seq = fasta_seqs.get(qseqid, "")
        
        attD_id = pairs.get(qseqid)
        attD_seq = fasta_seqs.get(attD_id, "") if attD_id else ""
        
        results.append({
            "genome": genome_id,
            "insertion_id": idx,
            "attA_id": attA_id,
            "attA_seq": attA_seq,
            "attD_id": attD_id,
            "attD_seq": attD_seq,
            "human_chr": row["sseqid"],         # This represents attH
            "human_start": row["sstart"],
            "human_end": row["send"],
            "pident": row["pident"],
            "evalue": row["evalue"],
            "bitscore": row["bitscore"]
        })
        
    df_out = pd.DataFrame(results)
    
    # Deduplicate in the event a sequence hits the exact same locus multiple times
    df_out = df_out.drop_duplicates(subset=["attA_id", "human_chr", "human_start", "human_end"])
    df_out.to_csv(args.out, sep='\t', index=False)

if __name__ == "__main__":
    main()