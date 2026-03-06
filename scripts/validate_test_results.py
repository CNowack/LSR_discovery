import pandas as pd
import os
import subprocess
import tempfile

# Color scheme for amino acid visualization
AA_COLORS = {
    'A': '#C8C8C8', 'R': '#145AFF', 'N': '#00D100', 'D': '#E60A0A',
    'C': '#E6E600', 'Q': '#00D100', 'E': '#E60A0A', 'G': '#EBEBEB',
    'H': '#8282D2', 'I': '#0F820F', 'L': '#0F820F', 'K': '#145AFF',
    'M': '#E6E600', 'F': '#3232AA', 'P': '#DC9682', 'S': '#FA9600',
    'T': '#FA9600', 'W': '#B45AB4', 'Y': '#3232AA', 'V': '#0F820F',
    '-': '#FFFFFF'
}

def load_fasta(fasta_path):
    seqs = {}
    if not os.path.exists(fasta_path):
        return seqs
    with open(fasta_path, 'r') as f:
        header = None
        seq = []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header: 
                    # Strip whitespace and trailing stop codons (*) from Prodigal outputs
                    seqs[header] = "".join(seq).upper().strip().rstrip('*')
                header = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if header: 
            seqs[header] = "".join(seq).upper().strip().rstrip('*')
    return seqs

def generate_alignment_html(aligned_dict, title, summary=""):
    """Generates a color-coded HTML visualization for an alignment, sorted by similarity to control."""
    html = [f"<h3>{title}</h3>"]
    
    if summary:
        html.append(f"<p style='font-family: sans-serif; font-size: 15px; margin-bottom: 8px; color: #333;'>{summary}</p>")
        
    html.append("<div style='font-family: monospace; white-space: pre; overflow-x: auto; padding: 10px; background: #f4f4f4; border: 1px solid #ddd; margin-bottom: 15px;'>")
    
    # Get max label length for padding to align the sequences nicely
    if not aligned_dict: return ""
    max_label = max(len(k) for k in aligned_dict.keys())
    
    # Separate Controls from other sequences
    control_keys = sorted([k for k in aligned_dict.keys() if "CONTROL" in k.upper()])
    other_keys = [k for k in aligned_dict.keys() if "CONTROL" not in k.upper()]
    
    # Sort subsequent sequences by similarity to the primary control
    if control_keys:
        ref_seq = aligned_dict[control_keys[0]]
        
        def get_similarity(k):
            target_seq = aligned_dict[k]
            # Count exact amino acid matches (ignoring gaps in the reference)
            return sum(1 for c, t in zip(ref_seq, target_seq) if c == t and c != '-')
            
        other_keys.sort(key=get_similarity, reverse=True)
    else:
        other_keys.sort()
        
    keys = control_keys + other_keys
    
    for label in keys:
        sequence = aligned_dict[label]
        row = [f"<span style='display:inline-block; width:{max_label + 2}ch; font-weight:bold;'>{label}</span>"]
        for char in sequence:
            color = AA_COLORS.get(char, "#FFFFFF")
            # Dim the color slightly if it's a gap for better visual separation
            text_color = "black" if char != '-' else "#ccc"
            row.append(f"<span style='background-color:{color}; color: {text_color}; padding: 0 1px;'>{char}</span>")
        html.append("".join(row) + "<br>")
    
    html.append("</div>")
    return "\n".join(html)

def run_mafft(sequences_dict):
    """Runs MAFFT on a dictionary of sequences and returns aligned sequences."""
    if len(sequences_dict) < 2:
        return sequences_dict

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_in:
        for label, seq in sequences_dict.items():
            tmp_in.write(f">{label}\n{seq}\n")
        tmp_in_path = tmp_in.name

    try:
        # Run MAFFT (assuming it's in the PATH via conda)
        result = subprocess.run(
            ['mafft', '--auto', '--quiet', tmp_in_path],
            capture_output=True, text=True, check=True
        )
        aligned_fasta = result.stdout
        
        # Parse output
        aligned_dict = {}
        current_id = None
        for line in aligned_fasta.splitlines():
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                aligned_dict[current_id] = ""
            elif current_id:
                aligned_dict[current_id] += line.strip()
        return aligned_dict
    except Exception as e:
        print(f"Warning: MSA failed: {e}")
        return sequences_dict
    finally:
        if os.path.exists(tmp_in_path):
            os.remove(tmp_in_path)

def main():
    metadata_path = snakemake.input.metadata 
    fasta_path = snakemake.input.fasta
    test_csv_path = snakemake.config.get("test_genome_csv")
    
    output_log = snakemake.output.log if hasattr(snakemake.output, 'log') else snakemake.output[0]
    output_html = snakemake.output.html if hasattr(snakemake.output, 'html') else snakemake.output[1]

    df_meta = pd.read_csv(metadata_path, sep='\t') if os.path.exists(metadata_path) else pd.DataFrame()
    fasta_seqs = load_fasta(fasta_path)
    df_test = pd.read_csv(test_csv_path)

    ACC_COL, LABEL_COL, EXPECTED_COL, SEQ_COL = 'accession', 'label', 'expected_lsr', 'positive_seq'

    html_reports = [
        "<html><head><title>LSR Validation Alignments</title></head>",
        "<body style='font-family: sans-serif; padding: 20px;'>",
        "<h1>LSR Discovery Validation: Alignments Report</h1>"
    ]

    with open(output_log, 'w') as out_file:
        def log_and_print(msg):
            print(msg); out_file.write(msg + "\n")

        log_and_print("--- LSR DISCOVERY VALIDATION REPORT ---")

        # Iterate by Species (label)
        for species_label, species_df in df_test.groupby(LABEL_COL):
            log_and_print(f"\n=== SPECIES: {species_label} ===")
            html_reports.append(f"<hr><h2>Species: {species_label}</h2>")
            
            species_seqs_to_align = {}

            # Iterate by Genome (accession) within the species
            for _, row in species_df.iterrows():
                accession = row[ACC_COL]
                expected = str(row[EXPECTED_COL]).lower() == 'true'
                
                # Extract control seq, clean whitespace, and strip terminal stop codons
                control_seq = str(row[SEQ_COL]).strip().replace('\n', '').replace('\r', '').upper().rstrip('*') if not pd.isna(row[SEQ_COL]) else None
                
                genome_id = f"{species_label}_{accession}"
                genome_hits = df_meta[df_meta['source_genome'] == genome_id] if not df_meta.empty and genome_id in df_meta['source_genome'].values else pd.DataFrame()
                
                detected = len(genome_hits) > 0
                status = "PASS" if expected == detected else "FAIL"
                log_and_print(f"{status}\t{accession}\texpected={expected}\tdetected={detected}\thits={len(genome_hits)}")

                genome_seqs_to_align = {}

                # Include control sequence if it exists
                if control_seq:
                    control_key = f"CONTROL_{accession}"
                    genome_seqs_to_align[control_key] = control_seq
                    species_seqs_to_align[control_key] = control_seq

                if detected:
                    html_reports.append(f"<h3>Genome: {genome_id}</h3>")
                    
                    for _, lsr_row in genome_hits.iterrows():
                        header = f"{lsr_row['lsr_id']}_{lsr_row['pair_id']}"
                        raw_discovered = fasta_seqs.get(header, "")
                        
                        genome_seqs_to_align[header] = raw_discovered
                        species_seqs_to_align[header] = raw_discovered

                        # Only do pairwise if there's a control sequence to compare against
                        if control_seq:
                            to_align_pair = {control_key: control_seq, header: raw_discovered}
                            aligned_pair = run_mafft(to_align_pair)
                            
                            # Calculate Pairwise Homology and SNPs
                            aligned_control = aligned_pair.get(control_key, "")
                            aligned_target = aligned_pair.get(header, "")
                            
                            matches = 0
                            snps = []
                            control_pos = 0
                            
                            for i in range(len(aligned_control)):
                                c_char = aligned_control[i]
                                t_char = aligned_target[i] if i < len(aligned_target) else '-'
                                
                                if c_char != '-':
                                    control_pos += 1
                                    
                                if c_char == t_char and c_char != '-':
                                    matches += 1
                                elif c_char != '-' and t_char != '-':
                                    # It's a SNP (Mismatch without gaps)
                                    snps.append(f"{c_char}{control_pos}{t_char}")
                                    
                            control_len = len([c for c in aligned_control if c != '-'])
                            homology = (matches / control_len) * 100 if control_len > 0 else 0
                            
                            if not snps:
                                snp_str = "None"
                            elif len(snps) > 20:
                                snp_str = ", ".join(snps[:20]) + f" ... (+{len(snps)-20} more)"
                            else:
                                snp_str = ", ".join(snps)
                                
                            summary_text = f"<b>Homology:</b> {homology:.1f}% | <b>SNPs relative to control:</b> {snp_str}"
                            
                            if raw_discovered == control_seq:
                                log_and_print(f"\t- {header}: EXACT MATCH")
                            else:
                                log_and_print(f"\t- {header}: Mismatch (view HTML)")

                            html_reports.append(generate_alignment_html(aligned_pair, f"Pairwise: Control vs {header}", summary=summary_text))
                        else:
                            log_and_print(f"\t- {header}: Discovered (No control to match)")

                    # Genome MSA
                    if len(genome_seqs_to_align) > 2 or (not control_seq and len(genome_seqs_to_align) > 1):
                        aligned_genome_msa = run_mafft(genome_seqs_to_align)
                        html_reports.append(generate_alignment_html(aligned_genome_msa, f"Genome MSA: All LSRs in {genome_id}"))
                        
                elif control_seq and not detected:
                    log_and_print("\t- No LSRs discovered to compare.")
            
            # Species MSA
            if len(species_seqs_to_align) > 1:
                aligned_species_msa = run_mafft(species_seqs_to_align)
                html_reports.append(generate_alignment_html(aligned_species_msa, f"Species MSA: All Discoveries in {species_label}"))
                
        log_and_print("\n---------------------------------------")

    html_reports.append("</body></html>")
    with open(output_html, 'w') as f:
        f.write("\n".join(html_reports))

if __name__ == "__main__":
    main()