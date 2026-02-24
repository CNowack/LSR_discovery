configfile: "config.yaml"

import pandas as pd

# ── Test mode: load genomes from CSV instead of querying RefSeq ──────────────

if config.get("test_mode"):
    _test_df = pd.read_csv(config["test_genome_csv"])
    GENOMES = _test_df["label"].tolist()  # MUST BE LABELS
    LABEL_TO_ACC = dict(zip(_test_df["label"], _test_df["accession"]))
    BATCHES = ["test_batch"]
else:
    # In production, these will be handled by a checkpoint
    GENOMES  = []   
    BATCHES  = []

FINAL_OUTPUTS = [
    "results/lsr_database.tsv",
    "results/att_sites.fa",
    "results/lsr_clusters_50pct.tsv",
    "results/specificity_predictions.tsv",
    "results/genome_targeting_candidates.tsv"
]

if config.get("test_mode"):
    FINAL_OUTPUTS.append("results/test_validation_report.txt")


rule all:
    input: FINAL_OUTPUTS


# ── Step 1: Download and group genomes by species ──────────────────────────

rule download_assembly_summary:
    output: "data/assembly_summary.txt"
    shell:
        "wget -O {output} {config[refseq_assembly_summary]}"

rule group_genomes_by_species:
    input: "data/assembly_summary.txt"
    output:
        batches_dir=directory("data/batches"),
        species_map="data/species_batches.tsv"
    params:
        batch_size=config["batch_size"]
    script:
        "scripts/group_by_species.py"

rule create_related_list:
    input:
        # This ensures the genomes are actually there before we list them
        done = "data/genomes/{batch_id}/download.done"
    output:
        "data/related/{batch_id}/{genome}_related.txt"
    run:
        import os
        # Look at the actual files on disk in the batch folder
        batch_dir = os.path.dirname(input.done)
        all_fna = [os.path.join(batch_dir, f) for f in os.listdir(batch_dir) if f.endswith('.fna')]
        
        # Exclude the current genome (the 'query') from the 'relatives' list
        current_fna = os.path.join(batch_dir, f"{wildcards.genome}.fna")
        related = [os.path.abspath(f) for f in all_fna if os.path.abspath(f) != os.path.abspath(current_fna)]
        
        with open(output[0], "w") as f:
            f.write("\n".join(related))

rule download_genome_batch:
    input: 
        # If test_mode, we don't need the species map from RefSeq
        species_map = config["test_genome_csv"] if config.get("test_mode") else "data/species_batches.tsv"
    output: 
        filenames = expand("data/genomes/{{batch_id}}/{genome}.fna", genome=GENOMES),
        done = "data/genomes/{batch_id}/download.done"
    params:
        outdir = "data/genomes/{batch_id}",
        acc_mapping = lambda w: ",".join([f"{k}:{v}" for k, v in LABEL_TO_ACC.items()])
    shell:
        """
        mkdir -p {params.outdir}
        python scripts/download_genomes.py \
            --batch {wildcards.batch_id} \
            --mapping "{params.acc_mapping}" \
            --outdir {params.outdir}
        touch {output.done}
        """


# ── Step 2: Annotate CDSs with Prodigal ───────────────────────────────────

rule prodigal:
    input: "data/genomes/{batch_id}/{genome}.fna"
    output:
        proteins="data/annotations/{batch_id}/{genome}.faa",
        gff="data/annotations/{batch_id}/{genome}.gff"
    conda: "envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input} \
            -a {output.proteins} \
            -f gff -o {output.gff} \
            -p meta -q
        """


# ── Step 3: Identify recombinase domain proteins with HMMER ───────────────

rule hmmer_search:
    input:
        proteins="data/annotations/{batch_id}/{genome}.faa",
        hmm=config["recombinase_hmm"]
    output: "data/hmmer/{batch_id}/{genome}.tblout"
    conda: "envs/hmmer.yaml"
    threads: config["threads"]
    shell:
        """
        hmmsearch --tblout {output} \
            --cpu {threads} \
            -E {config[hmmer_evalue]} \
            {input.hmm} {input.proteins} > /dev/null
        """

rule collect_lsr_genomes:
    """Identify genomes that contain at least one LSR candidate."""
    input:
        expand("data/hmmer/{batch_id}/{genome}.tblout",
               batch_id=BATCHES, genome=GENOMES)
    output: "data/lsr_genome_list.tsv"
    script: "scripts/collect_lsr_genomes.py"


# ── Step 4: MGEfinder - compare genomes to find insertion boundaries ───────

rule compute_ani:
    """Use FastANI to find related genomes at >=95% ANI (same species)."""
    input:
        query="data/genomes/{batch_id}/{genome}.fna",
        ref_list="data/species_genome_lists/{species}.txt"
    output: "data/ani/{batch_id}/{genome}.ani"
    conda: "envs/fastani.yaml"
    shell:
        """
        fastANI -q {input.query} \
            --rl {input.ref_list} \
            -o {output} \
            --minFraction 0.5
        """

rule mgefinder_wholegenome:
    """
    Compare LSR-containing genome to related genomes lacking the LSR.
    Uses MGEfinder 'wholegenome' as described in the paper.
    """
    input:
        lsr_genome="data/genomes/{batch_id}/{genome}.fna",
        related_genomes="data/related/{batch_id}/{genome}_related.txt",
        done="data/genomes/{batch_id}/download.done"
    output:
        boundaries="data/mgefinder/{batch_id}/{genome}/mge_boundaries.tsv",
        sequences="data/mgefinder/{batch_id}/{genome}/mge_sequences.fa"
    conda: "envs/mgefinder.yaml"
    shell:
        """
        mgefinder wholegenome \
            {input.lsr_genome} \
            --related-genomes {input.related_genomes} \
            --outdir data/mgefinder/{wildcards.batch_id}/{wildcards.genome}
        """


# ── Step 5: Reconstruct attB and attP attachment sites ────────────────────

rule reconstruct_att_sites:
    """
    From MGE boundaries (attL, attR), reconstruct attB and attP.
    Logic from paper:
        attB = flanking_left + duplication + flanking_right
        attP = terminal_right + duplication + terminal_left
    Extract 50bp window around attachment site center.
    """
    input:
        boundaries="data/mgefinder/{batch_id}/{genome}/mge_boundaries.tsv",
        genome="data/genomes/{batch_id}/{genome}.fna"
    output:
        att_sites="data/att_sites/{batch_id}/{genome}.tsv"
    script: "scripts/reconstruct_att_sites.py"


# ── Step 6: Quality control filtering ────────────────────────────────────

rule quality_filter_lsrs:
    """
    Apply all QC filters from the paper's Methods:
    - Same-species genomes at >=95% ANI
    - Att center length <= 20bp
    - <5% ambiguous nucleotides in att sites
    - LSR length 400-650 aa
    - LSR has Resolvase, Recombinase, or Zn_ribbon_recom Pfam domain
    - <5% ambiguous amino acids
    - MGE size < 200kb
    - LSR within 500nt of predicted att site
    """
    input:
        att_sites=expand("data/att_sites/{batch_id}/{genome}.tsv",
                         batch_id=BATCHES, genome=GENOMES),
        hmmer_results=expand("data/hmmer/{batch_id}/{genome}.tblout",
                             batch_id=BATCHES, genome=GENOMES)
    output:
        filtered="data/lsr_candidates_filtered.tsv",
        fasta="data/lsr_candidates.faa"
    params:
        min_len=config["lsr_min_length_aa"],
        max_len=config["lsr_max_length_aa"],
        max_att_center=config["att_center_max_bp"],
        max_mge_kb=config["max_mge_size_kb"],
        max_dist=config["max_dist_lsr_to_att_bp"]
    script: "scripts/quality_filter_lsrs.py"


# ── Step 7: Cluster LSRs ──────────────────────────────────────────────────

rule cluster_lsrs_90pct:
    input: "data/lsr_candidates.faa"
    output:
        clusters="results/lsr_clusters_90pct.tsv",
        rep_seqs="results/lsr_reps_90pct.faa"
    conda: "envs/mmseqs2.yaml"
    threads: config["threads"]
    shell:
        """
        mmseqs easy-cluster {input} \
            results/lsr_90pct \
            tmp_90pct \
            --min-seq-id {config[cluster_identity_90]} \
            --threads {threads} \
            -c 0.8
        mv results/lsr_90pct_cluster.tsv {output.clusters}
        mv results/lsr_90pct_rep_seq.fasta {output.rep_seqs}
        """

rule cluster_lsrs_50pct:
    input: "data/lsr_candidates.faa"
    output:
        clusters="results/lsr_clusters_50pct.tsv",
        rep_seqs="results/lsr_reps_50pct.faa"
    conda: "envs/mmseqs2.yaml"
    threads: config["threads"]
    shell:
        """
        mmseqs easy-cluster {input} \
            results/lsr_50pct \
            tmp_50pct \
            --min-seq-id {config[cluster_identity_50]} \
            --threads {threads} \
            -c 0.8
        mv results/lsr_50pct_cluster.tsv {output.clusters}
        mv results/lsr_50pct_rep_seq.fasta {output.rep_seqs}
        """


# ── Step 8: Phylogenetic tree ─────────────────────────────────────────────

rule align_lsr_reps:
    input: "results/lsr_reps_50pct.faa"
    output: "results/lsr_reps_50pct.aln"
    conda: "envs/mafft.yaml"
    threads: config["threads"]
    shell:
        "mafft --globalpair --maxiterate 1000 --thread {threads} {input} > {output}"

rule build_phylogeny:
    input: "results/lsr_reps_50pct.aln"
    output:
        tree="results/lsr_phylogeny.treefile",
        log="results/lsr_phylogeny.log"
    conda: "envs/iqtree.yaml"
    threads: config["threads"]
    shell:
        """
        iqtree2 -s {input} \
            -m TEST \
            -B 1000 \
            --prefix results/lsr_phylogeny \
            -T {threads}
        """


# ── Step 9: Predict target site specificity ───────────────────────────────

rule cluster_target_genes:
    """Cluster genes disrupted by MGE insertion at 50% identity."""
    input: "data/lsr_candidates_filtered.tsv"
    output: "results/target_gene_clusters.tsv"
    conda: "envs/mmseqs2.yaml"
    threads: config["threads"]
    script: "scripts/cluster_target_genes.py"

rule predict_specificity:
    """
    Classify LSRs as:
    - Site-specific: single target gene cluster (targeting 1-3 unique clusters)
    - Multi-targeting: >3 target gene clusters (often contain DUF4368)
    Uses 90% clusters for multi-targeting, 50% for site-specific analysis.
    """
    input:
        lsr_clusters_50="results/lsr_clusters_50pct.tsv",
        lsr_clusters_90="results/lsr_clusters_90pct.tsv",
        target_clusters="results/target_gene_clusters.tsv",
        filtered_lsrs="data/lsr_candidates_filtered.tsv"
    output:
        specificity="results/specificity_predictions.tsv",
        site_specific="results/site_specific_lsrs.tsv",
        multi_targeting="results/multi_targeting_lsrs.tsv"
    params:
        min_lsr_clusters_per_target=3  # minimum to call a target gene cluster
    script: "scripts/predict_specificity.py"

# ── Step pre-10: Collect all predicted att sites ────────

rule aggregate_att_sites_fasta:
    """Collects all predicted att sites into a single FASTA for BLASTing."""
    input:
        expand("data/att_sites/{batch_id}/{genome}.tsv", 
               batch_id=BATCHES, genome=GENOMES)
    output:
        "results/att_sites.fa"
    script:
        "scripts/aggregate_fasta.py"


# ── Step 10: Identify genome-targeting candidates (BLAST to human) ────────

rule blast_att_to_human:
    """
    BLAST all attB/attP sequences against human genome (GRCh38).
    Rename best-match att site to attA (acceptor), other to attD (donor).
    Target site in human genome = attH.
    """
    input:
        att_sites="results/att_sites.fa",
        human_genome=config["human_genome"]
    output:
        blast_hits="results/human_blast_hits.tsv",
        genome_targeting="results/genome_targeting_candidates.tsv"
    conda: "envs/blast.yaml"
    threads: config["threads"]
    params:
        evalue=config["blast_evalue"]
    shell:
        """
        # Make BLAST db if needed
        makeblastdb -in {input.human_genome} -dbtype nucl \
            -out resources/GRCh38_db

        blastn -query {input.att_sites} \
            -db resources/GRCh38_db \
            -out {output.blast_hits} \
            -evalue {params.evalue} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -num_threads {threads}

        # Filter and annotate attA/attD/attH designations
        python scripts/annotate_genome_targets.py \
            --blast {output.blast_hits} \
            --att-sites {input.att_sites} \
            --out {output.genome_targeting}
        """


# ── Final: Compile database ───────────────────────────────────────────────

rule compile_lsr_database:
    input:
        # These are already aggregated files
        filtered="data/lsr_candidates_filtered.tsv",
        specificity="results/specificity_predictions.tsv",
        clusters_50="results/lsr_clusters_50pct.tsv",
        genome_targeting="results/genome_targeting_candidates.tsv"
    output:
        database="results/lsr_database.tsv",
    script: "scripts/compile_database.py"

rule validate_test_results:
    input:
        database="results/lsr_database.tsv",
        test_csv=config["test_genome_csv"]
    output:
        report="results/test_validation_report.txt"
    run:
        import pandas as pd
        expected = pd.read_csv(input.test_csv)
        found    = pd.read_csv(input.database, sep="\t")

        lines = []
        for _, row in expected.iterrows():
            acc  = row["accession"]
            exp  = row["expected_lsr"]
            hits = found[found["source_genome"] == acc]
            detected = len(hits) > 0
            status = "PASS" if (detected == exp or exp == "Unknown") else "FAIL"
            lines.append(f"{status}\t{acc}\t{row['label']}\texpected={exp}\tdetected={detected}\thits={len(hits)}")

        with open(output.report, "w") as f:
            f.write("\n".join(lines) + "\n")
        
        print(open(output.report).read())