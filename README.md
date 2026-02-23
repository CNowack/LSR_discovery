 # LSR Discovery Pipeline

A reimplementation of the computational workflow described in:

> Durrant, M.G., et al. "Systematic discovery of recombinases for efficient integration of large DNA sequences into the human genome." *Nature Biotechnology* (2023). https://doi.org/10.1038/s41587-022-01494-w

see (https://github.com/bhattlab/GenomeSearch)

The original study mined 194,585 bacterial genomes to identify 6,207 large serine recombinases (LSRs) and their cognate attachment sites, expanding known LSR diversity more than 100-fold. This pipeline reproduces that discovery workflow using Snakemake.

---
Added test_mode to specify input genomes to validate pipeline.

---

## Scripts

**`scripts/group_by_species.py`**
Reads the RefSeq assembly summary and groups genomes by species using NCBI taxonomy. Outputs a TSV mapping each species to its list of genome accessions, and splits those into batches (50 for first batch, 20 thereafter). This drives all downstream wildcard expansion.

**`scripts/download_genomes.py`**
Takes a batch ID and the species map, fetches the corresponding genome FASTA files from NCBI FTP, and saves them to the appropriate batch directory. Also handles retries and skips already-downloaded files.

**`scripts/collect_lsr_genomes.py`**
Scans all HMMER output files and identifies which genomes contain at least one protein hit against PF07508. Produces a list of LSR-positive genomes to feed into MGEfinder — you only want to run the expensive comparative genomics step on genomes that actually have a candidate.

**`scripts/reconstruct_att_sites.py`**
The core biological logic. Takes MGEfinder boundary output and reconstructs attB and attP from the post-integration structure (B1-D-P1-MGE-P2-D-B2). Extracts the 50bp window around the dinucleotide core center for each candidate.

**`scripts/quality_filter_lsrs.py`**
Applies all eight QC filters from the paper's Methods section — LSR length, ambiguous residues, attachment site center length, MGE size, distance from LSR to att site, etc. This is where most candidates get dropped; the paper went from ~12,000 initial candidates to 6,207 after this step.

**`scripts/cluster_target_genes.py`**
Takes the genes disrupted by MGE insertion (identified via MGEfinder boundaries) and clusters them at 50% identity using MMseqs2. The number of distinct target gene clusters an LSR hits is what determines whether it's site-specific or multi-targeting.

**`scripts/predict_specificity.py`**
Uses the target gene cluster assignments to classify each LSR. If it hits 1–3 unique gene clusters it's site-specific; if it hits more than 3 diverse clusters it's multi-targeting. Also checks for DUF4368 domain presence as a corroborating signal for multi-targeting LSRs.

**`scripts/annotate_genome_targets.py`**
Post-processes the BLAST results against GRCh38. Assigns the attA/attD/attH nomenclature — the att site with the best human genome match becomes attA (acceptor), the other becomes attD (donor), and the human locus is attH. Filters by E-value and outputs the genome-targeting candidate table.

**`scripts/compile_database.py`**
Joins all intermediate outputs — filtered LSRs, att sites, specificity predictions, cluster assignments, and genome-targeting hits — into the final `lsr_database.tsv`, equivalent to Supplementary Table 1 in the paper.
 