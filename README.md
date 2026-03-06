 # LSR Discovery Pipeline

A reimplementation of the computational workflow described in:

> Durrant, M.G., et al. "Systematic discovery of recombinases for efficient integration of large DNA sequences into the human genome." *Nature Biotechnology* (2023). https://doi.org/10.1038/s41587-022-01494-w

see:

    https://github.com/bhattlab/MGEfinder

    https://github.com/bhattlab/GenomeSearch

The original study mined 194,585 bacterial genomes to identify 6,207 large serine recombinases (LSRs) and their cognate attachment sites, expanding known LSR diversity more than 100-fold. This pipeline reproduces that discovery workflow using Snakemake.

## Additions:
* Added test_mode to specify input genomes to validate pipeline.
* Added validation report to quickly check LSRs discovered by the pipeline against those provided in the author's supplementary data


## Overview

This repository contains a fully automated Snakemake pipeline designed to systematically discover, structurally validate, and categorize Large Serine Recombinases (LSRs) from bacterial genomes. The workflow processes whole-genome assemblies to identify novel LSRs for use in genetic engineering.

## Pipeline Architecture

**The computational Directed Acyclic Graph (DAG) executes the following sequential stages:**

1. **Protein Prediction (Prodigal):** Predicts Open Reading Frames (ORFs) from raw bacterial contigs.

2. **LSR Identification (HMMER):** Identifies recombinase candidates using the PF07508 Pfam profile.

3. **Structural Discovery (MGEfinder):** Aligns occupied vs. empty isolate genomes to find mobile genetic element (MGE) boundaries.

4. **Site Reconstruction:** Mathematically recovers pre-integration attB and attP attachment sites.

5. **Quality Filtering:** Enforces thresholds for protein length, MGE size, and distance to the insertion site. Extracts complete DNA and AA sequences.

6. **Family Clustering (MMseqs2):** Defines LSR families at 50% and 90% identity thresholds.

7. **Specificity Analysis:** Classifies clusters as Site-Specific or Multi-Targeting based on integration site diversity.

8. **Test Validation (MAFFT):** Compares discovered test candidates against known positive controls via Pairwise and Multiple Sequence Alignments.

## Script Descriptions

### Data Acquisition & Pre-processing

* `download_genomes.py`: Automates the retrieval of whole-genome assemblies from NCBI based on a provided sampling strategy.

* `group_by_species.py`: Organizes downloaded genomes taxonomically to prepare for diversity saturation subsampling.

* `collect_lsr_genomes.py`: Parses HMMER outputs to quickly filter and retain only the genomes containing a valid PF07508 domain hit.

### Structural Analysis

* `run_mgefinder.py`: Executes the whole-genome alignment and insertion discovery workflow. Note: Interfaces directly with MGEfinder CLI tools via subprocess to maintain stability across isolated Conda environments.

* `reconstruct_att_sites.py`: Parses MGEfinder coordinates and mathematically slices the flanking DNA from the whole-genome assembly to reconstruct the attachment sites (attL and attR).

### Filtering & Database Compilation

* `quality_filter_lsrs.py`: Enforces strict quality control parameters (MGE size, distance to insertion, core homology). [Updated] Now natively deduplicates redundant structural predictions from MGEfinder and automatically extracts both the raw DNA sequence and Amino Acid sequence for all surviving candidates.

* `cluster_target_genes.py`: Prepares the genomic environment surrounding the insertion site for MMseqs2 clustering to evaluate target diversity.

* `predict_specificity.py`: Analyzes the relationship between LSR clusters and their target gene clusters to categorize recombinases as "Site-Specific" or "Multi-Targeting".

* `compile_database.py`: Aggregates all structural, specificity, and clustering metadata into the final relational database. [Updated] Features robust error-handling, missing-header detection for clustering TSVs, and safe-merging for empty outputs.

### Validation & Testing

* `validate_test_results.py` **(NEW)**: Evaluates pipeline accuracy by comparing discovered candidates in test batches against known positive control sequences (provided via CSV). Automatically calculates homology, identifies SNPs, and generates a color-coded HTML report of Pairwise and Multiple Sequence Alignments (via MAFFT) for quick visual validation.

## Usage

**Running a Test Batch**

To run the pipeline in test mode (which forces the validation rule to execute), pass a test CSV containing the expected control sequences:

`snakemake -c 4 --use-conda --config test=resources/Dn29.csv`


***Outputs:*** `results/test_validation_report.txt` and `results/test_validation_alignment.html.`

**Full Dataset Execution**

To run the standard discovery pipeline across the primary genome dataset:

`snakemake -c <threads> --use-conda`


## Conda Environments

Dependencies are strictly isolated to prevent version conflicts. Ensure conda or mamba is installed before executing Snakemake.

* envs/hmmer.yaml

* envs/prodigal.yaml

* envs/mgefinder.yaml

* envs/mmseqs2.yaml

* envs/mafft.yaml