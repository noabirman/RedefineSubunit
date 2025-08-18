# RedefineSubunit - Comprehensive Protein Complex Prediction and Analysis Pipeline

This repository contains a suite of scripts designed to facilitate the preprocessing, Multiple Sequence Alignment (MSA) generation, AlphaFold 3 inference, and post-prediction analysis for protein complexes. The pipeline integrates tools like IUPred3 for protein disorder prediction and orchestrates AlphaFold 3 runs to predict complex structures.

## Overview

The project provides a comprehensive workflow for handling protein sequences, preparing them for AlphaFold 3, generating necessary MSAs, executing AlphaFold 3 predictions, and performing various checks and analyses on the resulting complex structures. It focuses on modularity, allowing different stages of the prediction pipeline to be run and monitored independently.

## Features

*   **IUPred3 Integration**: Utilizes IUPred3 for the prediction of protein intrinsic disorder and ANCHOR2 for protein binding regions. This helps in identifying disordered regions for intelligent sequence splitting.
*   **Flexible Preprocessing**:
    *   Converts AlphaFold 3 input JSON files to multi-FASTA format and then splits them into individual FASTA files per chain.
    *   Runs IUPred3 on each sequence to obtain disorder scores.
    *   **Intelligent Sequence Splitting**: Splits long sequences (default >1000 residues) only at predicted disordered regions, avoiding cuts within ordered protein blocks.
    *   **Chain ID Renaming**: Renames protein chain IDs to shorter alphabetical labels (e.g., A, B, AA) and generates a mapping file for traceability.
    *   **MSA Input Preparation**: Splits the processed AlphaFold 3 JSON into individual per-chain JSON files suitable for MSA generation.
*   **MSA Generation**: Orchestrates the generation of Multiple Sequence Alignments using AlphaFold 3's `run_alphafold.py` in data pipeline-only mode, supporting both initial MSA generation (`msa.sh`) and resuming for missing MSAs (`msa_missing.sh`).
*   **Pairwise MSA Creation**: Generates pairwise MSA files from individual chain MSAs, including self-pairs for original multi-chain proteins, facilitating inter-chain interaction prediction.
*   **AlphaFold 3 Inference Orchestration**: Scripts to run AlphaFold 3 inference on prepared input directories, including a script to continue runs for missing predictions (`run_af_missing.sh`).
*   **Complex Analysis and Validation**:
    *   **Missing Pairs and Shared Chains Detection**: Identifies pairs that have not been processed by AlphaFold 3 and finds common chains between different subunits.
    *   **TM-score Comparison**: Calculates and compares TM-scores and RMSD values for predicted complex structures against reference PDB files, and can optionally compare against external benchmark scores.
    *   **PAE to Subunit Clustering**: Script for clustering based on Predicted Aligned Error (PAE) matrices to define subunits.
*   **Utility Scripts**: Includes various helper scripts for tasks such as cleaning and organizing directories, updating JSON file names, and extracting specific data.

## Directory Structure

The project follows a modular directory structure:

*   `iupred3/`: Contains the core IUPred3 library and main script for disorder prediction.
*   `junk/`: A collection of utility scripts, possibly experimental or dataset-specific, for tasks like checking complexes, processing PAE data, and calculating TM-scores.
*   `preprocess/`: Scripts for preparing input data, including FASTA/JSON conversions, sequence splitting, MSA generation, and chain ID renaming.
    *   `preprocess/junk/`: Additional preprocessing utilities such as `change_split_names.py` for adjusting JSON filenames, `clean_folders.py` for managing duplicate prediction folders, and `synapse_extract_pairs.py`/`synapse_preprocess.py` for specific data extraction (e.g., from a "Synapse" dataset).
*   `tests/`: Scripts for testing various components and transformations, such as `check_transformation.py` for validating graph connectivity of subunits.
*   Top-level scripts: Various scripts for specific tasks like `cif_to_pdb.py`, `complex_graph.py`, `filter_high_subunits.py`, `merge_graphs.py`, `partial_fold.py`, `run_cf.sh`, `second_preprocess.py`, `tm_score_new.py`, and `vizualization_plots.py`.

## Getting Started

### Prerequisites

*   **Conda**: The pipeline extensively uses `conda` for environment management.
*   **AlphaFold 3**: The scripts are designed to interact with a local AlphaFold 3 installation, specifically leveraging its data pipeline and inference capabilities. Paths to AlphaFold 3 binaries (`run_alphafold.py`) and databases are hardcoded in some shell scripts.
*   **HMMER Tools**: `jackhmmer`, `hmmbuild`, and `hmmsearch` binaries from the HMMER suite are required for MSA generation.
*   **MMalign**: The `MMALIGN_PATH` variable in `junk/tm_score.py` points to an MMalign executable, which is used for calculating TM-scores and RMSD.

### Installation

1.  **Clone the Repository**:
    ```bash
    git clone <repository_url>
    cd <repository_name>
    ```
2.  **Conda Environment Setup**:
    The `.sh` scripts indicate the use of a `conda` environment named `alphafold3-conda`. Ensure this environment is set up and activated before running the scripts.
    ```bash
    # Example (adjust based on your AlphaFold3 installation)
    conda activate /cs/usr/bshor/sci/installations/af3_variations/deepmind/localalphafold3/alphafold3-conda
    ```
3.  **Path Configuration**:
    Verify and adjust the hardcoded paths within the `.sh` scripts (e.g., `run_alphafold.py`, `jackhmmer`, `db_dir`, `model_dir`, `MMALIGN_PATH`) to match your local installation of AlphaFold 3, HMMER, and MMalign.

## Usage

The general workflow for processing a protein complex often involves these steps:

1.  **Prepare Initial Input JSON**:
    Start with a JSON file containing protein sequences, typically in an AlphaFold 3 input format. This file should be placed in a directory for the specific complex (e.g., `my_complex/subunits_info.json`).

2.  **Run Preprocessing Pipeline**:
    The `preprocess.py` script automates several initial steps, including FASTA conversion, IUPred3 analysis, intelligent sequence splitting, and chain ID renaming.
    ```bash
    python3 preprocess/preprocess.py /path/to/my_complex/subunits_info.json
    ```
    This will generate:
    *   `sequences.fasta`: Full FASTA file.
    *   `input_fastas/`: Directory with individual FASTA files per chain.
    *   `iupred_outputs/`: Directory with IUPred3 prediction results for each chain.
    *   `iupred_split_sequences.fasta`: FASTA file with sequences potentially split based on disorder.
    *   `iupred_split_mapping.json`: Mapping of original sequences to their split fragments.
    *   `af3_input.json`: AlphaFold 3 input JSON based on split sequences.
    *   `af3_input_renamed.json`: AlphaFold 3 input JSON with renamed chain IDs.
    *   `chain_id_mapping.json`: Mapping of new chain IDs to original protein IDs and sequence ranges.
    *   `msa_inputs/`: Directory containing individual AlphaFold 3 input JSON files for each chain (with renamed IDs).

3.  **Generate MSAs (Multiple Sequence Alignments)**:
    Use the `msa.sh` script to run the AlphaFold 3 data pipeline for MSA generation on the `msa_inputs` directory created in the previous step.
    ```bash
    sbatch preprocess/msa.sh /path/to/my_complex/msa_inputs/
    ```
    This will create an `msa_output/` directory containing the MSA results for each chain. If you need to re-run or continue for missing MSAs, use `msa_missing.sh`.

4.  **Create Pairwise MSAs**:
    After MSA generation, run `msa_to_pairwise.py` to prepare pairwise MSA JSON files for AlphaFold 3 inference. These files define the pairs of chains for complex prediction.
    ```bash
    python3 preprocess/msa_to_pairwise.py /path/to/my_complex/msa_output/ /path/to/my_complex/chain_id_mapping.json /path/to/my_complex/subunits_info.json
    ```
    This creates an `msa_pairs/` directory with `_` formatted JSON files (e.g., `A_B.json`).

5.  **Run AlphaFold 3 Inference**:
    Execute AlphaFold 3 inference using `run_af.sh` on the `msa_pairs` directory.
    ```bash
    sbatch preprocess/run_af.sh /path/to/my_complex/msa_pairs/
    ```
    The predictions will be saved in an `af_pairs/` directory. For resuming runs, `run_af_missing.sh` can be used.

6.  **Perform Complex Analysis (Optional)**:
    *   **Check Missing Pairs and Shared Chains**:
        ```bash
        python3 junk/check_complexes.py /path/to/my_complex_root_directory/
        ```
        This script will generate `summary.json` in the specified root directory, detailing missing AlphaFold predictions and shared chains between subunits.
    *   **Calculate TM-scores**:
        ```bash
        python3 junk/tm_score.py /path/to/my_complex_root_directory/ /path/to/ben_scores.json # (optional: for external comparison)
        ```
        This script downloads reference PDBs, calculates TM-scores and RMSD for predicted models, and saves them to `tm_score.txt` within each complex directory. It can also compare against a provided JSON of benchmark scores.

## Citations

The IUPred3 component of this pipeline is based on:

*   **IUPred3**:
    Gábor Erdős, Mátyás Pajkos, Zsuzsanna Dosztányi
    *Nucleic Acids Research 2021, Submitted*
*   **IUPred2A**:
    Balint Meszaros, Gabor Erdos, Zsuzsanna Dosztanyi
    *Nucleic Acids Research 2018;46(W1):W329-W337.*
