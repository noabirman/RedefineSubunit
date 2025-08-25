# RedefineSubunit: Improved Subunit Definition for CombFold/Graph-based Complex Splitting for CombFold

This repository provides an alternative way to insert protein complexes into [CombFold](https://github.com/dina-lab3D/CombFold), not just as individual chains but as **subunits defined by high-confidence structural regions**.  

It takes AlphaFold model outputs (mmCIF + confidence JSON) and:
1. **Splits chains into subunits** using plDDT thresholds and gaps.
2. **Builds graphs** where subunits are nodes and edges reflect inter-subunit PAE interactions.
3. **Merges overlapping subunits** across models into final `subunits_info.json` for CombFold.

---

## Installation  

Clone the repository and install dependencies:

```bash
git clone https://github.com/username/RedefineSubunit.git
cd RedefineSubunit
pip install -r requirements.txt
```

Dependencies include:
- 'biopython' for structure parsing
- 'networkx' for graph construction
- 'numpy' for matrix handling




**CombFold: predicting structures of large protein assemblies using a combinatorial assembly algorithm and AlphaFold2**  
Ben Shor & Dina Schneidman-Duhovny  
[https://doi.org/10.1038/s41592-024-02174-0](https://www.nature.com/articles/s41592-024-02174-0)

### About This Pipeline

Traditional AlphaFold inputs sometimes define subunits as continuous chains, missing biologically relevant domain splits or disorder. This can reduce the accuracy of complex structure prediction in large assemblies.

**RedefineSubunit** enhances the preprocessing of AlphaFold/CombFold inputs by leveraging sequence disorder prediction (via IUPred3) to better define chain boundaries. This results in more accurate subunit definitions and improved performance in combinatorial folding of complexes.

The pipeline handles:

- Conversion of AlphaFold JSON inputs to per-chain FASTA sequences
- Disorder prediction using IUPred3
- Splitting long chains at predicted disordered regions
- Renaming subunits to compact identifiers
- Generating clean AlphaFold3-style JSON files
- Running Multiple Sequence Alignments (MSAs)
- Filtering and merging complex graphs

### Key Features

- **Improved subunit splitting** using disorder regions, to better guide CombFold inference
- Support for **AlphaFold3 JSON formats**
- Integration with **IUPred3** for intrinsic disorder prediction
- Optional MSA generation and chain filtering steps
- Scripting support for SLURM-based clusters

### Getting Started

1. Install dependencies (Python 3.8+, Biopython, NumPy, etc.)
2. Prepare AlphaFold JSON input files for each complex
3. Use `preprocess/preprocess.py` to run the main pipeline
4. Use `run_af.sh` or `run_af_missing.sh` to execute AlphaFold or CombFold predictions
5. Analyze or visualize predictions using available scripts

### Example Entry Point

```bash
python preprocess/preprocess.py \
    --input_dir path/to/json_inputs/ \
    --output_dir path/to/processed_outputs/
```

This will:

- Convert JSON to FASTA
- Run IUPred3 and split sequences based on disorder
- Regenerate JSON using updated splits
- Prepare data for MSA and CombFold

### Output

- Cleaned JSON-ready input files for CombFold
- Per-chain FASTA sequences
- Mappings between original and updated chain IDs
- Evaluation-ready structured predictions

### Reference

If you use this pipeline, please cite:

**CombFold: predicting structures of large protein assemblies using a combinatorial assembly algorithm and AlphaFold2**  
Ben Shor & Dina Schneidman-Duhovny  
https://www.nature.com/articles/s41592-024-02174-0

---

For more about CombFold, see the official repository:  
https://github.com/dina-lab3D/CombFold
