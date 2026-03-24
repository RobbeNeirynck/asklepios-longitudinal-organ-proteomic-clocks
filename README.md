# Organ Proteomic Aging Clocks

This repository contains the key methodological scripts written for our manuscript "Longitudinal dynamics of organ-specific proteomic aging clocks over a decade of midlife"

A preprint is available at https://www.biorxiv.org/content/10.64898/2026.02.17.706320v1

Developed and tested on Windows 11, Python Version 3.10 and R Version 4.3.0.

### Python Environment Setup
Expected installation time: ~10-15 minutes.

```bash
# Create the environment from the provided file
conda env create -f Setup/environment.yml

# Activate the environment
conda activate organ_clocks_env
```

### R Environment Setup
Expected installation time: ~15-20 minutes.
You can install the required R packages by running the provided setup script:
```bash
Rscript Setup/install_dependencies.R
```

### Data
Individual-level Asklepios data is not publicly available. To test scripts and verify data structure, a permuted subset of the original data is provided in the `Data/` folder.

#### Using Your Own Data
To run these scripts on your own dataset, format your files to match the structure in `Data/`:
1. **Metadata.csv**: Must contain `SampleID`, `Age`, `Sex`, `ParticipantID`, and `Round`.
2. **Protein Dataset**: Columns should be `SampleID` followed by protein abundances.
3. **Organ Mapping**: A CSV file mapping proteins to organs (`SeqId` and `Organ` columns).

Ensure these files are placed in the `Data/` folder or update the file paths at the top of each script.


#### SomaScan QC
The SomaScan CV values used for QC (`SomaScanCVs.csv`) are derived from:
*Candia, J., et al. Assessment of Variability in the Plasma 7k SomaScan Proteomics Assay. Sci Rep 7, 14248 (2017).*

#### GTEx Enrichment Analysis
 The GTEx enrichment analysis requires raw GTEx GCT files not included here

Download the **GTEx v8 gene reads** (GCT format) from the [GTEx Portal](https://gtexportal.org/). Place the downloaded `.gct` files in `Data/GTEx/`.

### Usage
All scripts are designed to be run from the `Scripts/` directory.

#### 1. Organ Clock Training
**Expected Runtime:** ~5-10 minutes (on demo data).
**Output:** Trained models and predictions in `Results/`.

To train the organ-specific aging clocks using the provided protein data and metadata:
```bash
cd Scripts
python Organ_Proteomic_Aging_Clocks.py
```

#### 2. Structural Equation Modeling (SEM)
The SEM analysis is a two-step process:

1. **Prepare Age Gaps:**
    *   **Goal:** Generate sex-adjusted age gaps from organ age predictions.
    *   **Runtime:** < 1 minute.
    *   **Output:** `Results/agegap_sex_adjusted.csv`. This file is the required input for the next step.

```bash
cd Scripts  # if not already in Scripts/
Rscript SEM_prepare_agegaps.R
```


2. **Run FGES Network:**
    *   **Goal:** Infer causal structure using the Fast Greedy Equivalence Search (FGES) algorithm.
    *   **Runtime:** ~1-5 minutes (depending on bootstrap iterations).
    *   **Output:** Network graph structure and edge tables in `Results/`.

```bash
python SEM_FGES_network.py
```

#### 3. GTEx Enrichment
**Expected Runtime:** ~5 minutes.
**Output:** Enrichment plots and tables in `Results/GTEx_Enrichment/`.

To run the tissue enrichment analysis:
```bash
cd Scripts # if not already in Scripts/
Rscript GTEx_SomaScan_Enrichment_Analysis.R
```

### License
This project is covered under the MIT License.

### Citation
If you use this code or data in your research, please cite our preprint:

```bibtex
@article{Neirynck2026,
  title={Longitudinal dynamics of organ-specific proteomic aging clocks over a decade of midlife},
  author={Neirynck, Robbe E. and Chirinos, Julio A. and Van Damme, Menno and Coussement, Louis and Segers, Patrick and De Buyzere, Marc L. and Rietzschel, Ernst R. and De Meyer, Tim},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.02.17.706320v1}
}
```


