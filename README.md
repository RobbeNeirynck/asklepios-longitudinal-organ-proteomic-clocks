# Organ Proteomic Aging Clocks

This repository contains the key methodological scripts used for our manuscript, "Longitudinal dynamics of organ-specific proteomic aging clocks over a decade of midlife"

A preprint is available at: https://www.biorxiv.org/content/10.64898/2026.02.17.706320v1

Developed and tested on Windows 11, Python Version 3.10 and R Version 4.3.0.

## Python Environment Setup
Expected installation time: ~15 minutes.

```bash
conda env create -f Setup/environment.yml
conda activate organ_clocks_env
```
## R Environment Setup
Expected installation time: ~15 minutes.

```bash
Rscript Setup/install_dependencies.R
```

## Data
Individual-level Asklepios data are not publicly available. To test the scripts and verify the expected data structure, a permuted subset of the original data is provided in the `Data/` folder.

### SomaScan QC
The SomaScan CV values used for QC (`SomaScanCVs.csv`) are derived from:
*Candia, J., et al. Assessment of Variability in the Plasma 7k SomaScan Proteomics Assay. Sci Rep 7, 14248 (2017).*

### GTEx Enrichment Analysis
The GTEx enrichment analysis requires raw GTEx GCT files not included here. **GTEx v8 gene reads** (GCT format) can be downloaded from the [GTEx Portal](https://gtexportal.org/). 

## Usage

### 1. Organ Clock Training
* **Expected Runtime:** ~1h 
* **Output:** Predictions (train and test) and model coefficients 

To train the organ-specific aging clocks using the provided protein data, metadata, and aptamer-organ map:

```bash
python Organ_Proteomic_Aging_Clocks.py
```

### 2. Structural Equation Modeling (SEM)
The SEM analysis is a two-step process that starts from the organ age predictions from the test set:

#### Step 1: Prepare Age Gaps
* **Runtime:** ~1 minute
* **Output:** Sex-adjusted organ age gaps

```bash
Rscript SEM_prepare_agegaps.R
```

#### Step 2: Run FGES Network
* **Runtime:** ~5 minutes 
* **Output:** Directed network graph image

```bash
python SEM_FGES_network.py
```

### 3. GTEx Enrichment
* **Runtime:** ~10 minutes
* **Output:** Organ-enrichment summary table for SomaScan 7K aptamers needed to create aptamer-organ map

To run the tissue enrichment analysis, first download the publicly available GTEx data and put the `.gct` files in `Data/GTEx/`. 

```bash
Rscript GTEx_SomaScan_Enrichment_Analysis.R
```

### Using Your Own Data
To run these scripts on your own dataset, format your protein and meta data to match the structure in `Data/`.

## License
This project is covered under the MIT License.

## Citation
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
