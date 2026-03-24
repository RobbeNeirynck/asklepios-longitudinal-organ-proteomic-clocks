# Organ Proteomic Aging Clocks

This repository contains the key methodological scripts written for our manuscript "Longitudinal dynamics of organ-specific proteomic aging clocks over a decade of midlife"

A preprint is available at https://www.biorxiv.org/content/10.64898/2026.02.17.706320v1

Developed and tested on Windows 11, Python Version 3.10 and R Version 4.3.0.

### Python Environment Setup

```bash
# Create the environment from the provided file
conda env create -f Setup/environment.yml

# Activate the environment
conda activate organ_clocks_env
```

### R Environment Setup
You can install the required R packages by running the provided setup script:
```bash
Rscript Setup/install_dependencies.R
```

### Data
Individual-level Asklepios data is not publicly available. To test scripts and verify data structure, permuted versions of the original data are provided in the `Data/` folder.

#### SomaScan QC
The SomaScan CV values used for QC (`SomaScanCVs.csv`) are derived from:
*Candia, J., et al. Assessment of Variability in the Plasma 7k SomaScan Proteomics Assay. Sci Rep 7, 14248 (2017).*

#### GTEx Enrichment Analysis
 The GTEx enrichment analysis requires raw GTEx GCT files not included here

Download the **GTEx v8 gene reads** (GCT format) from the [GTEx Portal](https://gtexportal.org/). Place the downloaded `.gct` files in `Data/GTEx/`.

### Usage

#### 1. Organ Clock Training
To train the organ-specific aging clocks using the provided protein data and metadata:
```bash
python Scripts/Organ_Proteomic_Aging_Clocks.py
```

#### 2. Structural Equation Modeling (SEM)
The SEM analysis is a two-step process:

1. **Prepare Age Gaps:** first, generate the sex-adjusted age gaps.
```bash
Rscript Scripts/SEM_prepare_agegaps.R
```
This script reads the organ age predictions and metadata, adjusts for sex, and prepares the data for SEM analysis.

2. **Run FGES Network:** then, run the SEM FGES network script.
```bash
python Scripts/SEM_FGES_network.py
```

#### 3. GTEx Enrichment
To run the tissue enrichment analysis:
```bash
Rscript Scripts/GTEx_SomaScan_Enrichment_Analysis.R
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


