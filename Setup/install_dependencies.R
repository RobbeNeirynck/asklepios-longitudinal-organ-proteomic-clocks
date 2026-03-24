# CRAN Packages
cran_pkgs <- c("tidyverse", "writexl", "broom", "glue", "devtools", "BiocManager")
install.packages(cran_pkgs[!(cran_pkgs %in% installed.packages()[,"Package"])], repos = "https://cloud.r-project.org")

# Bioconductor Packages
BiocManager::install(c("edgeR", "cmapR", "SomaScan.db"), update = FALSE, ask = FALSE)
