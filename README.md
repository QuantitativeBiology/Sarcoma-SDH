# Sarcoma-SDH

**Project Synopsis:**
This repository explores the expression patterns of succinate dehydrogenase (SDH) genes in sarcomas, with a focus on SDHB. Our findings suggest that **SDHB is overexpressed in undifferentiated pleomorphic sarcomas (UPS)** and that its overexpression potentially correlates with worse survival outcomes.

---

## Contents
This repository contains multiple files, including:
- **R scripts** for data processing, differential expression analysis, and visualization
- **Documentation** describing methods and results
- **Supporting files** for figures and tables

---

## Datasets
Two main cohorts were analyzed:
1. **IPO Cohort**: Institutional data (DDLPS, LMS, and UPS sarcoma subtypes).  
   - *Availability*: The IPO data can be provided upon request (due to privacy/ethical considerations).
2. **TCGA-SARC Cohort**: Publicly available data from The Cancer Genome Atlas for sarcoma (SARC) samples.  
   - DDLPS, LMS, and UPS subtypes were also analyzed in this cohort.

**Key Findings**:
- Overexpression of **SDH** genes (particularly SDHB) correlated with survival outcomes in the IPO cohort.
- Differential gene expression (DGE) analysis was confirmed in the **TCGA-SARC** cohort, supporting the role of SDH overexpression in sarcoma progression.

---

## Methods

### R Packages Used
A variety of **R** packages were utilized for data manipulation, statistical analysis, and visualization:

- **Data Manipulation & Cleaning**  
  `dplyr`, `readxl`, `tidyr`, `tidyverse`

- **Differential Expression Analysis**  
  `limma`, `edgeR`

- **Survival Analysis**  
  `survival`, `survminer`

- **Visualization**  
  `EnhancedVolcano`, `ggplot2`, `ggpubr`, `ComplexHeatmap`, `gplots`, `forestploter`, `enrichplot`

- **Gene Set Enrichment & Pathway Analysis**  
  `GSEABase`, `fgsea`, `msigdbr`, `clusterProfiler`

- **Other Utilities**  
  `corto`, `rempsyc`

*(Depending on your specific R environment, you may need to install and load all of these packages to reproduce the analyses.)*

---

## Usage
1. **Clone the Repository**  
   ```bash
   git clone https://github.com/YourUsername/Sarcoma-SDH.git
