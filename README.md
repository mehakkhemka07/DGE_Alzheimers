# Differential Gene Expression Analysis in Alzheimer's Disease

Identification of differentially expressed genes between Alzheimer's disease and control brain tissue samples using publicly available microarray data (GSE5281) and R.

---

## Dataset

| Field | Details |
|---|---|
| **Source** | [NCBI Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/) |
| **Accession** | GSE5281 |
| **Disorder** | Alzheimer's disease |
| **Organism** | *Homo sapiens* |
| **Sample type** | Post-mortem brain tissue |
| **Platform** | GPL570 — Affymetrix Human Genome U133 Plus 2.0 Array |
| **Samples** | 161 total (74 affected, 87 control) |

---

## Objectives

1. Retrieve and preprocess the microarray dataset GSE5281 from GEO.
2. Perform differential gene expression (DEG) analysis between Alzheimer's disease and control samples using R.
3. Identify significantly upregulated and downregulated genes associated with Alzheimer's disease.
4. Visualize DEG results with a volcano plot.

---

## Experimental Design

GSE5281 follows a **case–control** design. Brain tissue RNA from post-mortem subjects (Alzheimer's vs. cognitively normal) was hybridized onto Affymetrix arrays. Expression levels were normalized and samples grouped into disease and control categories for comparative analysis.

---

## Methods

### Tools & Packages

- **R** with Bioconductor
- [`GEOquery`](https://bioconductor.org/packages/GEOquery/) — dataset retrieval
- [`limma`](https://bioconductor.org/packages/limma/) — linear modeling and empirical Bayes moderation
- [`EnhancedVolcano`](https://bioconductor.org/packages/EnhancedVolcano/) — volcano plot visualization
- `ggplot2`, `ggrepel`, `dplyr` — data manipulation and plotting

### Analysis Pipeline

1. Download dataset from GEO using `getGEO("GSE5281")`
2. Extract expression matrix and phenotype data
3. Classify samples into `affected` and `control` groups
4. Build design matrix and fit linear model with `lmFit`
5. Define contrasts (`affected vs. control`) and apply empirical Bayes with `eBayes`
6. Extract DEGs using `topTable` (Benjamini–Hochberg FDR correction)
7. Filter by **adjusted p-value < 0.005** and **|log2 fold change| > 50**
8. Annotate results and generate volcano plot

---

## Outputs

| File | Description |
|---|---|
| `GSE5281_DEGs.csv` | Filtered list of significant DEGs |
| `annotated_DEG_results.csv` | DEGs annotated with gene symbols |
| Volcano plot | Visual summary of expression changes (log2FC vs. adjusted p-value) |

---

## Usage

### Prerequisites

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "EnhancedVolcano"))
install.packages(c("ggplot2", "ggrepel", "dplyr"))
```

### Run the Analysis

Open the R script in RStudio or run from the command line:

```bash
Rscript deg_analysis.R
```

To generate the volcano plot, load `annotated_DEG_results.csv` when prompted and the plot will render using `EnhancedVolcano`.

---

## Repository Structure

```
.
├── deg_analysis.R            # Main R script
├── annotations.csv           # Probe-to-gene symbol annotation file
├── GSE5281_DEGs.csv          # Filtered DEG results (generated)
├── annotated_DEG_results.csv # Annotated DEG results (generated)
└── README.md
```

---

## Significance Thresholds

| Parameter | Threshold |
|---|---|
| Adjusted p-value | < 0.005 |
| \|log2 fold change\| | > 50 |
| FDR method | Benjamini–Hochberg |

> **Note:** The log2FC threshold of 50 is intentionally stringent for this assignment and may differ from thresholds used in publication-grade analyses (typically 1–2).


**Mehak Khemka**
