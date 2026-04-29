# Differential Gene Expression Analysis in Alzheimer's Disease

> Identifying significantly dysregulated genes and biological pathways between AD and control brain tissue using GEO dataset GSE5281.

---

## Objective

To identify significantly upregulated and downregulated genes between Alzheimer's disease and control brain tissue samples, and explore their biological significance through functional enrichment analysis.

---

## Dataset

| Field | Details |
|---|---|
| **Source** | [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) |
| **Accession** | GSE5281 |
| **Total Samples** | 161 |
| **AD / Control** | 74 Alzheimer's / 87 Control |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array |

---

## Workflow

| Step | Description |
|---|---|
| 01 | Data retrieval using **GEOquery** |
| 02 | Preprocessing and normalization |
| 03 | Sample grouping (Alzheimer's vs Control) |
| 04 | Differential expression analysis using **limma** (adjusted p-value, log2FC filtering) |
| 05 | Visualization — Volcano plot, PCA, Heatmap |
| 06 | Functional enrichment — GO and KEGG |

---

## Results

### 🌋 Volcano Plot
Visualizes significantly upregulated and downregulated genes between AD and control samples.

### 📉 PCA Plot
Shows separation between Alzheimer's and control sample groups in principal component space.

### 🔥 Heatmap
Displays expression patterns of top differentially expressed genes across all samples.

### 🧬 GO Enrichment
Highlights enriched biological processes such as neuronal signalling and cellular stress response.

### 🧪 KEGG Pathways
Identifies disrupted pathways related to neurodegeneration and altered metabolism.

---

## Output Files

| File | Description |
|---|---|
| `DEGs.csv` | Filtered list of significant differentially expressed genes |
| `volcano.png` | Volcano plot |
| `PCA_plot.png` | PCA visualization |
| `heatmap.png` | Heatmap of top genes |
| `GO_barplot.png` | GO enrichment bar plot |
| `KEGG_dotplot.png` | KEGG pathway dot plot |

---

## Tools & Packages

- **R**
- **GEOquery** — GEO data retrieval
- **limma** — Differential expression analysis
- **EnhancedVolcano** — Volcano plot
- **pheatmap** — Heatmap visualization
- **clusterProfiler** — GO and KEGG enrichment

---

## How to Run

```bash
# Clone the repository
git clone https://github.com/yourusername/repo-name.git
cd repo-name
```

Then in R:

```r
source("main_script.R")
```

---

## Key Insight

The analysis identifies significant gene expression changes between Alzheimer's and control brain tissue. Enriched biological processes point to disrupted neuronal function and elevated cellular stress, while KEGG results highlight affected neurodegenerative and metabolic pathways as hallmarks of disease.
