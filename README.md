# Differential Gene Expression Analysis in Alzheimer's Disease

Identifying significantly dysregulated genes and biological pathways between AD and control brain tissue using GEO dataset GSE5281.

---

## Objective

To identify significantly upregulated and downregulated genes between Alzheimer's disease and control brain tissue, and explore their biological significance through functional enrichment analysis.

---

## Dataset

| Field | Details |
|---|---|
| **Source** | [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/) |
| **Accession** | GSE5281 |
| **Samples** | 161 (outliers removed post-QC) |
| **AD / Control** | 74 / 87 |
| **Platform** | Affymetrix Human Genome U133 Plus 2.0 Array |

---

## Workflow

| Step | Description |
|---|---|
| 01 | Data retrieval using GEOquery |
| 02 | Log2 transformation and preprocessing |
| 03 | Outlier removal (PCA-based, > 3 SD on PC1/PC2) |
| 04 | Differential expression using limma (adj. p < 0.05, log2FC > 1) |
| 05 | Volcano plot, PCA, Heatmap |
| 06 | Probe-to-gene mapping via hgu133plus2.db |
| 07 | GO and KEGG enrichment using top 500 DEGs |

---

## Results

**Volcano Plot** — 16,454 DEGs identified; net upregulation bias observed in AD tissue after log2 transformation.

**PCA Plot** — Moderate AD vs control separation after outlier removal, consistent with heterogeneous bulk brain tissue.

**Heatmap** — Top 50 DEGs shown with gene symbols. Key genes include MT-ND6, STMN1, OLFM3, and FGF13.

**GO Enrichment** — Top processes: synaptic vesicle transport, aerobic respiration, mitochondrial membrane organisation, and cellular response to interferon.

**KEGG Pathways** — Top hits: Alzheimer disease, Pathways of neurodegeneration, Oxidative phosphorylation, Synaptic vesicle cycle, Proteasome.

---

## Output Files

| File | Description |
|---|---|
| `DEGs.csv` | Significant DEGs (adj. p < 0.05, log2FC > 1) |
| `volcano.png` | Volcano plot |
| `PCA_plot.png` | PCA after outlier removal |
| `heatmap.png` | Heatmap of top 50 DEGs with gene symbols |
| `GO_barplot.png` | GO enrichment bar plot |
| `GO_results.csv` | Full GO results table |
| `KEGG_dotplot.png` | KEGG pathway dot plot |
| `KEGG_results.csv` | Full KEGG results table |

---

## Tools

R, GEOquery, limma, hgu133plus2.db, EnhancedVolcano, pheatmap, clusterProfiler, org.Hs.eg.db

---

## How to Run

```bash
git clone https://github.com/yourusername/repo-name.git
cd repo-name
```

```r
source("main_script.R")
```

---

## Key Insight

The analysis identifies 16,454 dysregulated genes in AD brain tissue. Enrichment results consistently point to synaptic vesicle dysfunction, mitochondrial energy failure, and neuroinflammation as core mechanisms. The Alzheimer disease KEGG pathway ranks among the top enriched terms, directly validating the experimental findings.
