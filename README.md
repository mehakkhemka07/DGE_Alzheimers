**Differential Gene Expression Analysis in Alzheimer's Disease**

Identifying significantly dysregulated genes and biological pathways between AD and control brain tissue using GEO dataset GSE5281.


**Objective**
To identify significantly upregulated and downregulated genes between Alzheimer's disease and control brain tissue samples, and explore their biological significance through functional enrichment analysis.

**Dataset**
FieldDetailsSourceGene Expression OmnibusAccessionGSE5281Total Samples161AD / Control74 Alzheimer's / 87 ControlPlatformAffymetrix Human Genome U133 Plus 2.0 Array

**Workflow**
StepDescription01Data retrieval using GEOquery02Preprocessing and normalization03Sample grouping (Alzheimer's vs Control)04Differential expression analysis using limma (adjusted p-value, log2FC filtering)05Visualization — Volcano plot, PCA, Heatmap06Functional enrichment — GO and KEGG

**Results**
🌋 Volcano Plot
Visualizes significantly upregulated and downregulated genes between AD and control samples.
📉 PCA Plot
Shows separation between Alzheimer's and control sample groups in principal component space.
🔥 Heatmap
Displays expression patterns of top differentially expressed genes across all samples.
🧬 GO Enrichment
Highlights enriched biological processes such as neuronal signalling and cellular stress response.
🧪 KEGG Pathways
Identifies disrupted pathways related to neurodegeneration and altered metabolism.

**Output Files**
FileDescriptionDEGs.csvFiltered list of significant differentially expressed genesvolcano.pngVolcano plotPCA_plot.pngPCA visualizationheatmap.pngHeatmap of top genesGO_barplot.pngGO enrichment bar plotKEGG_dotplot.pngKEGG pathway dot plot

**Tools & Packages**

R
GEOquery — GEO data retrieval
limma — Differential expression analysis
EnhancedVolcano — Volcano plot
pheatmap — Heatmap visualization
clusterProfiler — GO and KEGG enrichment


▶️ How to Run
bash# Clone the repository
git clone https://github.com/yourusername/repo-name.git
cd repo-name
Then in R:
rsource("main_script.R")

📌 Key Insight
The analysis identifies significant gene expression changes between Alzheimer's and control brain tissue. Enriched biological processes point to disrupted neuronal function and elevated cellular stress, while KEGG results highlight affected neurodegenerative and metabolic pathways as hallmarks of disease.
