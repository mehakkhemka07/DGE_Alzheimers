# ===============================
# 0. Install Packages (only if missing)
# ===============================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("dplyr", "ggplot2", "pheatmap")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

bioc_pkgs <- c("GEOquery", "limma", "EnhancedVolcano",
               "clusterProfiler", "org.Hs.eg.db",
               "hgu133plus2.db")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# ===============================
# 1. Load Libraries
# ===============================
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)

# ===============================
# 2. Load Data
# ===============================
gse <- getGEO("GSE5281", GSEMatrix = TRUE)
eset <- gse[[1]]

expr <- exprs(eset)
pheno <- pData(eset)

# ===============================
# 3. Define Groups
# ===============================
group <- factor(ifelse(grepl("control", pheno$title, ignore.case = TRUE),
                       "control", "affected"))

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# ===============================
# 4. Differential Expression
# ===============================
fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts(
  affected_vs_control = affected - control,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

deg <- topTable(fit2,
                coef = "affected_vs_control",
                number = Inf,
                adjust.method = "BH")

# ===============================
# 5. Filter DEGs (REALISTIC)
# ===============================
deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

write.csv(deg_filtered, "DEGs.csv")

# ===============================
# 6. Volcano Plot
# ===============================
png("volcano.png", 800, 600)

EnhancedVolcano(deg,
                lab = rownames(deg),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Volcano Plot')

dev.off()

# ===============================
# 7. PCA Plot
# ===============================
expr_t <- t(expr)

pca <- prcomp(expr_t, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Group = group
)

p <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA Plot")

ggsave("PCA_plot.png", plot = p)

# ===============================
# Heatmap (FIXED)
# ===============================
if (nrow(deg_filtered) > 0) {
  
  # Ensure valid genes
  valid_genes <- intersect(rownames(expr), rownames(deg_filtered))
  top_genes <- valid_genes[1:min(50, length(valid_genes))]
  
  heatmap_data <- expr[top_genes, ]
  
  # Scale safely
  heatmap_data <- t(scale(t(heatmap_data)))
  heatmap_data <- na.omit(heatmap_data)
  
  # Fix annotation
  annotation_col <- data.frame(Group = group)
  rownames(annotation_col) <- colnames(expr)
  
  png("heatmap.png", 800, 600)
  
  pheatmap(heatmap_data,
           annotation_col = annotation_col,
           main = "Top DEGs")
  
  dev.off()
  
} else {
  message("No genes available for heatmap.")
}

# ===============================
# 9. PROBE → GENE SYMBOL MAPPING (CRITICAL FIX)
# ===============================
probe_ids <- rownames(deg_filtered)

gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

gene_symbols <- na.omit(gene_symbols)

# ===============================
# 10. SYMBOL → ENTREZ
# ===============================
gene_ids <- bitr(gene_symbols,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# ===============================
# 11. GO ENRICHMENT
# ===============================
if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
  
  ego <- enrichGO(gene = gene_ids$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)
  
  if (!is.null(ego) && nrow(ego@result) > 0) {
    
    png("GO_barplot.png", 800, 600)
    barplot(ego, showCategory = 10)
    dev.off()
    
  } else {
    message("No GO enrichment terms found.")
  }
  
} else {
  message("Gene ID conversion failed.")
}

# ===============================
# 12. KEGG ENRICHMENT
# ===============================
if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
  
  kegg <- enrichKEGG(gene = gene_ids$ENTREZID,
                     organism = 'hsa')
  
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    
    png("KEGG_dotplot.png", 800, 600)
    dotplot(kegg)
    dev.off()
    
  } else {
    message("No KEGG pathways found.")
  }
  
} else {
  message("Gene ID conversion failed.")
}