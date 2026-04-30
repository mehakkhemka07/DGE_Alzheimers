if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c("dplyr", "ggplot2", "pheatmap")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

bioc_pkgs <- c("GEOquery", "limma", "EnhancedVolcano",
               "clusterProfiler", "org.Hs.eg.db",
               "hgu133plus2.db", "AnnotationDbi")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# 1. Load Libraries

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(AnnotationDbi)

# 2. Load Data

gse  <- getGEO("GSE5281", GSEMatrix = TRUE)
eset <- gse[[1]]
expr  <- exprs(eset)
pheno <- pData(eset)


# 3. Define Groups

group <- factor(ifelse(grepl("control", pheno$title, ignore.case = TRUE),
                       "control", "affected"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)


# 4. Differential Expression
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(
  affected_vs_control = affected - control,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
deg <- topTable(fit2,
                coef          = "affected_vs_control",
                number        = Inf,
                adjust.method = "BH")
# 5. Filter DEGs
# Full filtered list (for heatmap, volcano, CSV)
deg_filtered <- deg %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

message("Total DEGs after filtering: ", nrow(deg_filtered))
write.csv(deg_filtered, "DEGs.csv")

# Top 500 by adj p-value for enrichment
# (16k+ genes = ~46% of universe; enrichment needs a focused list)
deg_for_enrich <- deg_filtered %>%
  arrange(adj.P.Val) %>%
  head(500)

message("Genes submitted for enrichment: ", nrow(deg_for_enrich))

# 6. Volcano Plot
png("volcano.png", 800, 600)
EnhancedVolcano(deg,
                lab      = rownames(deg),
                x        = 'logFC',
                y        = 'adj.P.Val',
                pCutoff  = 0.05,
                FCcutoff = 1,
                title    = 'Volcano Plot: Affected vs Control')
dev.off()
message("volcano.png saved.")

# 7. PCA Plot
pca_result <- prcomp(t(expr), scale. = TRUE)

pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  Group = group
)

ggsave("PCA_plot.png",
       ggplot(pca_df, aes(PC1, PC2, color = Group)) +
         geom_point(size = 3, alpha = 0.8) +
         theme_minimal() +
         labs(title = "PCA Plot", x = "PC1", y = "PC2"))

message("PCA_plot.png saved.")

# 8. Heatmap
if (nrow(deg_filtered) > 0) {
  top_genes    <- intersect(rownames(expr), rownames(deg_filtered))[1:min(50, nrow(deg_filtered))]
  heatmap_data <- t(scale(t(expr[top_genes, ])))
  heatmap_data <- na.omit(heatmap_data)
  
  annotation_col <- data.frame(Group = group)
  rownames(annotation_col) <- colnames(expr)
  
  png("heatmap.png", 1000, 800)
  pheatmap(heatmap_data,
           annotation_col = annotation_col,
           show_colnames  = FALSE,
           main           = "Top 50 DEGs Heatmap")
  dev.off()
  message("heatmap.png saved.")
} else {
  message("No genes for heatmap.")
}

# 9. Probe -> Symbol -> Entrez
# --- Background universe: ALL expressed probes ---
all_symbols_raw <- mapIds(hgu133plus2.db,
                          keys      = rownames(expr),
                          column    = "SYMBOL",
                          keytype   = "PROBEID",
                          multiVals = "first")
all_symbols <- na.omit(all_symbols_raw)

universe_ids <- bitr(all_symbols,
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = org.Hs.eg.db)
message("Universe size (Entrez IDs): ", nrow(universe_ids))

# --- DEG list: top 500 focused gene set ---
deg_symbols_raw <- mapIds(hgu133plus2.db,
                          keys      = rownames(deg_for_enrich),
                          column    = "SYMBOL",
                          keytype   = "PROBEID",
                          multiVals = "first")
deg_symbols <- na.omit(deg_symbols_raw)
message("DEG symbols after probe mapping: ", length(deg_symbols))

gene_ids <- bitr(deg_symbols,
                 fromType = "SYMBOL",
                 toType   = "ENTREZID",
                 OrgDb    = org.Hs.eg.db)
message("DEG Entrez IDs for enrichment: ", nrow(gene_ids))

# 10. GO Enrichment
if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
  
  ego <- enrichGO(gene          = gene_ids$ENTREZID,
                  universe      = universe_ids$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  
  message("GO terms found: ", nrow(ego@result))
  
  if (!is.null(ego) && nrow(ego@result) > 0) {
    png("GO_barplot.png", 1000, 800)
    print(barplot(ego, showCategory = 15,
                  title = "GO Biological Process Enrichment"))
    dev.off()
    message("GO_barplot.png saved.")
    
    write.csv(as.data.frame(ego), "GO_results.csv", row.names = FALSE)
  } else {
    message("No significant GO terms found.")
  }
  
} else {
  message("Gene ID conversion failed — check probe mapping.")
}
# 11. KEGG Enrichment
if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
  
  kegg <- tryCatch(
    enrichKEGG(gene              = gene_ids$ENTREZID,
               universe          = universe_ids$ENTREZID,
               organism          = 'hsa',
               pvalueCutoff      = 0.05,
               qvalueCutoff      = 0.2,
               use_internal_data = FALSE),
    error = function(e) {
      message("KEGG online fetch failed: ", e$message,
              "\nRetrying with internal data...")
      enrichKEGG(gene              = gene_ids$ENTREZID,
                 universe          = universe_ids$ENTREZID,
                 organism          = 'hsa',
                 pvalueCutoff      = 0.05,
                 use_internal_data = TRUE)
    }
  )
  
  message("KEGG pathways found: ", nrow(kegg@result))
  
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    png("KEGG_dotplot.png", 1000, 800)
    print(dotplot(kegg, showCategory = 15,
                  title = "KEGG Pathway Enrichment"))
    dev.off()
    message("KEGG_dotplot.png saved.")
    
    write.csv(as.data.frame(kegg), "KEGG_results.csv", row.names = FALSE)
  } else {
    message("No significant KEGG pathways found.")
  }
  
} else {
  message("Gene ID conversion failed — check probe mapping.")
}

message("=== Analysis complete ===")