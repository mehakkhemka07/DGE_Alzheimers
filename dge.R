if (!require(&quot;BiocManager&quot;, quietly = TRUE))
install.packages(&quot;BiocManager&quot;)
BiocManager::install(c(&quot;GEOqurey&quot;,&quot;limma&quot;), ask = FALSE, update = FALSE)
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)
gse &lt;- getGEO(&quot;GSE5281&quot;, GSEMatrix = TRUE)
length (gse)
eset &lt;- gse[[1]]
expr &lt;- exprs(eset)
pheno &lt;-pData(eset)
dim(expr)
head(pheno)
colnames(pheno)
table(pheno$title) 
group &lt;- factor(ifelse(grepl(&quot;control&quot;, pheno$title, ignore.case =  TRUE), &quot;affected&quot;, &quot;control&quot;))
design &lt;- model.matrix(~ 0 + group)
colnames(design) &lt;- levels(group)
design
fit &lt;- lmFit(expr, design)
contrast.matrix &lt;- makeContrasts(
affected_vs_control = affected - control,
levels = design)
fit2 &lt;- contrasts.fit(fit, contrast.matrix) 
fit2 &lt;-eBayes(fit2)
deg &lt;- topTable(fit2,
coef = &quot;affected_vs_control&quot;,
number = Inf,
adjust.method = &quot;BH&quot;)
deg_filtered &lt;- deg[
deg$adj.P.Val &lt; 0.005 &amp; abs(deg$logFC) &gt; 50, ]
write.csv(deg_filtered, &quot;GSE5281_DEGs.csv&quot;)
deg_results &lt;- read.csv(&quot;GSE5281_DEGs.csv&quot;, stringsAsFactors = FALSE)
annotations &lt;- read.csv(&quot;annotations.csv&quot;, stringsAsFactors = FALSE)
library(dplyr)
colnames(deg_results)

colnames(annotations)
annotated_results &lt;- deg_results %&gt;%
left_join(annotations, by =c(&quot;X&quot; = &quot;ID&quot;))
write.csv(annotated_results, &quot;annotated_DEG_results.csv&quot;, row.names = FALSE)
library(ggplot2)
library(dplyr)
BiocManager::install(&quot;EnhancedVolcano&quot;)
library(EnhancedVolcano)
BiocManager::install(&quot;STRINGdb&quot;)
library(STRINGdb)
BiocManager::install(&quot;WGCNA&quot;, force = TRUE)
library(WGCNA)
data&lt;- read.csv(file.choose() , header = TRUE)
head(data)
EnhancedVolcano(data,
lab = data$Gene.Symbol, # columnwith gene names
x = &#39;logFC&#39;, #column with log2 fold change
y = &#39;adj.P.Val&#39;, #column with p-value
pCutoff = 0.05, # Adjust p-value threshold
FCcutoff = 2, #Adjust fold-change threshold
title = &#39;volvano plot of DEGs&#39;,
xlab = &#39;Log2 Fold Change&#39;,
ylab = &#39;p-value&#39;,
col = c(&#39;yellow&#39;,&#39;orange&#39;,&#39;red&#39;, &#39;green&#39;), #custoom colour
legendPosition = &#39;right&#39;)