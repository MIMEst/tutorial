# RNA-seq downstream analysis
# DESeq2
library(ggplot2)
library(tidyr)
library(dplyr)
library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(fgsea)
library(GseaVis)
library(purrr)
library(clusterProfiler)

# 1. Establishment of DEGs (Differentially Expressed Genes)
# expected_counts: merged_RSEM_expected_counts, metadata: information of each samples
expected_counts <- expected_counts[rowSums(expected_counts) >= 10,] 
dds <- DESeqDataSetFromMatrix(countData = expected_counts, colData = metadata, design = ~group)
dds <- DESeq(dds, fitType = "local", minReplicatesForReplace = Inf)
res <- results(dds, alpha=0.05, contrast=c("group", "numerator", "denominator")
summary(res)
res_df <- as.data.frame(res)
res_df$mgi_symbol <- mapIds(org.Mm.eg.db, keys=rownames(df), keytype="ENSEMBL", column="SYMBOL")


# 2. PCA plot (checking whether my samples are separated as we expected?!) 
# what is PCA plot?
# https://ddongwon.tistory.com/114 (PCA plot 설명)
mat <- assay(vst(dds, blind = FALSE))
pca <- prcomp(t(mat))
nudge <- position_nudge(y = -0.05)
autoplot(pca, data = metadata, colour = "group", size= 5, alph=0.75) + theme_classic()

# finding out top contributing genes for each PCA
loadings_df <- data.frame(
  ENSEMBL = rownames(pca$rotation),
  PC1 = pca$rotation[, "PC1"],
  row.names = NULL
) %>%
  mutate(
    abs_PC1  = abs(PC1),
    direction = ifelse(PC1 >= 0, "Positive", "Negative")
  )
sym <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = loadings_df$ENSEMBL,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
loadings_df$SYMBOL <- unname(sym)
loadings_df$gene <- ifelse(is.na(loadings_df$SYMBOL), loadings_df$ENSEMBL, loadings_df$SYMBOL)

top_n <- 30
plot_df <- loadings_df %>%
  arrange(desc(abs_PC1)) %>%
  slice_head(n = top_n) %>%
  mutate(gene = factor(gene, levels = rev(gene)))

ggplot(plot_df, aes(x = gene, y = PC1, fill = direction)) +
  geom_col(width = 0.8) +
  coord_flip() +
  labs(
    title = paste0("Top ", top_n, " PC1-contributing genes"),
    x = NULL,
    y = "PC1 loading"
  ) +
  theme_classic() +
  theme(legend.position = "top")

pc1_go <- enrichGO(gene = rownames(loadings_df)[1:50], OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",ont = "ALL")
barplot(pc1_go)

# 3. Hierarchical clustering
rld <- rlog(dds, blind=TRUE)
rld_cor <- cor(rld_mat)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color = heat.colors, border_color=NA, fontsize = 10,fontsize_row = 10, height=20)

# below code is wrong! you need change!!
sampleDists <- dist( t( assay(rld) ) )
as.matrix( sampleDists )[ 1:3, 1:3 ]
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment, 
   rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

# 4. Dispersion plot
plotDispEsts( dds, ylim = c(1e-6, 1e1))

# 5. p-value histogram
hist( res$pvalue, breaks=20, col="grey" )
  
# 6. volcano plot 
# plot showing differentially expressed genes (looks like exploding volcano!)
plot_volcano <- function(res_df, title_text) {
  df <- res_df
  # Replace NA padj with 1
  df$padj[is.na(df$padj)] <- 1
  # Replace 0 padj with floor(min nonzero padj / 100)
  min_nonzero <- min(df$padj[df$padj > 0], na.rm = TRUE)
  df$padj[df$padj == 0] <- min_nonzero / 100
  
  # Significance flag
  df$sig <- ifelse(df$padj < 0.05 & df$log2FoldChange > 1, "Up", ifelse(df$padj < 0.05 & df$log2FoldChange < -1, "Down", "NS"))
  
  # Filter genes to label
  label_df_up <- df %>% 
    filter(padj < 0.05 & log2FoldChange > 2) %>% 
    arrange(padj) %>% 
    slice_head(n = 50)
  label_df_down <- df %>%
    filter(padj < 0.05 & log2FoldChange < -2) %>%
    arrange(padj) %>% 
    slice_head(n= 50)
  label_df <- rbind(label_df_up,label_df_down)
  
  # Plot
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.7, size = 2.5) + 
    geom_text_repel(data = label_df, aes(label = symbol), size = 3.5, max.overlaps = 100, box.padding = 0.005) +
    scale_color_manual(values = c("Up" = "#7B3014", "Down" = "#26456E", "NS" = "gray80")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top") +
    labs(
      title = title_text,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance")
}
plot_volcano(res_df)

# 7. Heatmap
library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
     trace="none", dendrogram="column", 
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
     ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
        colData(rld)$treatment ] )
        
res_sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
sig_genes <- rownames(res_sig)
top_genes <- head(
  rownames(res_df[order(res$padj), ]),
  50
)
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))
pheatmap(
  mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8
)


# 8. Gene Set Enrichment Analysis
# https://happip-jh.tistory.com/3 (gsea 설명) 
# https://www.gsea-msigdb.org/gsea/login.jsp (gsea gmt file download)
pathway.hallmark <- read.gmt("msigdb.v2025.1.Mm.symbols.gmt")
res_df <- res_df[order(-res_df$log2FoldChange),]
res_df <- na.omit(res_df)
gsea <- res_df$log2FoldChange
names(gsea) <- res_df$symbol
gsea_res <- GSEA(geneList = gsea, TERM2GENE = pathway.hallmark, minGSSize = 1, pvalueCutoff = 1)
gseaNb(object = gsea_res,
       geneSetID = geneSetID, 
       subPlot = 2,
       termWidth = 35,
       legend.position = c(0.8,0.8),
       addGene = gene)
      
# 9. KEGG pathway analysis
# ranked gene list
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)
entrez <- bitr(
  rownames(res_sig),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)
kegg_enrich <- enrichKEGG(
  gene         = entrez$ENTREZID,
  organism     = "mmu",   # mouse
  pvalueCutoff = 0.05
)
kegg_gsea <- gseKEGG(
  geneList     = gene_list,
  organism     = "mmu",
  pvalueCutoff = 0.05
)
dotplot(kegg_enrich, showCategory = 20)
emapplot(kegg_enrich)
library(pathview)
pathview(
  gene.data  = gene_list,
  pathway.id = "mmu00010",  # Glycolysis
  species    = "mmu"
)

