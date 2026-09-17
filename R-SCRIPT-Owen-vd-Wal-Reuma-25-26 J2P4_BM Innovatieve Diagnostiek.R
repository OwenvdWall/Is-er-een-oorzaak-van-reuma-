
#========================
# WORKDIRECTORY
#========================
setwd("C:/Users/owenv/OneDrive - NHL Stenden/Project REUMA leerjaar 2 periode 4")

#========================
# PACKAGES 
#========================
required_packages <- c("DESeq2", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "pathview", "dplyr", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new_packages)
}

library(DESeq2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(pathview)
library(ggplot2)

#========================
# COUNT MATRIX INLADEN
#========================
counts <- read.table("count_matrix_RA.txt", header = TRUE, row.names = 1)
colnames(counts)

#========================
# METADATA 
#========================
treatment <- c("normal", "normal", "normal", "normal",
               "Rheuma", "Rheuma", "Rheuma", "Rheuma")

coldata <- data.frame(treatment = factor(treatment))
rownames(coldata) <- colnames(counts)

#========================
# DESEQ ANALYSE
#========================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ treatment)
dds <- DESeq(dds)

resultaten <- results(dds, contrast = c("treatment", "Rheuma", "normal"))
resultaten <- resultaten[order(resultaten$padj), ]

#========================
# RESULTATEN OPSLAAN
#========================
write.csv(as.data.frame(resultaten), file = "Reuma_results.csv")

#========================
# EXTRA: AANTAL SIGNIFICANTE GENEN
#========================
upregulated <- sum(resultaten$padj < 0.05 & resultaten$log2FoldChange > 1, na.rm = TRUE)
downregulated <- sum(resultaten$padj < 0.05 & resultaten$log2FoldChange < -1, na.rm = TRUE)

cat("Upregulated genen:", upregulated, "\n")
cat("Downregulated genen:", downregulated, "\n")

#========================
# VOLCANO PLOT (Aangepast voor overzichtelijkheid)
#========================

#=======================================================================
# VOLCANO PLOT (Automatisch opslaan in hoge resolutie)
#=======================================================================

# 1. Definieer de specifieke genen die je gelabeld wilt hebben
key_genes <- c("SRGN", "PTGFR", "BCL2A1", "DVL2", "ADAMDEC1", "IGHV1-69")

# 2. Open de PNG-writer (overschrijft je oude bestand met de schone versie)
png(
  "Volcano_Plot.png",
  width = 8,
  height = 8,
  units = "in",
  res = 300  # Dwingt de hoge 300 DPI resolutie af
)

# 3. Teken de plot in de PNG-writer
EnhancedVolcano(
  resultaten,
  lab = rownames(resultaten),
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano Plot",
  subtitle = "Rheuma versus Normaal",
  pCutoff = 0.05,
  FCcutoff = 1,
  selectLab = key_genes, 
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "black",
  labSize = 4.0,
  legendPosition = "right"
)

# 4. Sluit de writer af
dev.off()

# 5. Toon de plot ook direct in je RStudio venster
EnhancedVolcano(
  resultaten,
  lab = rownames(resultaten),
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano Plot",
  subtitle = "Rheuma versus Normaal",
  pCutoff = 0.05,
  FCcutoff = 1,
  selectLab = key_genes, 
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "black",
  labSize = 4.0,
  legendPosition = "right"
)

# ========================
# PCA PLOT
# ========================
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "treatment") + ggtitle("PCA Plot: Rheuma versus Normaal")

# Toon in RStudio
print(pca_plot)

# Sla op als PNG
png("PCA_plot_Rheuma_vs_Normaal.png", width = 8, height = 6, units = "in", res = 300)
print(pca_plot)
dev.off()

#========================
# SIGNIFICANTE GENEN
#========================
sig_genes <- resultaten %>%
  as.data.frame() %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write.csv(sig_genes, "Significante_genen.csv")

gene_list <- rownames(sig_genes)
if(length(gene_list) == 0){
  stop("Geen significante genen gevonden. Pas eventueel je cutoffs aan.")
}

#========================
# ID CONVERSIE (SYMBOL → ENTREZ)
#========================
gene_df <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

if(nrow(gene_df) == 0){
  stop("Geen gene mapping gevonden. Controleer of de counts matrix SYMBOLs (gennamen) gebruikt.")
}

# ========================
# GO ANALYSE
# ========================
go_results <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

write.csv(as.data.frame(go_results), "GO_results.csv")

# Plot aanmaken
go_plot <- dotplot(go_results, showCategory = 10, title = "Gene Ontology Enrichment Analyse")
print(go_plot)

png("GO_Enrichment.png", width = 10, height = 8, units = "in", res = 300)
print(go_plot)
dev.off()

# ========================
# KEGG ANALYSE
# ========================
kegg_results <- enrichKEGG(
  gene = gene_df$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(kegg_results), "KEGG_results.csv")

kegg_plot <- dotplot(kegg_results, showCategory = 10, title = "KEGG Pathway Enrichment Analyse")
print(kegg_plot)

png("KEGG_Enrichment.png", width = 10, height = 8, units = "in", res = 300)
print(kegg_plot)
dev.off()

#========================
# PATHVIEW 
#========================
# Zorg dat we alleen genen pakken die succesvol gekoppeld zijn naar Entrez
rownames(sig_genes) <- rownames(sig_genes) # Voor de zekerheid

# Maak een schone vector op basis van de gemapte tabel
gene_vector <- sig_genes[gene_df$SYMBOL, "log2FoldChange"]
names(gene_vector) <- gene_df$ENTREZID

# Verwijder eventuele NA's die ontstaan zijn
gene_vector <- gene_vector[!is.na(names(gene_vector))]

# Genereer de Pathview diagrammen (deze worden direct in je working directory opgeslagen als .png)
pathview(gene.data = gene_vector, pathway.id = "04670", species = "hsa") # Leukocyte transendothelial migration
pathview(gene.data = gene_vector, pathway.id = "04660", species = "hsa") # T cell receptor signaling pathway

