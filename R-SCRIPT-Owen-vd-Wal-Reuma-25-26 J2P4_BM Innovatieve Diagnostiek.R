#================================
#Working derectory
#================================
setwd("C:/Users/owenv/OneDrive - NHL Stenden/Project REUMA leerjaar 2 periode 4")

#================================
#Pakages 
#================================
install.packages('BiocManager')
BiocManager::install('Rsubread')
library(Rsubread)
#===========================================================
#het referentie genoom van de mens indexeren
#Met de volgende code kun je het referentiegenoom indexeren.
#===========================================================
buildindex(
  basename = 'ref_homosapian',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

#================================
#alignen van de samples
#================================
align.819_31_normal<-align(index = "ref_homosapian", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "819_31_normal.BAM")
align.820_15_normal<-align(index = "ref_homosapian", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "820_15_normal.BAM")
align.828_31_normal<-align(index = "ref_homosapian", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "828_31_normal.BAM")
align.831_42_normal<-align(index = "ref_homosapian", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "831_42_normal.BAM")
align.879_54_Rheuma<-align(index = "ref_homosapian", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "879_54_Rheuma.BAM")
align.980_66_Rheama<-align(index = "ref_homosapian", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "980_66_Rheama.BAM")
align.986_60_Rheuma<-align(index = "ref_homosapian", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "986_60_Rheuma.BAM")
align.988_59_Rheuma<-align(index = "ref_homosapian", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "988_59_Rheuma.BAM")

#===================================================================
# Laad Rsamtools voor sorteren en indexeren (dowloaden indien nodig)
#===================================================================
BiocManager::install('Rsamtools')
library(Rsamtools)

#================================
# Bestandsnamen van de monsters
#================================
samples <- c('819_31_normal', '820_15_normal', '828_31_normal', '831_42_normal', '879_54_Rheuma', '980_66_Rheama', '986_60_Rheuma', '988_59_Rheuma' )

#======================================================
# Voor elk monster: sorteer en indexeer de BAM-file
# Sorteer BAM-bestanden
#======================================================

lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})
# Indexeer de gesorteerde BAM-file
lapply(samples, function(s) {indexBam(file = paste0(s, '.sorted.BAM'))
})

library(Rsubread)
#============================================
#het maken van een count matrix voor 1 sample
#============================================
count_matrix <- featureCounts(
  files = "819_31_normal.BAM",
  annot.ext = "genomic.mens.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.featureType = "gene",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)
#================================
#voor meerdere samples
#================================

allsamples <- c("819_31_normal.BAM", "820_15_normal.BAM", "828_31_normal.BAM", "831_42_normal.BAM", "879_54_Rheuma.BAM", "980_66_Rheama.BAM", "986_60_Rheuma.BAM","988_59_Rheuma.BAM")
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "genomic.mens.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "gene", 
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE
)

allsamples <- c("819_31_normal.BAM", "820_15_normal.BAM", "828_31_normal.BAM", "831_42_normal.BAM", "879_54_Rheuma.BAM", "980_66_Rheama.BAM", "986_60_Rheuma.BAM", "988_59_Rheuma.BAM")
#================================
#count matrix opslaan onder counts
#================================
counts <- count_matrix$counts
head(counts)

#================================
#koppen veranderen van de matrix
#================================
colnames(counts) <- c("819_31_normal", "820_15_normal", "828_31_normal", "831_42_normal", "879_54_Rheuma", "980_66_Rheama", "986_60_Rheuma", "988_59_Rheuma")
head(counts)

#========================
# WORKDIRECTORY
#========================
setwd("C:/Users/owenv/OneDrive - NHL Stenden/Project REUMA leerjaar 2 periode 4")

#========================
# PACKAGES (Slimme installatie)
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

#========================
# OUTPUT CHECK
#========================
cat("GO terms gevonden:", nrow(as.data.frame(go_results)), "\n")
cat("KEGG terms gevonden:", nrow(as.data.frame(kegg_results)), "\n")