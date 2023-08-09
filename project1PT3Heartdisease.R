# --------------------------------
# Loading R packages


library(DESeq2)
library(pheatmap)
library(annotate)
library(vctrs)
library(tidyverse)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ensembldb)
library(pathview)
library(AnnotationHub)
library(biomaRt)
library(org.Mm.eg.db)
library(ggplot2)
library(cowplot)
library(utils)
library(utils, lib.loc = "C:/Program Files/R/R-4.1.1/library")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("DelayedArray")
BiocManager::install("GenomeInfoDb")
BiocManager::install("DESeq2")
aBiocManager::install("DOSE")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ensembldb")
BiocManager::install("pathview")
BiocManager::install("AnnotationHub")
BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("ggnewscale")
BiocManager::install("Rattus.norvegicus")


install.packages(pheatmap)


# --------------------------------
# Data description and importation
setwd("C:/Users/melle/Downloads")
data <- read.table("GSE183025_All.HTSeq.counts1.txt", header = TRUE, row.names = 1)
data <- as.matrix(data)
meta <- read.table("heartdiseaseinfosample.txt", header = TRUE)
rownames(meta) <- meta$labels


# Match the metadata and counts data
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

View(counts(dds))

## Create DESEq2 object

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ group)


## Generate the normalized counts

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)


write.table(normalized_counts, file=" normalized_counts.txt", sep="\t", quote=F,
            col.names=NA)

## QC for DE analysis using DESeq2
## Transform normalized counts using the rlog function

rld <- rlog(dds, blind=TRUE)


#Principal components analysis (PCA
plotPCA(rld, intgroup="group")

## Extract the rlog matrix from the object
rld_mat <- assay(rld)


rld_cor <- cor(rld_mat)

pheatmap(rld_cor)


## Differential expression analysis with DESeq2

dds <- DESeq(dds)

plotDispEsts(dds)

res <- results(dds)

summary(res)

DESeq2::plotMA(res, ylim=c(-2,2)) ###


# Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
all_genes <- as.character(rownames(res))

# Extract significant results
signif_res <- res[res$padj < 0.05 & !is.na(res$padj), ]
signif_genes <- as.character(rownames(signif_res))


##Annotation to DESeq2 results.


## Functional analysis with clusterProfiler
## Over-representation analysis with clusterProfiler
## Run GO enrichment analysis
ego <- enrichGO(gene = signif_genes,
                universe = all_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)

head(ego)



# Output results from GO analysis to a table
cluster_summary <- data.frame(ego)


## Visualizing clusterProfiler results
dotplot(ego, showCategory=50)
emapplot(ego, showCategory=50)



# To color genes by log2 fold changes
signif_res_lFC <- signif_res$log2FoldChange

cnetplot(ego,
         categorySize="pvalue",
         showCategory = 5,
         foldChange= signif_res_lFC,
         vertex.label.font=6)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters="ensembl_gene_id",
               attributes=c("ensembl_gene_id", "entrezgene_id"),
               values= all_genes,
               mart=mart)
indNA = which(is.na(genes$entrezgene_id))
genes_noNA <- genes[-indNA,]
indnodup = which(duplicated(genes_noNA$ entrezgene_id) == F)
genes_noNA_nodup <- genes_noNA[indnodup,]
lFC <- res$log2FoldChange[-indNA]
lFC <- lFC[indnodup]
names(lFC) <- genes_noNA_nodup$entrezgene_id
# Sort fold changes in decreasing order
lFC <- sort(lFC, decreasing = TRUE)


## Perform the GSEA using KEGG gene sets:
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- getBM(filters="ensembl_gene_id",
               attributes=c("ensembl_gene_id", "entrezgene_id"),
               values= all_genes,
               mart=mart)
indNA = which(is.na(genes$entrezgene_id))
genes_noNA <- genes[-indNA,]
indnodup = which(duplicated(genes_noNA$ entrezgene_id) == F)
genes_noNA_nodup <- genes_noNA[indnodup,]
lFC <- res$log2FoldChange[-indNA]
lFC <- lFC[indnodup]
names(lFC) <- genes_noNA_nodup$entrezgene_id


# Sort fold changes in decreasing order
lFC <- sort(lFC, decreasing = TRUE)
gseaKEGG <- gseKEGG(geneList = lFC,
                    organism = "mmu",
                    nPerm = 1000, # default number permutations
                    minGSSize = 5, # minimum gene set size
                    pvalueCutoff = 0.1, # padj cutoff value
                    verbose = FALSE)

# Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

##volcano Plot
pValue <- 0.05
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(), ylim=c(0,2), xlim=c(-1,1), main="VolcanoPlot", xlab="log2(fold-change)", ylab="-log10(adjusted p-value", pch=20, cex=0.3)





