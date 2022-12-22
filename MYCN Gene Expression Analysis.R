# Uncomment to Install Packages if Necessary
# install.packages('tidyverse')
# install.packages('clusterProfiler')
# install.packages('msigdbr')
# install.packages('enrichplot')
# install.packages('msigdbr')
# install.packages('pheatmap')
# install.packages('DESeq2')
# install.packages('ggplot2')
# install.packages('reshape')
# install.packages('DESeq2')
# install.packages("org.Hs.eg.db")
# install.packages('GenomicFeatures')
# install.packages('tximport')
# install.packages('diffEnrich')
# install.packages('GenomicFeatures')
# install.packages('tximport')
# install.packages('TxDb.Hsapiens.UCSC.hg38.knownGene')
# install.packages('heatmaply')

# Read in Packages Used
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(msigdbr)
library(pheatmap)
library(DESeq2)
library(ggplot2)
library(reshape)
library(DESeq2)
library("org.Hs.eg.db")
library(GenomicFeatures)
library(tximport)
library(diffEnrich)
library(GenomicFeatures)
library(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(heatmaply)
library(clusterProfiler)

# set working directory to Downloads.
path = "C:/Users/emmet/OneDrive/Desktop/Masters/Bio Principles/"
setwd(path)

# Read sample names from folder.
samples <- read.table("./samplelist.txt", header = FALSE)
colnames(samples) <- c("condition","run")

# Create  the path to read each file. 
files <- file.path(path,samples$run, "quant.sf")
names(files) <- samples$run

# Read a transcript to gene file
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k2 <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k2, keytype = "TXNAME", columns = "GENEID")
tx2gene <- tx2gene[!is.na(tx2gene$GENEID),] 
#Note, you can query the txdb by keys of a keytype and get columns back. 
#To see the available keytypes: keytypes(txdb) and coluns: columns(txdb)
#keytypes(txdb)
#columns(txdb)

symbols <- mapIds(org.Hs.eg.db, keys = tx2gene[,2], column="SYMBOL", keytype="ENTREZID", multiVals="first")
ensembles <-  mapIds(org.Hs.eg.db, keys = tx2gene[,2], column="ENSEMBL", keytype="ENTREZID", multiVals="first")

# We can update our tx2gene into tx2symb object.
# Counts with Entrez Id
txis_f <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE)

# tx import salmon files
# txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
# head(txi.salmon$counts)

## imports
# Make a DDs object
dds <- DESeq(DESeqDataSetFromTximport(txis_f, colData =  samples, ~ condition))
counts = assay(dds)
# Set condition
dds$condition <- factor(dds$condition, levels = c("TOT_0h","TOT_1h", "TOT_4h", "TOT_24h"))

# Check counts of data
rowSums(counts)
range(rowSums(counts))

br = c(0,100, 500, 1000, 1909973)
ranges = paste(head(br,-1), br[-1], sep=" - ")
freq   = hist(rowSums(counts(dds)), breaks=br, include.lowest=TRUE, plot=FALSE)
data.frame(range = ranges, frequency = freq$counts)

br = c(0,20, 40, 60, 80, 100, 1909973)
ranges = paste(head(br,-1), br[-1], sep=" - ")
freq   = hist(rowSums(counts(dds)), breaks=br, include.lowest=TRUE, plot=FALSE)
data.frame(range = ranges, frequency = freq$counts)


# Due to difference in large and small counts, remove counts less than 100
keep <- rowSums(counts(dds)) > 100
dds <- dds[keep,]
nrow(dds)

# Compute Results
res <- results(dds)
# Lookup rowname, symbol and ensembl to the results data
ens.str <- rownames(res)
res$symbol <- mapIds(org.Hs.eg.db,keys=ens.str,column="SYMBOL",keytype="ENTREZID",multiVals="first")
res$ensembl <- mapIds(org.Hs.eg.db,keys=ens.str,column="ENSEMBL",keytype="ENTREZID",multiVals="first")

# It is nice to have the most significant results on the top. 
# Filter for log2FoldChange greater than 1 and adjusted p-value < 0.05
# We sort the results by adjusted pvalue
resfil <- as.data.frame(res) %>% filter(log2FoldChange > 1 & padj < 0.05)
resfil <- resfil[order(resfil$padj),]

## Visualization
#Have a at the 7 genes based on p-value
# Produce 7 count plots on 1 page of normalised dds, with 3 rows and 3 columns 
# For visualization purposes is it usually beneficial to somehow normalize or standardize 
# the data. Here we use the variance shrinkage algorithm (see DESeq2 vignette for details)
par(mfrow = c(2,2))

rlogdds <- rlog(dds)
rlogdds$condition <- factor(rlogdds$condition, levels = c("TOT_0h","TOT_1h", "TOT_4h", "TOT_24h"))
assayrlogdds <- assay(rlogdds)
sigGenes <- as.character(rownames(resfil))
  
for(i in 1:length(sigGenes)){
  plot(rlogdds$condition, assayrlogdds[which(rownames(rlogdds)==sigGenes[i]),],main = sigGenes[i], pch = 19,xlab="Times", ylab="Normalized Counts")
  
}

# We also make a list of the top 25 most variable genes, and focus the visualization on that
#topVarGenes <- head(order(rowVars(assayrlogdds), decreasing = TRUE), 25)
topVarGenes <- rownames(resfil)
mat <- assay(rlogdds)[topVarGenes, ]
anno <- as.data.frame(colData(dds)[, c("condition","run")])
# Produce 1 plot on page
par(mfrow = c(1,1))
# Produce Heatmap
pheatmap(mat, legend = TRUE,annotation_col=anno)


# Examine coefficient overtime and produce heatmap
ddscoef <- as.data.frame(coef(dds))
ddscoef$gene <- rownames(ddscoef)
ddscoef <- ddscoef[order(-genelogdf$log2FoldChange),] %>%
  filter(gene %in% sigGenes) %>%
  select(-gene)
pheatmap(ddscoef, legend = TRUE)

# Conduct PCA on normalised data
plotPCA(rlogdds, intgroup = c("run","condition"))
plotPCA(rlogdds, intgroup = c("run","condition"), returnData = TRUE)

### we can do now some visualization from top differentially expressed
# Enrichment Analysis

# First create geneList named number vector

genelogdf <- as.data.frame(results(dds)) %>% 
  select(log2FoldChange)
genelogdf$gene <- rownames(genelogdf)
genelogdf <- genelogdf[order(-genelogdf$log2FoldChange),] 
geneList <- setNames(genelogdf$log2FoldChange,genelogdf$gene)
gene <- na.omit(names(geneList)[abs(geneList) > 2])

# Get all human gene sets:
m_df <- msigdbr(species = "Homo sapiens")

# View what collections there are
collections <- msigdbr_collections()

# Select (subset) the Hallmark collection and selct necessary columns
m_df_H <- m_df[m_df$gs_cat=="H",]
msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

# Enrichment by hypergeometric test implemented in the enricher function
res_enricher <- enricher(gene = gene, TERM2GENE = msigdbr_t2g)
res_enricher

# Gene set enrichment analysis is implemented in the GSEA funcion
res_GSEA <- GSEA(geneList, TERM2GENE = msigdbr_t2g,pvalueCutoff = 0.05)
res_GSEA


# Pathway enrichment by time

Kegg_c_up = list()
Kegg_c_down = list()

for (i in 1:length(coefs)){
  if (!is.null(res_coefs[[coefs[[i]]]])){
    gs = rownames(res_coefs[[coefs[i]]])
    up_genes = gs[res_coefs[[coefs[i]]]$log2FoldChange>1]
    down_genes = gs[res_coefs[[coefs[i]]]$log2FoldChange<1]
    Kegg_c_up[[coefs[i]]] = enrichKEGG(gene =up_genes, organism = 'hsa', pvalueCutoff = 0.05)
    Kegg_c_down[[coefs[i]]] = enrichKEGG(gene =down_genes, organism = 'hsa', pvalueCutoff = 0.05)
  }
  
}

gs
# Up and down regulated genes
up_genes
down_genes
# Examine up and down expressed pathways
Kegg_c_up
Kegg_c_down

