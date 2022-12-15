# Packages Used
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

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

library(GenomicFeatures)
library(tximport)

samples <- read.table("C:/Users/emmet/OneDrive/Desktop/Masters/Bio Principles/samplelist.txt", header = FALSE)
colnames(samples) <- c("condition","run")
# rownames(samples) <- samples$condition
samples 

files <- file.path("C:/Users/emmet/OneDrive/Desktop/Masters/Bio Principles/", samples$run, "quant.sf")
names(files) <- samples$run
all(file.exists(files))
files

# prepare txdb, create map from tx to gene
txdb <- makeTxDbFromGFF("C:/Users/emmet/OneDrive/Desktop/Masters/Bio Principles/Homo_sapiens.GRCh38.101.gtf.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

#Note, you can query the txdb by keys of a keytype and get columns back. 
#To see the available keytypes: keytypes(txdb) and coluns: columns(txdb)
#keytypes(txdb)
#columns(txdb)

# tx import salmon files
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
head(txi.salmon$counts)

## imports


# run the data through DESeq2
ddsSE <- DESeqDataSetFromTximport(txi.salmon, colData =  samples, ~ condition)
dds <- ddsSE

# Check counts of data
data.frame(sort(rowSums(counts(dds))))

data.frame(counts(dds))

range(rowSums(counts(dds)))

br = c(0,100, 500, 1000, 7500000)
ranges = paste(head(br,-1), br[-1], sep=" - ")
freq   = hist(rowSums(counts(dds)), breaks=br, include.lowest=TRUE, plot=FALSE)
data.frame(range = ranges, frequency = freq$counts)

br = c(0,20, 40, 60, 80, 100, 7500000)
ranges = paste(head(br,-1), br[-1], sep=" - ")
freq   = hist(rowSums(counts(dds)), breaks=br, include.lowest=TRUE, plot=FALSE)
data.frame(range = ranges, frequency = freq$counts)

# Set condition
dds$condition <- factor(dds$condition, levels = c("TOT_0h","TOT_1h", "TOT_4h", "TOT_24h"))
dds <- DESeq(dds)

# Due to difference in large and small counts, remove counts less than 100
keep <- rowSums(counts(dds)) >= 100
dds <- dds[keep,]
nrow(dds)

# Check resutls and compute p-values
res <- results(dds)
head(res)
summary(res)

## Filter the results for significance
# Filter the results by false dicovery rate
res.05 <- results(dds, alpha = 0.05)
# Filter the results false  discovery and by log2 fold change (LFC)
LFC1 <- results(dds, lfcThreshold=1, alpha=0.05)

summary(LFC1)

ens.str <- rownames(res)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

# It is nice to have the most significant results on the top. 
# We sort the results by to pvalue
resOrdered <- res[order(res$pvalue),]
head(resOrdered)

# export to table
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]

## Visualization

#Have a at the top genes:
# Produce 6 plots on one page, with 3 rows and 2 columns 
par(mfrow = c(3, 2))
plotCounts(dds, gene = "ENSG00000134323", intgroup=c("condition"))
plotCounts(dds, gene = "ENSG00000205670", intgroup=c("condition"))
plotCounts(dds, gene = "ENSG00000273590", intgroup=c("condition"))
plotCounts(dds, gene = "ENSG00000106560", intgroup=c("condition"))
plotCounts(dds, gene = "ENSG00000158270", intgroup=c("condition"))
plotCounts(dds, gene = "ENSG00000288349", intgroup=c("condition"))


# For visualization purposes is it usually beneficial to somehow normalize or standardize 
# the data. Here we use the variance shrinkage algorithm (see DESeq2 vignette for details)
min(rowSums(counts(dds)))
max(rowSums(counts(dds)))

rlogd <- rlog(dds, blind = TRUE)
anno <- as.data.frame(colData(rlogd)[, c("run","condition")])

# We also make a list of the top 15 most variable genes, and focus the visualization on that
topVarGenes <- head(order(rowVars(assay(rlogd)), decreasing = TRUE), 15)
mat <- assay(rlogd)[topVarGenes, ]

# Produce 1 plot on page
par(mfrow = c(1,1))

pheatmap(mat)

plotPCA(rlogd, intgroup = c("run","condition"))

pcaData <- plotPCA(rlogd, intgroup = c("run","condition"), returnData = TRUE)

# # gene enrichment
# # Use dds dataset
# geneiddf <- as.data.frame(unique(names(dds)))
# geneiddf$ID <- rownames(geneiddf)
# colnames(geneiddf) <- c("gene","ID")
# genes <-  head(names(rlogd)[order(rowVars(assay(rlogd)), decreasing = TRUE)],1000)
# geneiddf <- filter(geneiddf, gene %in% genes)
# gene1 <- as.character(geneiddf$ID) 
# 
# 
# ddsdf <- as.data.frame(counts(ddsSE))
# ddsdf$gene <- rownames(ddsdf)
# ddsdf <- filter(ddsdf, gene %in% genes)
# dds.genedf <- melt(ddsdf,id = c("gene"))
# dds.genedf <- merge(dds.genedf, geneiddf, by = "gene") %>% select(-gene, -value)
# colnames(dds.genedf) <- c("gs_name","entrez_gene")
# dds.genedf$entrez_gene <- as.integer(dds.genedf$entrez_gene)
# dds.genedf$gs_name <- as.character(dds.genedf$gs_name)
# 
# 
# str(dds.genedf)
# str(msigdbr_t2g)
# str(gene1)
# str(gene)
# # Enrichment by hypergeometric test implemented in the enricher function
# res_enricher <- enricher(gene = gene1, TERM2GENE = dds.genedf)
# 
# 
# # Gene set enrichment analysis is implemented in the GSEA funcion
# res_GSEA <- GSEA(ddsdf, TERM2GENE = dds.genedf)
# 
# # The GSEA result is a special data type: gseaResult. 
# # We can convert it into a dataFrame:
# res_GSEA_df <- as.data.frame(res_GSEA)
# 
# # Visualization
# # barplot(res_GSEA, showCategory=2)
# dotplot(res_GSEA, showCategory=2)
# 
# 
# mat_gene_exp <- 
# 
# res_enricher <- enricher(gene = gene, TERM2GENE = msigdbr_t2g)