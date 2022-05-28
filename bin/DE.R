#!/usr/bin/env Rscript

library( "optparse" )

parser <- OptionParser()

parser <- add_option(parser, c( "-i", "--input" ),  action = "store", 
                     type="character", dest = "input_filename", 
                     help = "Input file (counts data from featureCounts).")
parser <- add_option(parser, c( "-o", "--output-de" ), action = "store", 
                     type="character", dest = "de_filename", 
                     help = "Filename for output (DE genes with Padj values).")
parser <- add_option(parser, c( "-v", "--output-vulcano" ), action = "store", 
                     type="character", dest = "vulcano_filename", 
                     help = "Filename of eps file for Vulcano plot.")
options = parse_args(parser, positional_arguments = TRUE)

#setwd("/home/egle/Desktop/SystemBiology/Transcriptomics/snake/collibri-snakemake/outputs/STAR/all/Collibri")
#data = read.csv('counts.txt', header=TRUE, sep="\t", skip=1, 
#                row.names = "Geneid")
data = read.csv(options$options$input_filename, header=TRUE, sep="\t", skip=1, 
                row.names = "Geneid")
cts = as.matrix(data[,-c(1,2,3,4,5)])
samples = colnames(cts)
condition = vector()
for (i in (1:length(samples))){
    parts = strsplit(samples[i], split="\\.")
    condition[i] = parts[[1]][4]
    samples[i] = paste(parts[[1]][4], parts[[1]][5], parts[[1]][7], sep="-")
    title = parts[[1]][3]
}

# UHRR is a sample originating from several cancerous cell lines,
# HBR represents “normal” sample --> creating coldata based on that:

condition = gsub("^UHRR$", "cancerous", condition)
condition = gsub("^HBR$", "normal", condition)
colnames(cts) = samples
coldata = data.frame(condition, row.names = samples)
coldata$condition = factor(coldata$condition)

cts_file = gsub("/([A-Z]*[a-z]*)*\\.[a-z]*$", "/cts.csv", options$options$de_filename)
write.csv(cts,cts_file, row.names = TRUE)
coldata_file = gsub("cts.csv$", "coldata.csv", cts_file)
write.csv(coldata,coldata_file, row.names = TRUE)

# Make sure that sample names match in counts and column datatables:
if (! all(rownames(coldata) == colnames(cts))){
  cts = cts[, rownames(coldata)]
}

suppressPackageStartupMessages(library(DESeq2))
dds = DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# While it is not necessary to pre-filter low count genes before running 
# the DESeq2 functions, there are two reasons which make pre-filtering useful: 
# by removing rows in which there are very few reads, we reduce the memory 
# size of the dds data object, and we increase the speed of the transformation 
# and testing functions within DESeq2. Here we perform a minimal pre-filtering 
# to keep only rows that have at least 10 reads total. 

keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

# By default, R will choose a reference level for factors based on alphabetical 
# order. Then, if you never tell the DESeq2 functions which level you want to 
# compare against (e.g. which level represents the control group), the 
# comparisons will be based on the alphabetical order of the levels. A solution
# is to explicitly set the factors levels. In order to see the change of 
# reference levels reflected in the results names, you need to either run DESeq 
# or nbinomWaldTest/nbinomLRT after the re-leveling operation. Setting the 
# factor levels can be done using factor:

dds$condition = factor(dds$condition, levels = c("normal","cancerous"))
dds = DESeq(dds)
res = results(dds)

# Shrinkage of effect size (LFC estimates) is useful for visualization and 
# ranking of genes. We provide the dds object and the name or number of the 
# coefficient we want to shrink, where the number refers to the order of the 
# coefficient as it appears in resultsNames(dds).

# resultsNames(dds)
resLFC = lfcShrink(dds, coef=2, type="apeglm")
summary(resLFC)

#The output should be:
#  --  a listof DE genes wit Padj values. -> chose to output all result data.
# DE_genes_with_Padj_values = subset(res, select = c("padj") )
# DE_genes_with_Padj_values
resLFC = resLFC[order(resLFC$padj),]

write.csv(resLFC,options$options$de_filename, row.names = TRUE)

#  --  Vulcanplot.
colors = resLFC$log2FoldChange
colors[resLFC$log2FoldChange<=-1] = "blue"
colors[resLFC$log2FoldChange>-1] = "black"
colors[resLFC$log2FoldChange>=1] = "blue"

setEPS()
postscript(options$options$vulcano_filename, width = 9, height = 6)
plot(resLFC$log2FoldChange, -log10(resLFC$padj), col=colors, panel.first=grid(), 
     main=paste("Volcano plot", title, sep=" - "),
     xlab="Effect size: log2(FoldChange)", ylab="Significance: -log10(padj)", 
     pch=20, cex=0.6)
a = 0.1 # Threshold on the adjusted p-value - alpha
abline(v=0)
abline(h=-log10(a), col="red")

# Selecting values where log2FoldChange is bigger than 2.5 and adjusted p-value
# is lower than threshold
selected = abs(resLFC$log2FoldChange) > 2.5 & resLFC$padj < a 
# Adding names to these values
text(resLFC$log2FoldChange[selected], -log10(resLFC$padj)[selected], 
     lab=rownames(resLFC)[selected], cex=0.4)
dev.off()


