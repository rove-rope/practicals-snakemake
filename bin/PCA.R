#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library( "optparse" )

parser <- OptionParser()

parser <- add_option(parser, c( "-i", "--input1" ),  action = "store", 
                     type="character", dest = "input1", 
                     help = "Input file cts 1.")
parser <- add_option(parser, c( "-j", "--input2" ),  action = "store", 
                     type="character", dest = "input2", 
                     help = "Input file cts 2.")
parser <- add_option(parser, c( "-o", "--output" ), action = "store", 
                     type="character", dest = "pca_filename", 
                     help = "Filename of eps file for PCA plot.")
options = parse_args(parser, positional_arguments = TRUE)

cts1 = as.matrix(read.csv(args[1], header=TRUE, sep=",", row.names = "X"))

cts2 = as.matrix(read.csv(args[2], header=TRUE, sep=",", row.names = "X"))

coldata1 = read.csv(gsub("cts","coldata",args[1]), 
                    header=TRUE, sep=",", row.names = "X")

coldata2 = read.csv(gsub("cts","coldata",args[2]), 
                    header=TRUE, sep=",", row.names = "X")

suppressPackageStartupMessages(library(DESeq2))

dds1 = DESeqDataSetFromMatrix(countData = cts1,
                             colData = coldata1,
                             design = ~ condition)

dds2 = DESeqDataSetFromMatrix(countData = cts2,
                             colData = coldata2,
                             design = ~ condition)

keep = rowSums(counts(dds1)) >= 10
dds1 = dds1[keep,]

keep = rowSums(counts(dds2)) >= 10
dds2 = dds2[keep,]

dds1$condition = factor(dds1$condition, levels = c("normal","cancerous"))
dds1 = DESeq(dds1)

dds2$condition = factor(dds2$condition, levels = c("normal","cancerous"))
dds2 = DESeq(dds2)

dds1 = dds1[row.names(dds1) %in% rownames(dds2), ]
dds2 = dds2[row.names(dds2) %in% rownames(dds1), ]

vsd1 <- varianceStabilizingTransformation(dds1, blind=FALSE)
vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)

pcaData1 <- plotPCA(vsd1, intgroup=c("condition"), returnData=TRUE)
pcaData2 <- plotPCA(vsd2, intgroup=c("condition"), returnData=TRUE)

setEPS()
postscript(args[3], width = 9, height = 6)

plot(pcaData1$PC1, pcaData1$PC2, col=factor(pcaData1$condition), pch=16, xlab="PC1", ylab="PC2", main="PCA plot", cex=1.5)
points(pcaData2$PC1, pcaData2$PC2, col=factor(pcaData2$condition), pch=17, cex=1.5)
#text(pcaData1$PC1, pcaData1$PC2, lab=pcaData1$name, cex=0.4)
#text(pcaData2$PC1, pcaData2$PC2, lab=pcaData2$name, cex=0.4)
legend(1.3, 0.23, legend=c("Collibri-cancerous", "Collibri-normal", "KAPA-cancerous", "KAPA-normal"),
       col=c("red", "black", "red", "black"), pch=c(16,16,17,17), cex=1)

dev.off()

