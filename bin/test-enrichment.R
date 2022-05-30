#!/usr/bin/env Rscript

# Sources:
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://stephenturner.github.io/deseq-to-fgsea/

suppressMessages(library(fgsea))
suppressMessages(library(ggplot2))
suppressMessages(library( "optparse" ))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(tidyft))

parser <- OptionParser()
parser <- add_option(parser, c( "-i", "--input" ),  action = "store", 
                     type="character", dest = "input", 
                     help = "Input file (res.csv or resLFC.csv).")
parser <- add_option(parser, c( "-o", "--output" ), action = "store", 
                     type="character", dest = "output", 
                     help = "Filename for output (eps format).")
options = parse_args(parser, positional_arguments = TRUE)
set.seed(42)

# Input - clean rownames
res = read.csv(options$options$input)
res = transform(res, row = as.character(row))
og_row_names = res$row
edited_row_names = gsub("\\.[0-9]+$", "",og_row_names)
res$row = edited_row_names
# Find entrezid ids based on ensembl ids


ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=res$row, 
                                    columns="ENTREZID",
                                    keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)

res <- inner_join(res, ens2symbol, by=c("row"="ENSEMBL"))
res2 = res

# Keep only two rows
res2=res2[c("ENTREZID","log2FoldChange")]

# Drom na, change form
res2 = res2 %>% drop_na()
ranks <- deframe(res2)

# analysis
pathways <- reactomePathways(names(ranks))
fgseaRes <- fgsea(pathways, ranks, maxSize=500, nperm=1000)
fgseaRes2 <- fgseaRes[order(NES,pval),]
#fgseaRes2 = as.data.frame(fgseaRes2)
#fgseaRes2$leadingEdge <- sapply(fgseaRes2$leadingEdge, paste, collapse=",")
#write.table(fgseaRes2,options$options$output, row.names = TRUE, sep="/t")


# visualization 1
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#fname = gsub("\\.tsv", ".eps", options$options$output)
setEPS()
postscript(options$options$output, width = 10, height = 8)
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)
dev.off()

# visualization 2 - cancer only
fname = gsub("\\.eps", "-cancer.eps", options$options$output)
setEPS()
postscript(fname, width = 10, height = 2)
cancer_related = grep("(C|c)ancer", names(pathways), value=TRUE)
plotGseaTable(pathways[cancer_related], ranks, fgseaRes, 
              gseaParam=0.5)
dev.off()

#fgseaResTidy <- fgseaRes %>%
#  as_tibble() %>%
#  arrange(desc(NES))
# Show in a nice table:
#fgseaResTidy %>% 
#  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#  arrange(padj) %>% 
#  DT::datatable()

# visualisation 3 - NES and p-values for each pathway
#fname2 = gsub("\\.eps", "-NES-pval.eps", options$options$output)
#setEPS()
#postscript(fname, width = 4, height = 10)
#ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#  geom_col(aes(fill=pval<0.1)) +
#  coord_flip() +
#  labs(x="Pathway", y="Normalized Enrichment Score",
#       title="Reactome pathways NES") + 
#  theme_minimal()
#dev.off()
