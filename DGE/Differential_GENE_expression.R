setwd("C:/Users/Aashi/Desktop/DGE")
#####DGE#####



library(DESeq2)
library(apeglm)
BiocManager::install("apeglm")

#loading count matrix
gse = read.csv("GSE92945_Fibroblast_RNAseq_counts_1.csv", header = TRUE, row.names = 1)
info = read.table("ind.txt", header = T,sep ='\t') 
#creating a data matrix  
dds = DESeqDataSetFromMatrix(round(gse), info, ~conditions)

#removing lowly expressed genes
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

#main deseq
ddsDE = DESeq(dds)

#export normalized read counts
normCOunts = counts(ddsDE, normalized = T)
write.csv(normCOunts,"normalised_data.csv")

#DESEQ results
res = results(ddsDE, alpha = 0.05)

#output DESeq result
resOrdered = res[order(res$padj),]
write.csv(resOrdered, "resOrdered1.csv")

resOrdered







######plotting DGE8#####
library(ggplot2)
library(pheatmap)
library(cowplot)


#plotting PCA
vsat = vst(ddsDE, blind = FALSE)
plotPCA(vsat, intgroup = "conditions")

plotDispEsts(ddsDE)



normCOunt = read.csv('normalised_data.csv', row.names = 1)

deSeqRes = read.csv('resOrdered.csv',row.names = 1)

deSeqRes$sig =ifelse(deSeqRes$padj <= 0.05, 'yes', 'no') 

#plotMa
deSeqRes = na.omit(deSeqRes)
maplot = ggplot(deSeqRes, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) + geom_point()

#volcano plot 
volcano = ggplot(deSeqRes, aes(x = log2FoldChange, y = -log10(padj), color = sig)) + geom_point()
volcano = volcano + ylim(0,15)

plot_grid(maplot,volcano)
#pheatmap
signi = subset(deSeqRes, padj <= 0.05)
allsig= merge(normCOunt, signi, by = 0)

sigcounts = allsig[,2:11]
row.names(sigcounts) = allsig$Row.names 

#plotting and beautifying
pheat = pheatmap(log2(sigcounts + 1), scale = 'row', show_rownames = F, treeheight_row = 0, treeheight_col = 0)

plot_grid(volcano,pheat[[4]])
sessionInfo()


#####GSEA#####
#naming the first column as it was unnamed 
library(dplyr)
library(tidyverse)

rest <- read_csv("resOrdered1.csv")
head(rest)
rest = rest %>%
  rename(gene = ...1 )
#creating a mapping table
library(org.Hs.eg.db)
ens2symbol = AnnotationDbi::select(org.Hs.eg.db,keys = rest$gene, 
                                   columns="SYMBOL",
                                   keytype="ENSEMBL")
ens2symbol <- as_tibble(ens2symbol)
view(ens2symbol)

rest <- inner_join(rest, ens2symbol, by=c("gene"="ENSEMBL"))
rest

res2 <- rest %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2

#GSEA using fgsea
library(fgsea)
ranks <- deframe(res2)
head(ranks, 20)

#loading reference enrichment database 
library("msigdbr")

Hs.GOBP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
Hs.Reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
Hs.Hallmark <- msigdbr(species = "Homo sapiens", category = "H")

Hs.GOBP.Entrez <- Hs.GOBP %>% split(x = .$entrez_gene, f = .$gs_name)
Hs.Hallmark.Entrez <- Hs.Hallmark %>% split(x = .$entrez_gene, f = .$gs_name)
Hs.GOBP.Symbol <- Hs.GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
Hs.Hallmark.Symbol <- Hs.Hallmark %>% split(x = .$gene_symbol, f = .$gs_name)




fgseaRes.Hallmark <- fgsea::fgseaMultilevel(pathways=Hs.Hallmark.Symbol,stats = ranks, eps =0)
# Make the Results tidy

###making it tidy 
fgseaResTidy <- fgseaRes.Hallmark %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

head(fgseaResTidy,10)

#plotting
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
