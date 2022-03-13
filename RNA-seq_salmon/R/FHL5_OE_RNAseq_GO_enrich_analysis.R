## DEseq2_Final 

rm(list=ls())

library(tximport) 
library(tidyverse)
library(data.table)
library(AnnotationDbi)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(ggrepel)


## gathering files
dir<-"/Users/doriswong/Documents/RNA_seq_final_102221/"
setwd(dir)
dir_2<-paste0(dir, "salmon_outdir_2/")
files<-list.files(dir_2, pattern="BASAL")
quant_files<-paste0(dir_2, files, "/quant.sf")

txdb <- makeTxDbFromGFF("/Users/doriswong/Documents/gencode.v31.annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")

#importing salmon data with tximport
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
txi.salmon <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)
sampleTable <- data.frame(condition = factor(c(rep("FHL5", 3), rep("HA", 3), rep('NLS', 3))))
rownames(sampleTable) <- files

## running DESeq2
library(DESeq2)

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

vsd <- vst(dds, blind = TRUE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("HA_FHL5_PCA.pdf", width=4, height=3)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_classic()+scale_color_brewer(palette="Dark2")+
  theme(legend.position="top", 
        axis.text = element_text(size=12), 
        axis.title=element_text(size=14))
dev.off()

dds <- DESeq(dds)

res_NLS<-results(dds, contrast=c("condition", "NLS", "HA"))
res_FHL5<-results(dds, contrast=c("condition", "FHL5", "HA"))

## FHL5-NLS vs. HA analysis 
resOrdered_NLS <- res_NLS[order(res_NLS$padj), ]
HA_NLS_all_genes<-data.frame(resOrdered_NLS)
HA_NLS_all_genes<-data.frame(HA_NLS_all_genes)%>%
  rownames_to_column("gene_name")%>%
  separate(gene_name, sep='\\.', into=c('gene_1'))
HA_NLS_all_genes$symbol <- mapIds(org.Hs.eg.db, HA_NLS_all_genes$gene_1, column="SYMBOL", keytype="ENSEMBL")
HA_NLS_sig_only<-HA_NLS_all_genes%>%
  filter(padj < 0.05 )%>%
  filter(abs(log2FoldChange) > 0.6)%>%
  rownames_to_column("gene_name")%>% 
  drop_na(symbol)
write_delim(HA_NLS_sig_only, "HA_NLS_DE_sig_genes.txt", '\t')

## FHL5 vs. HA analysis 

HA_FHL5_all_genes<-data.frame(resOrdered_FHL5)
HA_FHL5_all_genes<-data.frame(HA_FHL5_all_genes)%>%
  rownames_to_column("gene_name")%>%
  separate(gene_name, sep='\\.', into=c('gene_1'))
HA_FHL5_all_genes$symbol <- mapIds(org.Hs.eg.db, HA_FHL5_all_genes$gene_1, column="SYMBOL", keytype="ENSEMBL")
HA_FHL5_sig_only<-HA_FHL5_all_genes%>%
  filter(padj < 0.05 )%>%
  filter(abs(log2FoldChange) > 0.6)%>%
  rownames_to_column("gene_name")%>% 
  drop_na(symbol)
write_delim(HA_FHL5_sig_only, "HA_FHL5_DE_sig_genes.txt", '\t')

## upset plot comparing up and downregulated genes in FHL5/NLS cells 

FHL5_upgenes<-HA_FHL5_sig_only[HA_FHL5_sig_only$log2FoldChange>0,]$symbol
FHL5_downgenes<-HA_FHL5_sig_only[HA_FHL5_sig_only$log2FoldChange<0,]$symbol

NLS_upgenes<-HA_NLS_sig_only[HA_NLS_sig_only$log2FoldChange>0,]$symbol
NLS_downgenes<-HA_NLS_sig_only[HA_NLS_sig_only$log2FoldChange<0,]$symbol

degenes<-list('FHL5_upgenes'=FHL5_upgenes,'FHL5_downgenes'=FHL5_downgenes, 'NLS_upgenes'=NLS_upgenes, 'NLS_downgenes'=NLS_downgenes )
  
pdf('FHL5-NLS_DEgenes_intersect.pdf', width=6, height=4)
upset(fromList(degenes), 
      order.by = "freq", text.scale = c(2, 2, 2, 2, 2, 2)
      )
dev.off()

## pathway enrichment 

library(clusterProfiler)

fhl5_BP_enrich<- enrichGO(gene=HA_FHL5_sig_only$symbol,
                          universe= HA_FHL5_all_genes$symbol,
                          keyType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)

fhl5_BP_enrich<-data.frame(fhl5_BP_enrich)%>%
  mutate(GO="BP")

fhl5_MF_enrich<- enrichGO(gene=HA_FHL5_sig_only$symbol,
                          universe= HA_FHL5_all_genes$symbol,
                          keyType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db,
                          ont = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)

fhl5_MF_enrich<-data.frame(fhl5_MF_enrich)%>%
  mutate(GO="MF")

fhl5_CC_enrich<- enrichGO(gene=HA_FHL5_sig_only$symbol,
                          universe= HA_FHL5_all_genes$symbol,
                          keyType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db,
                          ont = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)

fhl5_CC_enrich<-data.frame(fhl5_CC_enrich)%>%
  mutate(GO="CC")

FHL5_GO_enrich<-rbind(fhl5_BP_enrich, fhl5_MF_enrich, fhl5_CC_enrich)

write_delim(data.frame(FHL5_GO_enrich), "FHL5_DEgenes_GO_enrich.txt")

##

nls_BP_enrich<- enrichGO(gene=HA_NLS_sig_only$symbol,
                          universe= HA_NLS_all_genes$symbol,
                          keyType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
nls_BP_enrich<-data.frame(nls_BP_enrich)%>%
  mutate(GO="BP")

NLS_MF_enrich<- enrichGO(gene=HA_NLS_sig_only$symbol,
                          universe= HA_NLS_all_genes$symbol,
                          keyType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db,
                          ont = "MF",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
NLS_MF_enrich<-data.frame(NLS_MF_enrich)%>%
  mutate(GO="MF")

NLS_CC_enrich<- enrichGO(gene=HA_NLS_sig_only$symbol,
                          universe= HA_NLS_all_genes$symbol,
                          keyType = "SYMBOL", 
                          OrgDb = org.Hs.eg.db,
                          ont = "CC",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = FALSE)
NLS_CC_enrich<-data.frame(NLS_CC_enrich)%>%
  mutate(GO="CC")

NLS_GO_enrich<-rbind(nls_BP_enrich, NLS_MF_enrich, NLS_CC_enrich)

write_delim(data.frame(NLS_GO_enrich), "NLS_DEgenes_GO_enrich.txt")


