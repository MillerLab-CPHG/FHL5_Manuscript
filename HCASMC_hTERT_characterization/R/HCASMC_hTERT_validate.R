## comparing transcriptome of HCASMC-hTERT with HCASMC, ENCODE cell lines
rm(list=ls())

library(tximport) 
library(rhdf5)
library(tidyverse)
library(data.table)
library(AnnotationDbi)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(ggrepel)

dir<-"/Users/doriswong/Documents/ENCODE_RNAseq_salmon_out/"

setwd(dir)

files<-list.files(dir)

quant_files<-paste0(dir, files, "/quant.sf")

txdb <- makeTxDbFromGFF("/Users/doriswong/Documents/gencode.v31.annotation.gtf")

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

txi.salmon <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

sampleTable <- data.frame(celline=c('HUVECS', 'HUVECS',  'HUVECS', 
                          'K562', 'HEPG2', 'K562', 'HEPG2', 
                          'Cosmo', 'Cosmo', 'Cosmo', 'HCASMC', 'HCASMC',
                          'HCASMC'))

rownames(sampleTable)<-files


library(DESeq2)

dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~celline)

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

vsd <- vst(dds, blind = TRUE)

pcaData <- plotPCA(vsd, intgroup=c("celline"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("Cosmo_Validate_PCA.pdf", width=4, height=4)
ggplot(pcaData, aes(PC1, PC2, color=celline)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+theme_classic()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=14),
        legend.position = "top")+
  scale_color_brewer(palette="Dark2")
dev.off()

        
## HCASMC_htERT vs. HCASMC (2105)
rm(list=ls())

library(tximport) 
library(rhdf5)
library(tidyverse)
library(data.table)
library(AnnotationDbi)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(ggrepel)

dir<-"/Users/doriswong/Downloads/"

setwd(dir)

files<-list.files(dir, pattern="HCASMC")

quant_files<-paste0(dir, files, "/quant.sf")

txdb <- makeTxDbFromGFF("/Users/doriswong/Documents/gencode.v31.annotation.gtf")

k <- keys(txdb, keytype = "TXNAME")

### make sure you don't have tidyverse read in 
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")

txi.salmon <- tximport(quant_files, type = "salmon", tx2gene = tx2gene)

tpm<-txi.salmon$abundance

colnames(tpm)<-c('HCASMC_2105', 'HCASMC-hTERT_2105')

tpm_df<-data.frame(tpm)%>%
  rownames_to_column("ENID")%>%
  separate(ENID, sep="\\.", into=c('ENID'))

tpm_df$gene<- mapIds(org.Hs.eg.db, tpm_df$ENID, column="SYMBOL", keytype="ENSEMBL")

tpm_df%>%
  mutate(HCASMC_2105=log(HCASMC_2105, 2))%>%
  mutate(HCASMC.hTERT_2105=log(HCASMC.hTERT_2105, 2))%>%
  ggplot(aes(x=HCASMC_2105, y=HCASMC.hTERT_2105))+
  geom_point(alpha=0.2, color="black")+
  geom_smooth(method = "lm", col = "red") +
  theme_classic()

tpm_df%>%
  dplyr::filter(HCASMC_2105 !=0 )%>%
  dplyr::filter(HCASMC.hTERT_2105 != 0)%>%
  mutate(HCASMC_2105=log(HCASMC_2105, 2))%>%
  mutate(HCASMC.hTERT_2105=log(HCASMC.hTERT_2105, 2))%>%
  ggplot(aes(x=HCASMC_2105, y=HCASMC.hTERT_2105))+
  geom_point(alpha=0.2, color="black")+
  geom_smooth(method = "lm", col = "red", se = TRUE) +
  theme_classic()+xlab("log2(HCASMC_2105-TPM")+ylab("log2(HCASMC-hTERT-TPM)")+
  ggtitle("r=0.934")


cor.test(tpm_trans$HCASMC_2105, tpm_trans$HCASMC.hTERT_2105, method = "spearman", conf.level = 0.95)
