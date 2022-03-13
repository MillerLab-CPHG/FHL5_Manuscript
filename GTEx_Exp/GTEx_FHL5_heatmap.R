#GTEx_FHL5_heatmap.R
rm(list=ls())
pkgs<-c('tidyverse', 'scales', 'RColorBrewer','biomaRt' )
lapply(pkgs, library, character.only=T)


#get TSS to order genes by distance to rs948
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attributes = listAttributes(ensembl, page = "structure")

tss <- getBM(attributes = c("transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             filters = "external_gene_name", 
             values = fhl5_locus,
             mart = ensembl)

gene_order<-tss%>%
  group_by(external_gene_name)%>%
  dplyr::slice(which.min(transcription_start_site))%>%
  mutate(dist=abs(96612248-transcription_start_site))%>%
  dplyr::select(c(external_gene_name, dist))%>%
  arrange(dist)%>%
  pull(external_gene_name)

gene_order<-factor(gene_order, levels=gene_order)

dat<-read_delim('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.txt', '\t')

dat_fhl5locus<-dat%>%
  dplyr::filter(Description %in% gene_order)%>%
  gather(key=key, value=value, 
         `Adipose - Subcutaneous`:`Whole Blood`)%>%
  mutate(`log2(TPM)`=log2(value))

dat_fhl5locus$Description<-factor(dat_fhl5locus$Description, 
                                  levels=rev(gene_order))

dat_fhl5locus_artery<-dat_fhl5locus%>%
  separate(key, sep=" - ", into=c("Tissue_Group", "Tissue"))%>%
  dplyr::filter(grepl("Artery", Tissue_Group))

dat_fhl5locus_rest<-dat_fhl5locus%>%
  separate(key, sep=" - ", into=c("Tissue_Group", "Tissue"))%>%
  group_by(Tissue_Group)%>%
  arrange(-value)%>%
  dplyr::slice(1)

dat_fhl5locus_all<-rbind(dat_fhl5locus_rest,dat_fhl5locus_artery)

dat_fhl5locus_all$key<-paste0(dat_fhl5locus_all$Tissue_Group, "-", dat_fhl5locus_all$Tissue)  

library(RColorBrewer)

p<-ggplot(dat_fhl5locus_all,aes(key,Description,fill=`value`))+
  geom_tile(color= "white",size=0.1) + 
  scale_fill_gradient2()


gene_order

p <-p + theme(legend.position = "top")+
  theme_classic()+
  theme(plot.title=element_text(size = 14))+
  theme(axis.text.y=element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(plot.title=element_text(hjust=0))+
  theme(axis.text=element_text(size=10))+
  theme(legend.title=element_text(size=10))+
  theme(legend.text=element_text(size=12))+
  xlab("")+ylab("")

p

pdf("FHL5_locus_GTEx_exp.pdf", width=8, height=4)
p
dev.off()

