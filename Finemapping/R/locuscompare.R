#locuscompare plots to show FHL5 eQTL /MI GWAS colocalizations
pkgs<-c('tidyverse', 'locuscomparer')
lapply(pkgs, library,character.only=T)

dir="/Users/doriswong/Documents/"

fhl5_eqtl<-read_delim(paste0(dir,"FHL5_AOR_cis_rsid_p.txt"), 
                      '\t', col_names=F)
colnames(fhl5_eqtl)<-c('rsid', 'pval')

mi_rs948<-read_delim(paste0(dir, "FHL5_MI_EHJ_locuszoom.txt"), 
                     '\t', col_names=T)
colnames(mi_rs948)<-c('rsid', 'pval')


pdf("MI_FHL5_eQTL_locuscompare.pdf", width=6, height=3)
locuscompare(in_fn1 = mi_rs948, 
             in_fn2 = fhl5_eqtl, 
             title1 = 'MI GWAS', 
             title2 = 'FHL5 eQTL (AOR)')
dev.off()
