#FHL5_locus_coloc.R
pkgs<-c('coloc', 'tidyverse', 'data.table')
lapply(pkgs, library, character.only=T)

setwd("/Users/doriswong/Documents")

eqtl_dir<-c("/Users/doriswong/Desktop/FHL5_Manuscript/Finemapping/DAtasets/GTEx")
gwas_dir<-c('/Users/doriswong/Desktop/FHL5_Manuscript/Finemapping/DAtasets/GWAS')
gwas_meta<-fread('/Users/doriswong/Documents/GWAS_meta.txt')

gwas_files<-list.files(gwas_dir)

## gathering all case control gwas files 
cc_gwas<-gwas_files[c(1,2,5,7,10)]
cc_n<-c(317636,639221,452264,361141,297415)
case_prop<-c(0.03,0.10,0.27,0.03,0.12)
cc_gwas_files<-paste0(gwas_dir,"/", cc_gwas)
cc_gwas_dat_list<-lapply(cc_gwas_files, fread)
names(cc_gwas_dat_list)<-cc_gwas

#keep unique rsid, drop NA
filt_snps<-function(x){
  df2<-x%>%
    distinct(SNP, .keep_all=TRUE)%>%
    drop_na(SNP)%>%
    dplyr::filter(freq > 0 &&freq <1 )
  x<-df2
}

cc_gwas_dat_list<-lapply(cc_gwas_dat_list, filt_snps)

## formatting GTEx eQTL data
gtex_header<-c('gene_id', 'variant_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf',
               'pval_nom', 'slope', 'slope_se', 'rsid')

aor_eqtl<-fread(paste0(eqtl_dir,"/", "FHL5_locus_GTEx_AOR_eQTL.txt"))
colnames(aor_eqtl)<-gtex_header

gtex_tissues<-c("AOR", 'COR', 'TIB')
sample_size<-c(387,213,584)

eqtl_dat_filt_split<-split(aor_eqtl, f=aor_eqtl$gene_id)

filt_snps_eqtl<-function(x){
  df2<-x%>%
    distinct(rsid, .keep_all=TRUE)%>%
    drop_na()%>%
    dplyr::filter(maf > 0 &&maf <1 )
  x<-df2
}

eqtl_dat_filt_split<-lapply(eqtl_dat_filt_split, filt_snps_eqtl)
names(eqtl_dat_filt_split)

dat_all_cc=data.frame()

## running coloc, looping across all GWAS files for all gene
for (j in 1:length(cc_gwas_dat_list)){
  print (names(cc_gwas_dat_list)[j])
  
  for (t in 1:length(eqtl_dat_filt_split)){
    
    print(names(eqtl_dat_filt_split)[t])
    
    dat2<-coloc.signals(
      dataset1=list(snp=cc_gwas_dat_list[[j]]$SNP,
                    N=cc_n[j],type="cc", 
                    s=case_prop[j], 
                    MAF=cc_gwas_dat_list[[j]]$freq, beta=cc_gwas_dat_list[[j]]$b, 
                    varbeta=(cc_gwas_dat_list[[j]]$se^2)), 
      dataset2=list(snp=eqtl_dat_filt_split[[t]]$rsid, 
                    type="quant", 
                    MAF=eqtl_dat_filt_split[[t]]$maf,
                    N=387, 
                    beta=eqtl_dat_filt_split[[t]]$slope,
                    varbeta=(eqtl_dat_filt_split[[t]]$slope_se^2)),
      p1 = 1e-04,
      p2 = 1e-04,
      p12= 1e-05,
      maxhits = 3)
    
    print("coloc done") 
    
    df_3<-data.frame(dat2$summary)%>%
      mutate(Gene=eqtl_dat_filt_split[[t]]$gene_id[1])%>%
      mutate(file=cc_gwas[j])
    
    dat_all_cc<-rbind(dat_all_cc, df_3)
  
  }
}


## quant_gwas 
gwas_files<-list.files(gwas_dir)
quant_gwas<-gwas_files[c(3,4,8,9)]
quant_n<-c(26909,459777,459777,459777)
quant_gwas_files<-paste0(gwas_dir,"/", quant_gwas)
quant_gwas_dat_list<-lapply(quant_gwas_files, fread)
names(quant_gwas_dat_list)<-quant_gwas

filt_snps<-function(x){
  df2<-x%>%
    distinct(SNP, .keep_all=TRUE)%>%
    drop_na()%>%
    dplyr::filter(freq > 0 &freq <1 )
  x<-df2
}


quant_gwas_dat_list<-lapply(quant_gwas_dat_list, filt_snps)

dat_all_quant<-data.frame()

for (j in 1:length(quant_gwas_dat_list)){
  print (names(quant_gwas_dat_list)[j])
  
  for (t in 1:length(eqtl_dat_filt_split)){
    
    print(names(eqtl_dat_filt_split)[t])
    
    dat2<-coloc.signals(
      dataset1=list(snp=quant_gwas_dat_list[[j]]$SNP,
                    N=quant_n[j],
                    type="quant", 
                    MAF=quant_gwas_dat_list[[j]]$freq, 
                    beta=quant_gwas_dat_list[[j]]$b, 
                    varbeta=(quant_gwas_dat_list[[j]]$se^2)), 
      dataset2=list(snp=eqtl_dat_filt_split[[t]]$rsid, 
                    type="quant", 
                    MAF=eqtl_dat_filt_split[[t]]$maf,
                    N=387, 
                    beta=eqtl_dat_filt_split[[t]]$slope,
                    varbeta=(eqtl_dat_filt_split[[t]]$slope_se^2)),
      p1 = 1e-04,
      p2 = 1e-04,
      p12= 1e-05,
      maxhits = 3)
  
    df_3<-data.frame(dat2$summary)%>%
      mutate(Gene=eqtl_dat_filt_split[[t]]$gene_id[1])%>%
      mutate(file_1=paste0(quant_gwas[j]))
    
    dat_all_quant<-rbind(dat_all_quant, df_3)
    
  }
}

colnames(dat_all_quant)<-colnames(dat_all_cc)

dat_all<-rbind(dat_all_quant,dat_all_cc)
write_delim(dat_all, "FHL5_eQTL_AOR_GWAS_coloc.txt", '\t')



