### description study pop phewas ###
library(dplyr)
library(data.table)
library(ggplot2)

## load data 
meta.data<-fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/meta.data.txt')
des.bin<-fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS/hpo_binary_vars.txt')
des.cont<-fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS/hpo_continue_vars.txt')

### merge files
meta.hpo.bin<-merge(meta.data,des.bin,all=T)
str(meta.hpo.bin)

meta.hpo.bin.finale<- meta.hpo.bin %>% filter(!is.na(N.cases)) %>% 
  mutate(freq_cases=paste0(N.cases,' (',round(100*N.cases/N,2),')'), 
        freq_controls=paste0(N.controls,' (',round(100*N.controls/N,2),')')) %>%
  select(name,pheno_text,freq_cases,freq_controls,cat_text,cui_id,ukb_id)


meta.hpo.cont<-merge(meta.data,des.cont)
str(meta.hpo.cont)

meta.hpo.cont.finale<- select(meta.hpo.cont,name,pheno_text,Mean:N,cat_text,cui_id,ukb_id)

### save data 
fwrite(meta.hpo.bin.finale,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/desc_bin_vars.txt')
fwrite(meta.hpo.cont.finale,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/desc_cont_vars.txt')

### merging rs id with snp id in the ukbb for CHD gws from aldo
single.snp.results<- fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/single.snp.metadata.txt')


rsid<- fread('/Users/catherine.t/Documents/diabetes_chd/GWAS_summarystat/rsid.select_chr.txt')

str(rsid)

colnames(rsid)<- c('ID', 'rsid', 'chromosome', 'position', 'alleleA', 'alleleB')

### merge rsid and smp result

#  extract snp list 
snp.list<- unique(single.snp.results$ID)

rsid.2<- rsid %>% filter(ID%in%snp.list) %>% select(ID,rsid) 

snp.result<-merge(rsid.2,single.snp.results,all=T) 
snp.result<- snp.result %>% filter(!is.na(P),!is.na(BETA), !is.na(SE)) %>% select(ID:ALT_CTRL_CT,ukb_id:cat_text)
snp.result <- snp.result %>% group_by(outcome) %>% mutate(rsid=ifelse(is.na(rsid), ID,rsid)) %>% 
  mutate(label=ifelse(P==min(P) & min(P)< 10e-8,paste0(rsid,'(',locus_name,")"),"")) %>% ungroup()

result.snp.forplot<- snp.result %>% mutate(hpo_label=paste0(pheno_text,' (',outcome,')',', p=',P)) %>% select(rsid,locus_name,outcome,P,hpo_label)

fwrite(snp.result,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/significant.snp.result.withrsid.txt')

fwrite(result.snp.forplot,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/result.for.plot.txt',sep='\t')
