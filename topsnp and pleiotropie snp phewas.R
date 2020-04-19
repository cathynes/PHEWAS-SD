library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(metafor)
library(rmeta)

#### extracting info from clinvar data ####

single_snp<-fread('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/phewas_SD.final.results.txt')

##### table of the top signal by phenotype

top.signal<- single_snp %>% filter(P<3.16e-07) %>% group_by(outcome) %>% filter(P==min(P)) %>% ungroup()

snp.top.signal <- c('rs11066309',  'rs11917587',   'rs2035936',   'rs3821710', 'rs589668',   'rs6040076',  
                    'rs9852128','rs143997339')

top.signal<- single_snp %>% filter(P<5e-05, rsid%in%snp.top.signal) 
fwrite(top.signal,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/result.top.signal.txt',sep='\t')

#### table of variant that reach genome wide significance for at least 2 phenotype per syndrom
## snp with pleitropie
pleiotropie.snp<- single_snp %>% filter(P<3.2e-07) %>%  select(ID,rsid,locus_name,outcome,BETA,SE,OR,IC,P) %>% filter(duplicated(ID))

snp.with.pleio<-pleiotropie.snp[!duplicated(pleiotropie.snp$ID),'ID']

pleiotropie.snp<- single_snp %>% filter(ID%in%snp.with.pleio,P<1e-04)
pleiotropie.snp2<- pleiotropie.snp %>% filter(!duplicated(ID))

#### write the results 
fwrite(pleiotropie.snp,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/pleiotropie.snp.txt',sep='\t')
fwrite(pleiotropie.snp2,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/pleiotropie.withoutdup.txt',sep='\t')


####  extraction of result by syndrome
gene.result<-fread('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/phewas.SD.co.ra.variants.all_loci.all_phenotypes.skat.txt.meta_data.txt')

gene.result2<- gene.result %>% group_by(locus_name) %>% filter(weight%in%c("maf","weight_cadd"),set=='all') %>% 
  mutate(pheno_text=ifelse(outcome=='HP0000146','ovarian cyst',pheno_text)) %>% 
  mutate(pheno_text=ifelse(outcome=='HP0000486','Strabism',pheno_text)) %>% 
  mutate(cat_text=ifelse(outcome=='HP0000821','DIGEORGE syndrome, Noonan syndrome',cat_text)) %>% 
  mutate(cat_text=ifelse(cat_text=='','DIGEORGE syndrome, Noonan syndrome, Marfan syndrome, Alagille syndrome',cat_text)) %>% 
  mutate(cat_text=ifelse(outcome=='HP0000486','Alagille syndrome',cat_text)) %>%
  mutate(label=paste0(pheno_text,' (',outcome,')','; ',locus_name,';',' p=', formatC(p_value, format = "e", digits = 2))) %>%
  mutate(pheno_text=ifelse(outcome=='HP0000077','abnormality of the kidney',pheno_text)) %>% 
  select(outcome,pheno_text,locus_name,set,weight,n1,n2,n_marker_test,p_value,cat_text,label)

#### value for weight CADD
gene.result.cadd <- gene.result2 %>% filter(weight=='weight_cadd',set=='all') %>% 
  rename(P.weightCADD='p_value',label.cadd='label') %>%  select(-weight)

gene.result.maf <- gene.result2 %>% filter(weight=='maf',set=='all') %>% 
  rename(P.MAF='p_value',label.maf='label') %>%  select(-weight)

### merge both 
all.gene.result<-merge(gene.result.cadd,gene.result.maf,all=T)

#### extract result by syndrom for weight methods
AS.generesult<- gene.result.cadd %>%  filter(grepl("Alagille",cat_text), locus_name%in%c('JAG1','NOTCH2')) %>%  
                                      rename(pval='P.weightCADD',SNP='locus_name',phenotype='label') %>% 
                                      select(SNP,phenotype,pval)

MS.generesult<- gene.result.cadd %>%  filter(grepl("Marfan",cat_text), locus_name%in%c('FBN1')) %>%    
                                      rename(pval='P.weightCADD',SNP='locus_name',phenotype='label') %>% 
                                      select(SNP,phenotype,pval)

NS.generesultptpn11<- gene.result.cadd %>%  filter(grepl("Marfan",cat_text), locus_name%in%c('PTPN11')) %>%    
  rename(pval='P.weightCADD',SNP='locus_name',phenotype='label') %>% 
  select(SNP,phenotype,pval)

NS.gene<-c('SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS','PTPN11',
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B')

NS.generesult<- gene.result.cadd %>%  filter(grepl("Noonan",cat_text), locus_name%in%NS.gene) %>%   mutate(SNP='RASopathie') 
                               

gene<-c('PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B','FBN1',"JAG1",'NOTCH2')
DS.generesult<- gene.result.cadd %>%  filter(grepl("DIGEORGE",cat_text), !locus_name%in%gene, !is.na(P.weightCADD)) %>%    
                                     mutate(SNP='22q11') 

#### save all result
fwrite(all.gene.result,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/all.gene.result.txt',sep='\t')

fwrite(AS.generesult,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/AS.generesult.txt',sep='\t')
fwrite(MS.generesult,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/MS.generesult.txt',sep='\t')
fwrite(DS.generesult,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/DS.generesult.txt',sep='\t')
fwrite(NS.generesult,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/NS.generesult.txt',sep='\t')

fwrite(NS.generesultptpn11[,c('SNP','phenotype','pval')],'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/NS.generesultptpn11.txt',sep='\t')


### extraction of results by phenotype for fuma

## diastolic BP
diastolic.bp<- single_snp %>% filter(outcome=='HP0005117') %>% select(ID,`#CHROM`,POS,REF,ALT,ALT_FREQ,OBS_CT,BETA,SE,P)
fwrite(diastolic.bp,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/diastolic.bpresult.txt',sep='\t')

## systolic BP
systolic.bp<- single_snp %>% filter(outcome=='HP0004421') %>% rename(CHROM='#CHROM') %>%  select(ID,CHROM,POS,REF,ALT,ALT_FREQ,OBS_CT,BETA,SE,P)
fwrite(systolic.bp,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/systolic.bpresult.txt',sep='\t')

### abnormal body height 
body.height<- single_snp %>% filter(outcome=='HP0000002') %>% rename(CHROM='#CHROM') %>%  select(ID,CHROM,POS,REF,ALT,ALT_FREQ,OBS_CT,BETA,SE,P)
fwrite(body.height,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/body.height.txt',sep='\t')

### abnormal body height 
body.height<- single_snp %>% filter(outcome=='HP0000002') %>% rename(CHROM='#CHROM') %>%  select(ID,CHROM,POS,REF,ALT,ALT_FREQ,OBS_CT,BETA,SE,P)
fwrite(body.height,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/body.height.txt',sep='\t')

### cutanous adipocyte tissues 
cutanous.adip.tissues<- single_snp %>% filter(outcome=='HP0003758') %>% rename(CHROM='#CHROM') %>%  select(ID,CHROM,POS,REF,ALT,ALT_FREQ,OBS_CT,BETA,SE,P)
fwrite(cutanous.adip.tissues,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/cutanous.adip.tissues.txt',sep='\t')

### list of significant snp at 5x10e-07
snp.lower.p<- single_snp %>% filter(P<5e-07) %>%  select(ID,rsid,locus_name) %>% filter(!duplicated(ID))
 
table(snp.lower.p$locus_name) 
write.table(diastolic.bp[,'ID'],'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/all.snp.in.phewas.txt',
            sep='\t',quote=F,col.names=F,row.names=F)



