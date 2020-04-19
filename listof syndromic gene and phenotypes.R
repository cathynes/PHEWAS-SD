#### list of gene and phenotypes ####
setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/phenotypes_genotypes_files/')

###list of syndromic diseases with gene label
syndromic_diseases<- read.delim('syndromic_disease_monarch.tsv',header=T,sep='\t',fill=T,comment.char = '')

## remove duplicated to have exact number of gene and diseases
SD.nodup_syndromes<- syndromic_diseases[!duplicated(syndromic_diseases$object_label),c('subject_label','object_label','source')] ##1590 syndrome
SD.nodup_gene<- syndromic_diseases[!duplicated(syndromic_diseases$subject_label),c('subject_label','object_label','source')] ### 2011 phenotype

## load a diseases phenotype files ###
diseases_phenotypes<-read.delim('disease_phenotype.all.tsv',header=T,sep='\t',fill=T,comment.char = '')
### select only phenotypes associated with syndrome 
phenotypes_syndromes<- subset(diseases_phenotypes,diseases_phenotypes$subject_label%in%SD.nodup_syndromes$object_label)

### chek if I have all the syndrome
phenotypes_syndromes.nodup<-phenotypes_syndromes[!duplicated(phenotypes_syndromes$subject_label),] ## missing 211 syndromes

## list of SD that have HPO terms
syndrome1<-data.frame(mondo=phenotypes_syndromes.nodup$subject,syndrome=phenotypes_syndromes.nodup$subject_label,hpolabel='yes')
syndrome2<-data.frame(gene=SD.nodup_syndromes$subject_label,syndrome=SD.nodup_syndromes$object_label)
### merge 
all_syndrome<- merge(syndrome1,syndrome2, all=T)

missing<-all_syndrome[is.na(all_syndrome$mondo),] ## list of syndrome without hpo label

###check syndrome names with expression due to 

SAwithdueto<- diseases_phenotypes[grep('due to',diseases_phenotypes$subject_label),]


AS<- diseases_phenotypes[grep('Alagille',diseases_phenotypes$subject_label),]

### checking of some syndrome that might not have hpo terms
Cst1<- diseases_phenotypes[grep('Cockayne syndrome',diseases_phenotypes$subject_label),]
ac<- diseases_phenotypes[grep('acute coronary syndrome',diseases_phenotypes$subject_label),] ## is not a syndromic diseases but just a bunch of phenotypes


syndromic_diseases1<- syndromic[!duplicated(syndromic$object_label),c('subject_label','object_label','source')] ##1590 syndrome
syndromic_diseases2<- syndromic[!duplicated(syndromic$subject_label),c('subject_label','object_label','source')] ### 2011 phenotype

#### create a table with a set of syndrome and phenotype associated with phenotype
syndromic_diseases1<- syndromic[duplicated(syndromic$object_label),c('subject_label','object_label','source')]

