#### creation final table for gene phenotype and syndrome

## 1. load syndromic diseases data from different source ####
##1.a. data from monarchinitiative
setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/phenotypes_genotypes_files/')

###list of syndromic diseases with gene label
syndromic_diseases<- read.delim('syndromic_disease_monarch.tsv',header=T,sep='\t',fill=T,comment.char = '')

## remove duplicated to have exact number of gene and diseases
SD.nodup_syndromes<- syndromic_diseases[!duplicated(syndromic_diseases$object_label),c('subject_label','object_label','source','object')] ##1590 syndrome
SD.nodup_gene<- syndromic_diseases[!duplicated(syndromic_diseases$subject_label),c('subject_label','object_label','source','object')] ### 2011 phenotype
### gene with HPO terms 
## load a diseases phenotype files ###
diseases_phenotypes<-read.delim('disease_phenotype.all.tsv',header=T,sep='\t',fill=T,comment.char = '')
### select only phenotypes associated with syndrome 
phenotypes_syndromes<- diseases_phenotypes[grep('syndrome',diseases_phenotypes$subject_label),]
### remove duplicated
phenotypes_syndromes.nodup<-phenotypes_syndromes[!duplicated(phenotypes_syndromes$subject_label),] ## missing 211 syndromes

## list of SD that have HPO terms
syndrome1<-data.frame(mondo=phenotypes_syndromes.nodup$subject,syndrome=phenotypes_syndromes.nodup$subject_label,hpolabel='yes')
syndrome2<-data.frame(gene=SD.nodup_syndromes$subject_label,syndrome=SD.nodup_syndromes$object_label,mondo=SD.nodup_syndromes$object)
### merge 
all_syndrome<- merge(syndrome1,syndrome2, all=T)
missing<-all_syndrome[is.na(all_syndrome$hpolabel) | is.na(all_syndrome$gene),] ## list of syndrome without hpo label

##select only gene with know syndrome
SD_gene_pheno<-all_syndrome[!is.na(all_syndrome$hpolabel) & !is.na(all_syndrome$gene),] ## list of syndrome without hpo label

### extract from the data syndrome that have a due to in the names in order to capture all thier phenotype ##
otherSD<-all_syndrome[grep('due',all_syndrome$syndrome) ,]
otherSD$subject_label<-as.vector(sapply(strsplit(as.character(otherSD$syndrome),'due'),'[',1))
otherSD<-otherSD[!duplicated(otherSD$subject_label),]
##check those SD in the phenotype files
As<-phenotypes_syndromes.nodup[grep('Alagille syndrome',phenotypes_syndromes.nodup$subject_label),] ## present in different names

AS<-phenotypes_syndromes.nodup[grep('Alagille syndrome|Angelman syndrome',phenotypes_syndromes.nodup$subject_label),] ## present in different names

## syndrom that may be missing in the gene_syndrome data
otherSD.3<- syndrome1[grep(paste0(as.vector(otherSD$subject_label),collapse = '|'),syndrome1$syndrome),] ## present in different names

SD_gene_phenov2<-merge(SD_gene_pheno,otherSD.3,all=T)

### merge to have the final data 
SD_gene_phenov2$syndrome_names<- sapply(strsplit(as.character(SD_gene_phenov2$syndrome),'due'),'[',1)

## exact number of syndrome 1428 syndromes
SD_gene_phenov3<-SD_gene_phenov2[!duplicated(SD_gene_phenov2$syndrome_names),]

## file with entire list of gene and syndrome
syndrome_gene_file<-data.frame(syndrome=syndromic_diseases$object_label,gene_related=syndromic_diseases$subject_label,mondo=syndromic_diseases$object)

all_data<-merge(SD_gene_phenov2,syndrome_gene_file,all=T)
all_data$gene_related<-ifelse(is.na(all_data$gene_related),'',as.character(all_data$gene_related))

### delete syndrome without gene
### syndrome without associated gene
cc<-new.Data[all_data$gene_related=='',]
### delete syndrome without related gene

all_data<-all_data[all_data$gene_related!='',]
### delete a syndrome with unidentified name
all_data<-all_data[!is.na(all_data$syndrome_names),]

#### create the list of syndrome with corresponding gene
all_data$group<- as.factor(all_data$syndrome_names) ## attribute a number id to each syndrome as a group
group_syndrome<-function(x) paste(x,collapse = ',') # function to aggregare value of interest 
new.Data<-all_data[!duplicated(all_data$group),]
new.Data$gene_synd<-aggregate(gene_related~group,all_data,group_syndrome)


### creation of syndrome phenotype list ###

pheno_list<- phenotypes_syndromes[grep(paste0(as.vector(new.Data$syndrome_names),collapse = '|'),phenotypes_syndromes$subject_label),]

#### create the list of syndrome with corresponding gene

synd<-pheno_list[!duplicated(pheno_list$subject_label),]
syndrome<-rle(as.character(pheno_list$subject_label))[[1]] ## create a group id number for duplicated variable
pheno_list$group<- as.factor(pheno_list$subject_label) ## attribute a number id to each syndrome(group)
group_syndrome<-function(x) paste(x,collapse = ',') # function to aggregare value of interest 
new.Data_pheno<-pheno_list[!duplicated(pheno_list$group),]
new.Data_pheno$hpo_terms<-aggregate(object~group,pheno_list,group_syndrome)

new.Data$subject_label<-new.Data$syndrome_name

##merge geno and geno
pheno_geno<- merge(new.Data,new.Data_pheno)

####listof gene with gene and associated phenotype
gene_syndrome.2<- data.frame(syndrome=gene_syndrome_final$syndrome,gene_related=gene_syndrome_final$gene_synd$gene_related,
                             hpo_terms=gene_syndrome_final$hpo_terms$object)

### list of all hpo terms for syndromic disease
hpo_list_SD<-pheno_list[!duplicated(pheno_list$object),c('object','object_label')]

gene_SD<-data.frame(gene_list=SD.nodup_gene[,'subject_label'])

### save all the file
##1. syndrome list file
write.table(gene_syndrome.2,'syndromic_diseases_list.txt',quote=F,sep='\t',row.names=F,col.names=T)
write.table(hpo_list_SD,'hpo_list_SD.txt',quote=F,sep='\t',row.names=F,col.names=T)
write.table(gene_SD,'gene_list_SD.txt',quote=F,sep='\t',row.names=F,col.names=T)





