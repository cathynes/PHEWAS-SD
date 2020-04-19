#!/bin/bash
## SKAT test for all gene and all phenotype ###
module load R
## loop for running all the SKAT test in one #####
genelist = `cat /oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/gene_list.txt`

for gene in genelist 
do

echo $gene

echo "
library(data.table)
library(dplyr)
library(SKAT)


setwd('/oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants')

pheno_SD<- fread('phenotype.txt')

str(pheno_SD)
## change ID in pheno data
pheno_SD\$IID<-pheno_SD$app13721
pheno_SD\$app13721<-NULL 


str(pheno_SD)

pheno.name<-names(pheno_SD)
###genotyping data

$gene<-readRDS('/oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/extract/$gene/$gene.add.dosages.rds')

str($gene)

geno.name<-names($gene)

### before running SKAT, create a new data with the exact same number of Ind as in the genotype

#### merge data to have the same number of ind

all.data<-merge(pheno_SD, $gene,by='IID')

str(all.data)
## split for each process

pheno<-all.data %>% select(pheno.name)

geno<-all.data %>% select(geno.name)

### check dimension of each data 
dim(pheno)
dim(geno)

### get binary phenotype and contine phenotype
pheno_desc<- fread('/oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/pheno_desc.txt')

bin.pheno <- pheno_desc %>% filter(type=='D') %>% select(name)


### loop the skat test binary trait ####

pvalue<- data.frame()

for ( i in 1:nrow(bin.pheno)){

  outcome = bin.pheno[[i]]

  paste(bin.pheno[[i]])

obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='D')

SKAT.test<- SKAT(as.matrix(geno), obj, kernel = 'linear.weighted')

pvalue[paste0(outcome),'pval']<-SKAT.test\$p.value

}

# get continue phenotype
cont.pheno <- pheno_desc %>% filter(type=='C') %>% select(name)

#### loop skat for continue trait


#pvalue2<- data.frame()

#listhpo2<- as.list(names(pheno[,9:16]))


for ( i in 9:16){

  outcome = cont.phenotype[[i]]
  paste(outcome)

	obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='C')

	SKAT.test<- SKAT(as.matrix(geno), obj, kernel = 'linear.weighted')

pvalue[paste0(outcome),'pval']<-SKAT.test\$p.value

}


pvalue$HPO<- rownames(pvalue)

### merge pvalue 1 et 2

#pvalue_all<-merge(pvalue,pvalue2,all=T)

write.table(pvalue,'/oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/skat_test_results/pvalueskat_$gene.txt',row.names=F)
" > /oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/script_skattest/skat.test.$gene 
echo /oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/script_skattest/skat.text.$gene
#cd /oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/skat_test_log
#sbatch -p jpriest --qos=normal -N 1 -n 1 -t 05:00:00 --mem=100000  -J $gene.skat -o $gene.skat.log -e $gene.skat.err  R CMD BATCH  /oak/stanford/groups/jpriest/catherine/phewas.SD.co.ra.variants/script_skattest/skat.test.$gene 
done;

