## loop for reading gene file and write a SKAT script for each gene  #####
cat ../genes.list.txt | while read gene
do

echo $gene

echo "#!/bin/bash/R
library(data.table)
library(dplyr)
library(SKAT)


setwd('/oak/stanford/groups/jpriest/catherine/phewas_SD_revision_plos_genetic/phewas_SD_hpo_and_biomarkers')

pheno_SD<- fread('phenotype.phewas.revision_plosgen.txt')

str(pheno_SD)
## change ID in pheno data
pheno_SD\$IID<-pheno_SD\$app13721
pheno_SD\$app13721<-NULL 


str(pheno_SD)

#pheno.name<-names(pheno_SD)
###genotyping data

$gene<-fread('/scratch/users/tcheandj/revision_phewas_plos_genetic/data/dosage.$gene.raw')

str($gene)

geno.name<-names($gene)

### before running SKAT, create a new data with the exact same number of Ind as in the genotype

#### merge data to have the same number of ind

all.data<-merge(pheno_SD, $gene,by='IID')

str(all.data)
## split for each process

pheno<-all.data[,c(1:79)]
str(pheno)

geno<-all.data %>% select(geno.name)
str(geno)

### check dimension of each data 
dim(pheno)
dim(geno)

### get binary phenotype and contine phenotype
pheno_desc<- fread('pheno_desc.phewas.revision_plosgen.txt')

bin.pheno <- pheno_desc %>% filter(type=='D') %>% select(name)


### loop the skat test binary trait ####

pvalue<- data.frame()

for ( i in 39:93){


  outcome = colnames(pheno)[i]
 
	null_model_odj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='D')

	SKAT.test<- SKAT_CommonRare(as.matrix(geno), 
               null_model_odj,
               is_dosage=TRUE,
               weights.beta.rare=c(1, 25),
               weights.beta.common=c(0.5, 0.5),
               CommonRare_Cutoff=NULL,
               missing_cutoff=0.05,
               estimate_MAF=1)

	pvalue[paste0(outcome),'n_markers.test']<- SKAT.test\$param\$n.marker.test
	pvalue[paste0(outcome),'n_common']<- SKAT.test\$n.common
	pvalue[paste0(outcome),'n_rare']<- SKAT.test\$n.rare
	pvalue[paste0(outcome),'pval_RaCo']<-SKAT.test\$p.value

	nul.odj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='D')

	pvalue[paste0(outcome),'pval_raco_burden'] <- SKAT_CommonRare(as.matrix(geno), nul.odj, r.corr.rare=1, r.corr.common=1)

	nul.odj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='D')

	pvalue[paste0(outcome),'pval_Ra']<- SKAT_CommonRare(as.matrix(geno),nul.odj, test.type='Rare.Only')
}

# get continue phenotype
cont.pheno <- pheno_desc %>% filter(type=='C') %>% select(name)

#### loop skat for continue trait

for ( i in 10:38){

  outcome = colnames(pheno)[i]
  
	null_model_obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='C')

	SKAT.test<- SKAT_CommonRare(as.matrix(geno), null_model_obj, method='AR',
               is_dosage=TRUE, weights.beta.rare=c(1, 25),
               weights.beta.common=c(0.5, 0.5), CommonRare_Cutoff=NULL,
               missing_cutoff=0.05, estimate_MAF=1)

	pvalue[paste0(outcome),'n_markers.test']<- SKAT.test\$param\$n.marker.test
	pvalue[paste0(outcome),'n_common']<- SKAT.test\$n.common
	pvalue[paste0(outcome),'n_rare']<- SKAT.test\$n.rare
	pvalue[paste0(outcome),'pval_RaCo']<-SKAT.test\$p.value

	nul.obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='C')

	pvalue[paste0(outcome),'pval_raco_burden'] <- SKAT_CommonRare(as.matrix(geno), nul.obj, r.corr.rare=1, r.corr.common=1)

	#null_model_obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='C') 

	pvalue[paste0(outcome),'pval_Ra']<- SKAT_CommonRare(as.matrix(geno), nul.obj,test.type='Rare.Only')
}


pvalue\$HPO<- rownames(pvalue)

##write output 
write.table(pvalue,'/scratch/users/tcheandj/revision_phewas_plos_genetic/skat_test_results/pvalueskat_$gene.txt',row.names=F)
" > skat.test.$gene.R 
echo skat.text.$gene.R

echo skat.text.$gene.R


 sbatch -p jpriest,normal,owners  --qos=normal -N 1 -n 1 -t 05:00:00 --mem=100000  -J $gene.skat -o log/$gene.skat.log -e log/$gene.skat.err  R CMD BATCH  skat.test.$gene.R 
done
