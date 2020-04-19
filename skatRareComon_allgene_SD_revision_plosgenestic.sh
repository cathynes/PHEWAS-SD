#!/bin/bash

###1. create raw files with SNPs to be used in SKAT
while read f; do 
	echo $f
	grep $f snp.in.phewas.with.genes.annnot.txt | awk 'FS=OFS="\t" {print $8}' > $f.snp.list.txt ;  
done < genes.list.txt 


### creat the genes list with chr names 
awk 'FS=OFS="\t" {print $4,$7}' snp.in.phewas.with.genes.annnot.txt | sort | uniq > genes.lis_and_chrt.txt 

##### create a script to export each chr
#### script to extract gene resion for binbio
for i in $(seq 1 26); do 
  ### extract value of interest in each line of my table file
  var_line=$(sed -n ${i}p ../genes.lis_and_chrt.txt ); 
  chr=$(echo $var_line | cut -d' ' -f1);
  echo $chr;
  genes=$(echo $var_line | cut -d' ' -f2);
  echo $genes;

### write the corresponding script 

  echo "#!/bin/sh

  ### extraction of subset of SNP in a VCF file 
  /home/users/tcheandj/phewas/plink2 \
  --extract /scratch/users/tcheandj/revision_phewas_plos_genetic/$genes.snp.list.txt \
  --pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr${chr}_v3.mac1 \
  --remove /oak/stanford/groups/jpriest/catherine/phewas_SD_revision_plos_genetic/phewas_SD_hpo_and_biomarkers/pheno/phewas_phewas_SD_hpo_and_biomarkers.remove \
  --mach-r2-filter 0.6 2.0 \
  --export A \
  --out /scratch/users/tcheandj/revision_phewas_plos_genetic/data/dosage.$genes " > export_dosage.$genes.sh

chmod u= rwx export_dosage.$genes.sh
#sbatch -p jpriest,normal,owners --qos=normal -J plinK_export.$genes -t 24:00:00 -N 1 -n 1 -o log/plinK_export.$genes.log -e log/plinK_export.$genes.err --mem=22000 --open-mode=append export_dosage.$genes.sh

done

####. create the results data frame
mkdir -p /scratch/users/tcheandj/revision_phewas_plos_genetic/skat_test_results/

cd scratch/users/tcheandj/revision_phewas_plos_genetic/script/

###SKAT test for all gene and all phenotype ###
module load R/3.4.0
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

$gene<-fread('/scratch/users/tcheandj/revision_phewas_plos_genetic/data/$gene.add.dosages.rds')

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

for ( i in 25:79){

  outcome = colnames(pheno)[i]

	obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='D')

	SKAT.test<- SKAT_CommonRare(as.matrix(geno), qobj,
               null_model_obj,
               is_dosage=TRUE,
               method="C", test.type="Joint",
               weights.beta.rare=c(1, 25),
               weights.beta.common=c(0.5, 0.5),
               CommonRare_Cutoff=NULL,
               missing_cutoff=0.05,
               estimate_MAF=1
            )

pvalue[paste0(outcome),'pval']<-SKAT.test\$p.value

}

# get continue phenotype
cont.pheno <- pheno_desc %>% filter(type=='C') %>% select(name)

#### loop skat for continue trait

for ( i in 10:24){

  outcome = colnames(pheno)[i]
  
	obj<-SKAT_Null_Model(get(outcome) ~ age+sex+Array+PC1+PC2+PC3+PC4+PC5, data=pheno, out_type='C')

	SKAT.test<- SKAT_CommonRare(as.matrix(geno), qobj,
               null_model_obj,
               is_dosage=TRUE,
               method="C", test.type="Joint",
               weights.beta.rare=c(1, 25),
               weights.beta.common=c(0.5, 0.5),
               CommonRare_Cutoff=NULL,
               missing_cutoff=0.05,
               estimate_MAF=1
            )

pvalue[paste0(outcome),'pval']<-SKAT.test\$p.value

}


pvalue\$HPO<- rownames(pvalue)

##write output 
write.table(pvalue,'/oak/stanford/groups/jpriest/catherine/phewas_SD_revision_plos_genetic/phewas_SD_hpo_and_biomarkers/skat_test_results/pvalueskat_$gene.txt',row.names=F)
" > skat.test.$gene.R 
echo skat.text.$gene.R


# sbatch -p jpriest,normal,owners  --qos=normal -N 1 -n 1 -t 05:00:00 --mem=60000  -J $gene.skat -o log/$gene.skat.log -e log/$gene.skat.err  R CMD BATCH  skat.test.$gene.R 
done
