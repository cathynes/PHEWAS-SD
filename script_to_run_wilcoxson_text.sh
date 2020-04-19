### gene level association using EPACTS 


## 1. Annotate my VCF files 

cd /home/storage/catherine_storage/subset.ukbb.for.phewas
###
/home/tcheandj/EPACTS/bin/epacts anno --in all.genes.subset.vcf.gz --out all.genes.subset.anno.vcf.gz 


### create a marker group id

/home/tcheandj/EPACTS/bin/epacts make-group --vcf all.genes.subset.anno.vcf.gz  --out gene_marker_id --format [epacts, annovar, chaos or gatk] --nonsyn


### create ped style phenotype for gene level association 
library(data.table)
library(dplyr)

pheno<-fread("/home/storage/catherine_storage/subset.ukbb.for.phewas/phewas_phewas_SD_hpo_and_biomarkers.pheno")
covar<- fread("/home/storage/catherine_storage/subset.ukbb.for.phewas/phewas_phewas_SD_hpo_and_biomarkers.covar")

str(covar)
all<- merge(covar,pheno,sort = F)

str(all)

all<- all %>% filter(FID>1) %>% mutate(`#FAM_ID`=FID,IND_ID=paste0(IID,"_",IID),FAT_ID=0,MOT_ID=0,SEX=sex.l2) %>% select(`#FAM_ID`,IND_ID,FAT_ID,MOT_ID,SEX,age:PC5,Alanine_aminotransferase_adjusted:HP0002107)

## saVE THE PHENOTYPE
write.table(all,"/home/storage/catherine_storage/subset.ukbb.for.phewas/pheno_phewas_SD_hpo_and_biomarkers.ped",sep='\t',na="NA",quote=F,col.names = T,row.names = F)

### test a single phenotype

### create a  single ped files for the phenotypes and covariate and loop throught for the analysis 

/home/tcheandj/EPACTS/bin/epacts group --vcf all.genes.subset.anno.vcf.gz  \
  --groupf SD.gene.markers.id  --out skat_SD_test \
  --ped pheno_phewas_SD_hpo_and_biomarkers.ped  --max-maf 0.05 \
  --pheno Alanine_aminotransferase_adjusted --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov age \
  --cov Array --cov sex.l2 --test skat --skat-o --run 2


/home/tcheandj/EPACTS/bin/epacts group \
	--vcf /home/storage/catherine_storage/subset.ukbb.for.phewas/all.genes.subset.anno.vcf.gz  \
  --groupf /home/storage/catherine_storage/subset.ukbb.for.phewas/SD.gene.markers.id  \
  --out /home/tcheandj/phewas_SD_revision_plos_genetics/result_wilcox.Alanine_aminotransferase_adjusted \
  --ped /home/storage/catherine_storage/subset.ukbb.for.phewas/pheno_phewas_SD_hpo_and_biomarkers.ped  \
  --pheno Alanine_aminotransferase_adjusted --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov age \
  --max-maf 0.5  --test q.wilcox  --run 4  

/home/tcheandj/EPACTS/bin/epacts group \
	--vcf /home/storage/catherine_storage/subset.ukbb.for.phewas/all.genes.subset.anno.vcf.gz  \
  --groupf /home/storage/catherine_storage/subset.ukbb.for.phewas/SD.gene.markers.id  \
  --out /home/tcheandj/phewas_SD_revision_plos_genetics/result_wilcox.Alanine_aminotransferase_adjusted \
  --ped /home/storage/catherine_storage/subset.ukbb.for.phewas/pheno_phewas_SD_hpo_and_biomarkers.ped  \
  --pheno Alanine_aminotransferase_adjusted --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov age \
  --max-maf 0.5  --test emmaxVT  --run 4 


  ### create a looping for  the  analysis
  ## wilcoxon on CC phenpotypes

for pheno in HP0000002	HP0000821	HP0005117	HP0001507	HP0003758	HP0045081	HP0001518	HP0003124	HP0004421	Albumin_adjusted	HP0012603	HP0007018	HP0000077	HP0000120; do

	echo "!/bin/bash

	### wilcoxson test 

	/home/tcheandj/EPACTS/bin/epacts group \
	--vcf /home/storage/catherine_storage/subset.ukbb.for.phewas/all.genes.subset.anno.vcf.gz  \
  --groupf /home/storage/catherine_storage/subset.ukbb.for.phewas/SD.gene.markers.id  \
  --out /home/tcheandj/phewas_SD_revision_plos_genetics/result_wilcox.$pheno \
  --ped /home/storage/catherine_storage/subset.ukbb.for.phewas/pheno_phewas_SD_hpo_and_biomarkers.ped  \
  --max-maf 0.5 --pheno $pheno  --cov PC1 --cov PC2 --cov PC3 --cov PC4 --cov PC5 --cov age \
  --test q.wilcox --run 3 " > wilcoxson.$pheno.sh ;

  chmod u=rwx wilcoxson.$pheno.sh;

  ## sudmit 
  nohub ./wilcoxson.$pheno.sh & ;
  done

  #### 


