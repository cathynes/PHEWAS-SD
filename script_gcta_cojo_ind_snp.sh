## create gwas sum stat for linear regression 

for files in  results_chr*.glm.linear; do awk 'OFS=FS="\t" { print $3,$5,$4,$7,$10,$11,$12,$9}' $files > independent_snp_gcta/$files.gcta.txt ; done 

### gcta files for logistic regression 
for files in  results_chr*.glm.logistic.hybrid; do awk 'OFS=FS="\t" { print $3,$5,$4,$9,$12,$13,$14,$11}' $files > independent_snp_gcta/$files.gcta.txt ; done 

### remove .glm.linear and replace by .txt
cd 
rename rename glm.linear.gcta.txt txt *glm.linear.gcta.txt

rename glm.logistic.hybrid.gcta.txt txt *glm.logistic.hybrid.gcta.txt
### create sum stat for logistic regression


## gcta cojo

for i in $seq(1 11); do
	### run gcta cojo toidentify the set of ind snp in each genes
	var_line=$(sed -n$ {i}p ../phenotypes_with_sign_ass_for_gcta.txt);
	
	chr=$(echo $var_line | cut -d' ' -f1);
	
	echo $chr;
	
	pheno=$(echo $var_line | cut -d' ' -f2);
	
	echo $pheno;

	echo "#!/bin/sh

	/scratch/users/tcheandj/gcta_1.93.0beta/gcta64 \
	--bfile /scratch/users/tcheandj/subset.ukbb.for.clumping/ukbb.subset.for.prs.ld.chr${chr} \
	--chr $chr \
	--cojo-file /scratch/users/tcheandj/revision_phewas_plos_genetic/result.single.level.assos/independent_snp_gcta/results_chr${chr}.$pheno.txt \
	--cojo-slct --out /scratch/users/tcheandj/revision_phewas_plos_genetic/result.single.level.assos/independent_snp_gcta/ind.snp.chr${chr}.$pheno 
	" > gcta.cojo.chr${chr}.$pheno.sh 
	#sbatch -p jpriest,normal,owners --qos=normal -J gcta.cojo.chr${chr}.$pheno -t 24:00:00 -N 1 -n 1 -o gcta.cojo.chr${chr}.$pheno.log -e gcta.cojo.chr${chr}.$pheno.err --mem=22000 --open-mode=append export_vcf_subset.$col1.sh;
done


cd /scratch/users/tcheandj/revision_phewas_plos_genetic/result.single.level.assos/independent_snp_gcta/

module load R/3.4.0

R
library(data.table)
library(dplyr)

### import results 
path='/scratch/users/tcheandj/revision_phewas_plos_genetic/result.single.level.assos/independent_snp_gcta/'

### binary phenotypes
files<-list.files(path=path,'*.jma.cojo',full.names=T)

files<- files[!files%in%files[grep(".hybrid.id",files)]]

filesnames<-list.files(path, pattern = '*.jma.cojo')

phenotypes<- lapply(strsplit(filesnames,'\\.'),'[',4)


ind.snp<-lapply(files,fread,header=T,fill=TRUE)
names(ind.snp)<- phenotypes
for (i in 1:length(phenotypes))
  ind.snp[[i]]$pheno=phenotypes[[i]]

ind.snp<-do.call(rbind,ind.snp)

str(ind.snp)

indep_snp<- merge(all.genes,ind.snp,by.y=c("Chr","bp","SNP"),by.x=c("chromosome","position","snp.id"))

### compile all the results from gcta cojo
fwrite(indep_snp,"/scratch/users/tcheandj/revision_phewas_plos_genetic/Ind.snp.set.all.pheno.txt")

