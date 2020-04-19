#### running association on all the phenotypes ###
## for ICD10 code and hpo terms ###
### icd10 code
sbatch -p normal --qos=normal -J plinK_single_snp_assos -t 24:00:00 -N 1 -n 1 -o plinK_single_snp_assos.log -e plinK_single_snp_assos --mem=22000 \
--wrap="/home/users/tcheandj/phewas/plink2 \
--covar /scratch/users/tcheandj/phewas/icd10_phewas/pheno/phewas_icd10_phewas.covar \
  --covar-name PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Age_At_Recruitment Sex.l2 \
  --extract /scratch/users/tcheandj/phewas/icd10_phewas/extract/JAG1/JAG1_chr20_filter_extract.snplist \
  --glm hide-covar firth-fallback  \
  --input-missing-phenotype -999999999 \
  --mac 1 \
  --mach-r2-filter 0.8 2.0 \
  --max-maf 0.05 \
  --memory 22000 \
  --out /scratch/users/tcheandj/phewas/icd10_phewas/single_variant/JAG1_all_pheno/JAG1 \
  --pfile /oak/stanford/projects/ukbb/genotypes/pgen_app1372_hrc/ukb_imp_chr20_v2.mac1.hrc \
  --pheno /scratch/users/tcheandj/phewas/icd10_phewas/pheno/phewas_icd10_phewas.pheno \
  --remove /scratch/users/tcheandj/phewas/icd10_phewas/pheno/phewas_icd10_phewas.remove \
  --threads 2 "

### merge all the result into one file ##
cd /home/users/tcheandj/phewas/single_variant/icd10_result/

### add phenotype in the file as a columns #####

for i in JAG1.*.glm.logistic.hybrid;
do awk ' { x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  ORorBETA      SE      T_STAT  P"  } 
else { gsub("JAG1.","",FILENAME); gsub(".glm.logistic.hybrid","",FILENAME); print FILENAME"\t"$0 }  }  ' $i > ${i}.temp ; mv ${i}.temp ${i}; done

### merge all file 
awk 'FNR>1' JAG1.D61.glm.logistic.hybrid  JAG1.*.glm.logistic.hybrid  > all_result_icd10_jag1.txt

## add header names ##
echo 'outcome  CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  OR      SE      T_STAT  P' > header

cat header all_result_icd10_jag1.txt > logit_icd10_snp_jag1.txt

## linear model ####

for i in JAG1.*.glm.linear;
do awk ' { x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  BETA    SE      T_STAT  P"  } 
else { gsub("JAG1.","",FILENAME); gsub(".glm.linear","",FILENAME); print FILENAME"\t"$0 }  }  ' $i > ${i}.temp ; mv ${i}.temp ${i}; done

## cat all the result none file without header

awk 'FNR>1' JAG1.BMI.glm.linear JAG1.*.glm.linear > all_result_icd10_jag1.txt

echo 'outcome  CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  ORorBETA      SE      T_STAT  P' > header


cat header all_result_icd10_jag1.txt > linear_icd10_snp_jag1.txt

### merge logistic and linear results ###
module load R
R
logit<- read.table('logit_icd10_snp_jag1.txt',header=T,sep='')
linear<- read.table('linear_icd10_snp_jag1.txt',header=T,sep='')

result<- merge(logit,linear,all=T)

write.table(result, 'icd10_jag1_snp_results.txt',quote=F,col.names=T,row.names=F)
q()
n


### analyse for HPO terms hpotems 
cd /scratch/users/tcheandj/phewas/JAG1_phewas/single_variant/JAG1_all_pheno

sbatch -p normal --qos=normal -J plinK_single_snp_assos -t 10:00:00 -N 1 -n 1 -o plinK_single_snp_assos.log -e plinK_single_snp_assos --mem=22000 \
--wrap="/home/users/tcheandj/phewas/plink2 \
--covar /scratch/users/tcheandj/phewas/JAG1_phewas/pheno/phewas_JAG1_phewas.covar \
  --covar-name PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 Age_At_Recruitment Sex.l2 \
  --extract /scratch/users/tcheandj/phewas/JAG1_phewas/extract/JAG1/JAG1_chr20_filter_extract.snplist \
  --glm hide-covar firth-fallback  \
  --input-missing-phenotype -999999999 \
  --mac 1 \
  --mach-r2-filter 0.8 2.0 \
  --max-maf 0.05 \
  --memory 22000 \
  --out /scratch/users/tcheandj/phewas/JAG1_phewas/single_variant/JAG1_all_pheno/JAG1 \
  --pfile /oak/stanford/projects/ukbb/genotypes/pgen_app1372_hrc/ukb_imp_chr20_v2.mac1.hrc \
  --pheno /scratch/users/tcheandj/phewas/JAG1_phewas/pheno/phewas_JAG1_phewas.pheno \
  --remove /scratch/users/tcheandj/phewas/icd10_phewas/pheno/phewas_icd10_phewas.remove \
  --threads 2 "

##### merge all result into one file 
## add file names as colums to each file
cat JAG1.M84.glm.logistic.hybrid > JAG1.M84a.glm.logistic.hybrid 

awk ' { x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  OR      SE      T_STAT  P"  } 
else { gsub("JAG1.","",FILENAME); gsub(".glm.logistic.hybrid","",FILENAME); print FILENAME"\t"$0 }  }  ' JAG1.M84a.glm.logistic.hybrid > JAG1.M84a.glm

for i in JAG1.*.glm.logistic.hybrid;
do awk ' { x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  OR      SE      T_STAT  P"  } 
else { gsub("JAG1.","",FILENAME); gsub(".glm.logistic.hybrid","",FILENAME); print FILENAME"\t"$0 }  }  ' $i > ${i}.temp ; mv ${i}.temp ${i}; done


### append result without the first line

awk 'FNR>1' JAG1.H26.glm.logistic.hybrid JAG1.*.logistic.hybrid > all_result_hpo_jag1.txt

## add header names ##
echo 'CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  OR      SE      T_STAT  P' > header

cat header all_result_hpo_jag1.txt > result_hpo_snp_jag1.txt

## linear model ####

for i in JAG1.*.glm.linear;
do awk ' { x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  ORorBETA      SE      T_STAT  P"  } 
else { gsub("JAG1.","",FILENAME); gsub(".glm.linear","",FILENAME); print FILENAME"\t"$0 }  }  ' $i > ${i}.temp ; mv ${i}.temp ${i}; done

## cat all the result none file without header

awk 'FNR>1' JAG1.BMI.glm.linear JAG1.*.glm.linear > all_result_icd10_jag1.txt

echo 'outcome  CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  ORorBETA      SE      T_STAT  P' > header


cat header all_result_icd10_jag1.txt > linear_icd10_snp_jag1.txt

### merge logistic and linear results ###
module load R
R
logit<- read.table('logit_icd10_snp_jag1.txt',header=T,sep='')
linear<- read.table('linear_icd10_snp_jag1.txt',header=T,sep='')

result<- merge(logit,linear,all=T)

write.table('icd10_jag1_snp_results.txt',quote=F,col.names=T,row.names=F)
q()
n


### transfer file tomy directory

scp tcheandj@login.sherlock.stanford.edu://home/users/tcheandj/phewas/single_variant/icd10_result/icd10_jag1_snp_results.txt /Users/catherine.t/Documents/PHWAS_doc/result/

scp tcheandj@login.sherlock.stanford.edu://scratch/users/tcheandj/phewas/JAG1_phewas/single_variant/JAG1_all_pheno/result_hpo_snp_jag1.txt /Users/catherine.t/Documents/PHWAS_doc/result/




