####PTPTN11 gene

###1.1 logistic model ###

for i in single.variant.PTPN11.*.glm.logistic.hybrid;
do awk '{ x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  ORorBETA      SE      T_STAT  P"  } 
else { gsub("single.variant.PTPN11.","",FILENAME); gsub(".glm.logistic.hybrid","",FILENAME); print FILENAME"\t"$0 }  }' $i > ${i}.temp ; mv ${i}.temp ${i}; done

## 1.2 linear model ####

for i in single.variant.PTPN11.*.glm.linear;
do awk ' { x=1 ; if ( x == NR ) { print "Outcome CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  ORorBETA      SE      T_STAT  P"  } 
else { gsub("PTPN11.","",FILENAME); gsub(".glm.linear","",FILENAME); print FILENAME"\t"$0 }  }  ' $i > ${i}.temp ; mv ${i}.temp ${i}; done

#### 2. append all result 

###2.1 logistic model ###

awk 'FNR>1'  single.variant.PTPN11.*.logistic.hybrid > all_result_SD_PTPN11.txt #append without headr

echo 'outcome  CHROM  POS     ID      REF     ALT     A1      FIRTH  TEST    OBS_CT  OR      SE      T_STAT  P' > header ## create header

cat header all_result_SD_PTPN11.txt > logit_SD_snp_PTPN11.txt ## append the header


## 2. 2 linear model

awk 'FNR>1' single.variant.PTPN11.*.glm.linear > all_linear_SD_PTPN11.txt

## select few line 

cat all_linear_SD_PTPN11.txt| cut  -f 3,7-19 >  all_linear_SD_PTPN11v2.txt

echo 'outcome  CHROM  POS     ID      REF     ALT     A1      TEST    OBS_CT  ORorBETA      SE      T_STAT  P' > header

cat header all_linear_SD_PTPN11v2.txt> linear_SD_snp_PTPN11.txt

### 3. merge logistic and linear results ###
module load R
R
logit<- read.table('logit_SD_snp_PTPN11.txt',header=T,sep='')
linear<- read.table('linear_SD_snp_PTPN11.txt',header=T,sep='')

result<- merge(logit,linear,all=T)

write.table('SD_PTPN11_snp_results.txt',quote=F,col.names=T,row.names=F)
q()
n
