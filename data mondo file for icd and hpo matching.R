library(data.table)
library(dplyr)

setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/phenotypes_genotypes_files')

##############################################
#### file creation for matching in SORTA #####
##############################################

### 1. load my HPO list for syndromic diseases ###
hpo_SD<-fread('hpo_list_SD.txt')
str(hpo_SD)
names(hpo_SD)=c('HPO','HPO_label')

## hpo_SD_file for sorta ###
hpo_SD_sorta<-hpo_SD %>%
  mutate(Names=HPO_label, ID=HPO) %>%
  select(Names,ID)

fwrite(hpo_SD_sorta,'hpo_SD_sorta.csv',sep=';')


### 2. identification subgroup of ICD10 code for SORTA matching ######
icd<-fread('icd10_ukbb_coding19.txt')
### extract label from ID####

icdname<- sapply(strsplit(as.character(icd$Names),' '),'[',-1)

for (i in 1:length(icdname)){
  icd[i,'icdname']<- paste0(icdname[[i]],collapse=' ')
  
}
### select subgroup of icd code
ICD_sub<-icd[nchar(as.vector(icd$ID))==3,]
new_icd<-data.frame(Name=ICD_sub$icdname,ID=ICD_sub$ID)

fwrite(new_icd,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/phenotypes_genotypes_files/ICD_subgroup.csv',sep=';')

#############################################################
###### extraction result of matching HPO and ICD 10 code ####
#############################################################

### 1. matching file from MONDO##
mond_ICD_hpo<-fread('MONDO_ontobio_align_hpo_ICD10.tsv')

str(mond_ICD_hpo)

### icd and hpo extraction

mondo_mapping_ICD_HPO<- mond_ICD_hpo %>%
	mutate(ICD10=sapply(strsplit(as.character(left),'_'),'[',2),ICD10_label=left_label,
	HPO=right,HPO_label=right_label) %>%
	select(ICD10,ICD10_label,HPO,HPO_label)

str(mondo_mapping_ICD_HPO)

## check for duplicated ICD10 code)
mondo_mapping_ICD_HPO.dup<-mondo_mapping_ICD_HPO[duplicated(mondo_mapping_ICD_HPO$HPO),]


### merge my list and mondo matchin
ICD10_hpo_SD<-merge(mondo_mapping_ICD_HPO,hpo_SD)

str(ICD10_hpo_SD)
##check for duplicated in my list
ICD10_hpo_SD.dup<-ICD10_hpo_SD[duplicated(ICD10_hpo_SD$HPO),]


## function to capture ICD10 code that have same HPO terms
#reorder HPO term to have duplicated together
ICD10_hpo_SD<- ICD10_hpo_SD[order(ICD10_hpo_SD$HPO,decreasing = F),]
ICD10_hpo_SD$ICDforHPO<- NA
### function to write ICD code duplicated HPO terms ######
for (i in 1:length(ICD10_hpo_SD$HPO)){
  ICD10_hpo_SD[i,'ICDforHPO']<- ifelse(ICD10_hpo_SD[i,'HPO']==ICD10_hpo_SD[i+1,'HPO'],
  	                                 paste0(ICD10_hpo_SD[i,'ICD10'],',',ICD10_hpo_SD[i+1,'ICD10']),ICD10_hpo_SD[i,'ICDforHPO'])
}

ICD10_hpo_SD$ICDforHPO<- ifelse(is.na(ICD10_hpo_SD$ICDforHPO),ICD10_hpo_SD$ICD10,ICD10_hpo_SD$ICDforHPO)
## remove duplicated
ICD10_hpo_SD.nodup<-ICD10_hpo_SD[!duplicated(ICD10_hpo_SD$HPO),]

##### 2. matching from sorta #####
#### 2.1. matching of hpo_SD to icd10 code 
hpoSDtoicd<-fread('match-result_2018-03-08 00_00_03.csv')

names(hpoSDtoicd)<-c('HPO_names','HPO','ICD_names','ID','score_match','validated')

## 2.2. matching of subgroup icd10 to hpo terms #####

icdsubgptohpo<-fread('match-result_2018-03-13 21_58_40.csv')

## extract hpo terms from the matching 

icdsubgptohpo$HPO_label<- sub('http://purl.obolibrary.org/obo/','',icdsubgptohpo$ontologyTermIRI,fix=T)
icdsubgptohpo$HPO<- sub('_',':',icdsubgptohpo$HPO_label,fix=T)
icdsubgptohpo.2<-icdsubgptohpo[,c('Name','ID','ontologyTermName','HPO','score')]

### merge of the two data matching

sortamatching<-merge(hpoSDtoicd,icdsubgptohpo.2,all=T)

### extract only subgroup of icd10 code ##
#sortamatching.2<-sortamatching[nchar(as.vector(sortamatching$ID))==3,]

### retained matching with 80% score
sorta_ICD_HPO_match<- sortamatching[which(as.numeric(sortamatching$score_match)>=80 | as.numeric(sortamatching$score)>=70),]

### merge with Mondo matching 
all_match<-merge(ICD10_hpo_SD,sorta_ICD_HPO_match,all=T)

all_match.dup<-all_match[duplicated(all_match$HPO),]

### chechk how many HPO terms bilong to syndromic diseases #####
all_match.nodup<-all_match[!duplicated(all_match$HPO),]

match_hpo_SD<-merge(hpo_SD,all_match.nodup,all=T)
### remove row with hpo label missing cause those are not SD
match_hpo_SD<-match_hpo_SD[!is.na(match_hpo_SD$HPO_label),]

### capture possible matching from sorta
match_hpo_SD$ICD10<- ifelse(is.na(match_hpo_SD$ICD10),match_hpo_SD$ID,match_hpo_SD$ICD10)

### remove unmatch HPO  terms
match_hpo_SD<-match_hpo_SD[!is.na(match_hpo_SD$ICD10),]

### check the min number of ICD10 code that I currently have
match_hpo_SD.dup<-match_hpo_SD[duplicated(match_hpo_SD$ICD10),]

##### check the matching rate for usher Syndrome
usher_S<-fread("~/Documents/PheWas_syndromic_diseases/usher_S.csv")

match<-data.frame(HPO=match_hpo_SD$HPO,HPO_label=match_hpo_SD$HPO_label,ICD10=match_hpo_SD$ICD10,labelicd=match_hpo_SD$ICD10_label)
usher_S<-merge(usher_S,match)
