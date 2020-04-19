### creation of meta data and phenotype description file for pheas SD ####
### creation final pheno
library(data.table)
library(dplyr)

setwd("/oak/stanford/groups/jpriest/catherine")

pheno.SD<- fread('UKBB.SDpheno.in.hpocode.txt')

str(pheno.SD)

## verif valeur HP0000120
table(pheno.SD$HP0000120)

### creation phenotype desc file

pheno.desc<- data.frame(name=names(pheno.SD),type=c(0,"C","D",rep("C",6),rep("P",15),rep("B",54)))


### meta data ####
annotation.file<- fread('hpoSD.and.icd.ukfield.bieng.match.txt')

### correspondance with current file ####

annotation.file$name<- annotation.file$hpoid

annotation.phenodes.file<- merge(pheno.desc,annotation.file,all=T)

str(meta.file)
### remove covariate (they are missing for HPOid)

meta.file<- annotation.phenodes.file %>%
			filter(!is.na(hpoid)) %>%
			mutate(pheno_text=hpolabel2,cui_id=bioen.icd10.id,ukb_id=code.origine,
				cat_text=syndrome,cat_index_1=1,cat_index_2=2) %>%
			select(name,cui_id,ukb_id,pheno_text,cat_text,cat_index_1, cat_index_2)


#### save each data #####


fwrite(pheno.desc,'/oak/stanford/groups/jpriest/catherine/new.phewas.SD/pheno_desc.txt',sep='\t')

fwrite(meta.file,'/oak/stanford/groups/jpriest/catherine/new.phewas.SD/meta.data.txt',sep='\t')


fwrite(spheno.SD,'/oak/stanford/groups/jpriest/catherine/new.phewas.SD/phenotype.txt',sep='\t')








