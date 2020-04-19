#### phenotype file for Phewas ######
library(data.table)
library(dplyr)

setwd("/oak/stanford/groups/jpriest/catherine")

set_pheno1<- fread('Pheno_phewas.selectS_set1.txt')

set_pheno1$ID<-NULL
set_pheno1$app15860<-NULL


set_pheno2<- fread('Pheno_phewas.selectS_idbbeng.txt')

##id correspondance

id_manu<- read.table('ukb24983_13721_mapping.tsv',sep='\t',header=F)

names(id_manu)<- c("IID","app13721")

set_pheno2<- merge(id_manu,set_pheno2)

str(set_pheno2)

### merge both files ####
final_pheno_set1<-merge(set_pheno2[,1:11],set_pheno1)


str(final_pheno_set1)

mydata <- final_pheno_set1 %>% 
     mutate_if(is.integer,as.factor) %>%
	 mutate_at(vars(HP0001513:HP0000707), funs(recode(., '1'='0', '2'='1')))

str(mydata)
 
### recode binary phenotypes in the set_pheno2 

mydata2 <- set_pheno2 %>% 
     mutate_if(is.integer,as.factor) %>%
	 mutate_at(vars(HP0000023:HP0100753), funs(recode(., '1'='0', '2'='1'))) 
str(mydata2)

### combine the duplicated variable into one 

mydata2$HP0002107 <- ifelse(mydata2$HP0002107='1' | mydata2$HP0002107.1='1',"1",0)

### verif 
table(mydata2$HP0002107)
table(mydata2$HP0002107,mydata2$HP0002107.1)

## delete the extra varible
mydata2$HP0002107.1<-NULL


## merge all data
finalphenodata<- merge(mydata,mydata2)

str(finalphenodata)

### convert HP0003081 into numeric
finalphenodata <- finalphenodata %>%
		mutate(HP0003081=as.numeric(HP0003081), HP0000120=as.numeric(HP0000120),
			sex=as.integer(sex), Array=as.integer(Array))

str(finalphenodata)

### write the final table

write.table(finalphenodata,'SD_pheno_data.txt',quote=F,sep='\t',row.names=F)


