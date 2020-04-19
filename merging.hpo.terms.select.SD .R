#### creation metata data and pheno desc data for phewas pipeline ####
library(data.table)
library(dplyr)
library(tidyr)
library(VennDiagram)
library(colorfulVennPlot)
####create a phenotype list for the select SD

setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS')

#### HPO terms are some time duplicated from the data download in HPO website
Alagille<-fread('HPO_Alagille.csv',header=T,na.strings = '')
names(Alagille)<- c('hpoid','hpolabel','Alagille')

### remove duplicated hpo
Alagille<-Alagille %>% filter(!duplicated(hpoid))

Digeorge<-fread('HPO_Digeorge.csv',header=T)
names(Digeorge)<- c('hpoid','hpolabel','Digeorge')
### remove duplicated hpo
Digeorge<-Digeorge %>% filter(!duplicated(hpoid))


Noonan<-fread('HPO_Noonan.csv',header=T,na.strings = '')
names(Noonan)<- c('hpoid','hpolabel','Noonan')
### remove duplicated hpo
Noonan<-Noonan %>% filter(!duplicated(hpoid))


Marfan<-fread('HPO_Marfan.csv',header=T,na.strings = '')
names(Marfan)<- c('hpoid','hpolabel','Marfan')
### remove duplicated hpo
Marfan<-Marfan %>% filter(!duplicated(hpoid))


### merge of the 4 data
hpo_SD1 <- merge(Alagille, Digeorge,all=T)
hpo_SD2 <- merge(Noonan, Marfan,all=T)

hpo_SD<-merge(hpo_SD1,hpo_SD2,all=T)

#Put all the syndromic disease together
## replace all NA by empty space
#hpo_SD <- data.frame(lapply(hpo_SD, function(x) { x[is.na(x)] <- "" ; x })) no need 

hpo_SD$syndrome=apply(hpo_SD[,3:6], 1, function(x) paste(x[!is.na(x) & x != "No"], collapse = ","))
table(hpo_SD$syndrome)
## cluster
hpo_SD<- hpo_SD%>% mutate(Alagille1=ifelse(is.na(Alagille),0,1),Digeorge1=ifelse(is.na(Digeorge),0,1),
                          Noonan1=ifelse(is.na(Noonan),0,1),Marfan1=ifelse(is.na(Marfan),0,1))

aa<-hpo_SD %>% mutate(counts = rowSums(.[8:11])) %>% select(Alagille1:counts)
row.names(aa)<-hpo_SD$hpoid

##### try a ven diagramm base on HPO terms
## create the list of hpo
HPOSD<- list(Alagille=Alagille$hpolabel,Marfan=Marfan$hpolabel,Noonan=Noonan$hpolabel,Digeorge=Digeorge$hpolabel)


venplot<- venn.diagram(HPOSD,'venplotas.tiff',scaled=T,ext.text=T, fill=c('orange','blue','green','red'),
                       main='ven plot HPO terms among syndrome')
overlap<-calculate.overlap(HPOSD)

common<-data.frame(hpoid=overlap$a6)

venn.plot(venplot)

overlaphpo<-merge(hpo_SD,common)

### remove in initial hpo data all the term that 'abnormality of"

hpo_SD2<- hpo_SD %>% filter(!grepl('Abnormality of',hpolabel))

HPOSD2<- list(Alagille3=select(Alagille,hpolabel%in%hpo_SD2$hpolabel),
              Marfan=Marfan$hpolabel%in%hpo_SD2$hpolabel,
              Noonan=Noonan$hpolabel%in%hpo_SD2$hpolabel,
              Digeorge=Digeorge$hpolabel%in%hpo_SD2$hpolabel)

venplot2<- venn.diagram(HPOSD2,'venplotas2.tiff',scaled=T,ext.text=T, fill=c('orange','blue','green','red'),
                       main='ven plot HPO terms among syndrome2')


#### upload clean version of matching

set_pheno2<-fread("clean_HPO_bioenginematching2.txt")

### merge the clean file with hpo_SDv2
hpoSD.and.match<- merge(hpo_SD,set_pheno2,all=T)

### fill mising hpolabel
hpoSD.and.match2<- hpoSD.and.match %>% 
  mutate(hpolabel2=ifelse(is.na(hpolabel),as.character(HPO.label),as.character(hpolabel))) %>% 
  select(hpoid,hpolabel2,syndrome,bioen.icd10.id,bioeng.or.icd.label,code.origine)

###write the final table
fwrite(hpoSD.and.match2,'hpoSD.and.icd.bieng.match.txt',sep='\t')
fwrite(hpo_SD,'hpo_AS_MF_NS_DS.txt',sep='\t')



