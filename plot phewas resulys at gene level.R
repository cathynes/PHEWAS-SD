### plot PHEWAS results by syndromes
library(dplyr)
library(ggplot2)
library(PheWAS)
library(grid)
library(lattice)
library(qqman)
library(data.table)

##### load data set ######
setwd('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/select_hpoS/phewas_results')


### single gene result
single_gene<-fread('phewas_SD.all_loci.all_phenotypes.skat.txt.meta_datav0723.txt',header = T,sep='\t',fill=TRUE)


single_gene<-filter(single_gene, set=='all',weight=='maf',!is.na(p_value))

## QQplot  gene level all results 
qq(single_gene$p_value,xlim = c(0,10), ylim = c(0,10),main='QQplot gene level all results (Lambda = 1.34)')


single_gene$label<-ifelse(single_gene$p_value < 1e-4, paste(as.character(single_gene$outcome),''),'')

phwas.plot<- ggplot(single_gene,aes(x=locus_name, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=4,colour='red') +
  labs(title='Association between selected gene and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=single_gene$label),size=3,vjust=-1.1,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


### single gene result
single_gene<-fread('phewas_SD.all_loci.all_phenotypes.skat.txt.meta_datav0723.txt',header = T,sep='\t',fill=TRUE)


single_gene<-filter(single_gene, set=='all',weight=='weight_cadd',!is.na(p_value))

## QQplot  gene level all results 
qq(single_gene$p_value,xlim = c(0,10), ylim = c(0,10),main='QQplot gene level all results (Lambda = 1.34)')


single_gene$label<-ifelse(single_gene$p_value < 1e-4, paste(as.character(single_gene$outcome),''),'')

phwas.plot<- ggplot(single_gene,aes(x=locus_name, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=4,colour='red') +
  labs(title='Association between selected gene and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=single_gene$label),size=3,vjust=-1.1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot
#### ALagille syndromes

##SKAT test 
AS<- single_gene %>% filter(locus_name=='JAG1' |locus_name=='NOTCH2')

## qqplot 
### inflation 
chisq <- qchisq(1-AS$p_value,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 


AS$label<-ifelse(AS$p_value < 1e-4,paste(as.character(AS$locus_name),as.character(AS$outcome),sep='*'),'')


phwas.plot<- ggplot(AS,aes(x=outcome, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=4,colour='red') +
  labs(title='Association between SNPs in JAG1 and NOTCH2, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=AS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

#### Marfan syndromes

MF<- single_gene %>% filter(locus_name=='FBN1' )

## qqplot 
### inflation 
chisq <- qchisq(1-MF$p_value,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(MF$p_value,xlim = c(0,8),main='QQplot SNPs level Associations MS (FBN1) (Lambda = 1.18)')

MF$label<-ifelse(MF$P < 10e-20,paste(as.character(MF$ID),as.character(MF$outcome),sep='*'),'')


phewas.plot<- ggplot(MF,aes(x=outcome, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Association between SNPs in FBN1, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=MF$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot


##### Noonan syndrome 

### nonan S genes
gene<-c('PTPN11', 'SOS1', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B')

NS<- single_gene %>% filter(locus_name%in%gene )

## qqplot 
### inflation 
chisq <- qchisq(1-NS$p_value,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(NS$p_value,xlim = c(0,10),main='QQplot SNPs level Associations NS gene group (Lambda = 1.12)')

NS$label<-ifelse(NS$p_value < 10e-8,paste(as.character(NS$locus_name),as.character(NS$outcome),sep='*'),'')


phewas.plot<- ggplot(NS,aes(x=outcome, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=5,colour='red') + 
  labs(title='Association between SNPs Noonan, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot

### nonan resticted on 3 gene
gene<-c('PTPN11', 'RASA2', 'MAP2K1', 'KAT6B')

NS<- single_gene %>% filter(locus_name%in%gene )

## qqplot 
### inflation 
chisq <- qchisq(1-NS$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(NS$P,xlim = c(0,10),main='QQplot SNPs level Associations NS gene group (Lambda = 1.25)')

NS$label<-ifelse(NS$P < 10e-20,paste(as.character(NS$locus_name),as.character(NS$outcome),sep='*'),'')


phewas.plot<- ggplot(NS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') +
  labs(title='Association between SNPs in Noonan S gene group, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot

##### Digeoge syndrome

DS<- single_gene %>% filter(!locus_name%in%gene)

## qqplot 
### inflation 
chisq <- qchisq(1-DS$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(DS$P,xlim = c(0,8),main='QQplot SNPs level Associations DS gene in 22p11 deletion (Lambda = 1.03)')

DS$label<-ifelse(DS$p_value < 10e-8,paste(as.character(DS$locus_name),as.character(DS$outcome),sep='*'),'')


phewas.plot<- ggplot(DS,aes(x=outcome, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Association between SNPs in Digeorge S gene group, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot


