### plot PHEWAS results by syndromes
library(dplyr)
library(ggplot2)
library(PheWAS)
library(grid)
library(lattice)
library(qqman)
library(data.table)


##### load data set ######
setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS/phewas_results')


### single snp result
single_snp<-fread('phewas_SD_single_variant.all_loci.all_phenotypesv0723.txt',header = T,sep='\t',fill=TRUE)



### single snp result
single_snp<-fread('phewas_SD_single_variant.all_loci.all_phenotypesv0723.txt',header = T,sep='\t',fill=TRUE)

#### extract pheno that have missing value for p
single_snp2<- single_snp %>% filter(is.na(P)) %>% select(locus_name:BETA)
## reharmonized header for single_snp2
names(single_snp2)<- c("locus_name" , "outcome","#CHROM","POS","ID","REF","ALT" ,
                       "ALT_FREQ","MACH_R2","OBS_CT","BETA","SE","P" )

single_snp3<- single_snp %>% filter(!is.na(P))
## remerge data
single_snp<- merge(single_snp2,single_snp3,all=T)

###manhattan and QQ plot of snp level result
single_snp<- single_snp %>% filter(!is.na(single_snp$P))

### manhattan plot snp level

single_snp<- single_snp[order(single_snp$`#CHROM`,single_snp$POS,decreasing = F),]

single_snp$rownames<-rep(1:nrow(single_snp))

single_snp$label<-ifelse(single_snp$P < 5e-30,paste(as.character(single_snp$ID),as.character(single_snp$outcome),sep='*'),'')

phwas.plot<- ggplot(single_snp,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=outcome),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  geom_text(aes(label=single_snp$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


max<-single_snp %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))

### sort by chr pos and outcome 
single_snp<- single_snp[order(single_snp$`#CHROM`,single_snp$POS,single_snp$outcome,decreasing = F),]

single_snp$rownames<-rep(1:nrow(single_snp))

single_snp$label<-ifelse(single_snp$P < 5e-8,paste(as.character(single_snp$ID),as.character(single_snp$outcome),sep='*'),'')

phwas.plot<- ggplot(single_snp,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  geom_text(aes(label=single_snp$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


max<-single_snp %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))

qq(single_snp$P)

#### ALagille syndromes

AS<- single_snp %>% filter(locus_name=='JAG1' |locus_name=='NOTCH2' )

## qqplot 
### inflation 
chisq <- qchisq(1-AS$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(AS$P,xlim = c(0,8), ylim = c(0,16),main='QQplot SNPs level associations AS (JAG1, NOTCH2) (Lambda = 1.32)')

AS$label<-ifelse(AS$P < 10e-12,paste(as.character(AS$ID),as.character(AS$outcome),sep='*'),'')


phwas.plot<- ggplot(AS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,17)+
  labs(title='Association between SNPs in JAG1 and NOTCH2, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=AS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) +
  scale_color_manual(breaks = c(levels(as.factor(AS$locus_name))), values=c("orange", "blue"))
phwas.plot

#### Marfan syndromes

MF<- single_snp %>% filter(locus_name=='FBN1' )

## qqplot 
### inflation 
chisq <- qchisq(1-MF$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(MF$P,xlim = c(0,8),main='QQplot SNPs level Associations MS (FBN1) (Lambda = 1.18)')

MF$label<-ifelse(MF$P < 10e-20,paste(as.character(MF$ID),as.character(MF$outcome),sep='*'),'')


phewas.plot<- ggplot(MF,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Association between SNPs in FBN1, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=MF$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot


##### Noonan syndrome 

### nonan S genes
gene<-c('PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B')

NS<- single_snp %>% filter(locus_name%in%gene )

## qqplot 
### inflation 
chisq <- qchisq(1-NS$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(NS$P,xlim = c(0,10),main='QQplot SNPs level Associations NS gene group (Lambda = 1.12)')

NS$label<-ifelse(NS$P < 10e-20,paste(as.character(NS$ID),as.character(NS$outcome),sep='*'),'')


phewas.plot<- ggplot(NS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Association between SNPs Noonan, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot

### nonan resticted on 3 gene
gene<-c('PTPN11', 'RASA2', 'MAP2K1', 'KAT6B')

NS<- single_snp %>% filter(locus_name%in%gene )

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

DS<- single_snp %>% filter(`#CHROM`==22)

## qqplot 
### inflation 
chisq <- qchisq(1-DS$P,1)
lambda = median(chisq)/qchisq(0.5,1)
lambda 

qq(DS$P,xlim = c(0,8),main='QQplot SNPs level Associations DS gene in 22p11 deletion (Lambda = 1.03)')

DS$label<-ifelse(DS$P < 10e-8,paste(as.character(DS$locus_name),as.character(DS$outcome),sep='*'),'')


phewas.plot<- ggplot(DS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Association between SNPs in Digeorge S gene group, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phewas.plot


