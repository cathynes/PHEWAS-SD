### plot for results for phewas on syndromic diseases
library(dplyr)
library(ggplot2)
library(PheWAS)
library(grid)
library(lattice)
library(qqman)
library(data.table)


setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS/phewas_results')

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

##manhattan plot for select 
par(mfrow=c(1,1))
manhattan(single_snp,chr='#CHROM',bp='POS',P='P',snp="ID",chrlabs = c('1','2','3','7','10','11','12',
                                                                     '14','19','20','22'))

single_snp<- single_snp[order(single_snp$`#CHROM`,single_snp$POS ,decreasing = F),]

single_snp$rownames<-rep(1:nrow(single_snp))

phwas.plot<- ggplot(single_snp,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


max<-single_snp %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))

qq(single_snp$P)


#@##plot fo HP0000821

HP0000821<- single_snp %>% filter(outcome=='HP0000821')
HP0000821<- HP0000821[order(HP0000821$`#CHROM`,HP0000821$POS ,decreasing = F),]

HP0000821$rownames<-rep(1:nrow(HP0000821))

phwas.plot<- ggplot(HP0000821,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

max<-HP0000821 %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))



manhattan(HP0000821,chr='#CHROM',bp='POS',P='P',snp="ID")
qq(HP0000821$P,xlim = c(0,10), ylim = c(0,10))

HP0005117<- single_snp %>% filter(outcome=='HP0005117')
manhattan(HP0005117,chr='#CHROM',bp='POS',P='P',snp="ID")
qq(HP0005117$P)

HP0005117<- HP0005117[order(HP0005117$`#CHROM`,HP0005117$POS ,decreasing = F),]

HP0005117$rownames<-rep(1:nrow(HP0005117))

phwas.plot<- ggplot(HP0005117,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

max<-HP0005117 %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))



HP0000002<- single_snp %>% filter(outcome=='HP0000002')
manhattan(HP0000002,chr='#CHROM',bp='POS',P='P',snp="ID")
qq(HP0000002$P)


HP0000002<- HP0000002[order(HP0000002$`#CHROM`,HP0000002$POS ,decreasing = F),]

HP0000002$rownames<-rep(1:nrow(HP0000002))

phwas.plot<- ggplot(HP0000002,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

max<-HP0000002 %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))


### small gestational age

HP0001518<- single_snp %>% filter(outcome=='HP0001518')
manhattan(HP0001518,chr='#CHROM',bp='POS',P='P',snp="ID")
qq(HP0001518$P)

HP0001518<- HP0001518[order(HP0001518$`#CHROM`,HP0001518$POS ,decreasing = F),]

HP0001518$rownames<-rep(1:nrow(HP0001518))

phwas.plot<- ggplot(HP0001518,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

max<-HP0001518 %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))

### growth abnormality

HP0001507<- single_snp %>% filter(outcome=='HP0001507')
manhattan(HP0001507,chr='#CHROM',bp='POS',P='P',snp="ID")
qq(HP0001507$P)

HP0001507<- HP0001507[order(HP0001507$`#CHROM`,HP0001507$POS ,decreasing = F),]

HP0001507$rownames<-rep(1:nrow(HP0001507))

phwas.plot<- ggplot(HP0001507,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

max<-HP0001507 %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))


####PTPN11 and all phenotype


PTPN11<- single_snp %>% filter(locus_name=='PTPN11')
qq(PTPN11$P)
PTPN11$label<-ifelse(PTPN11$P < 10e-10,paste0(as.character(PTPN11$ID)),'')


phwas.plot<- ggplot(PTPN11,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=PTPN11$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

####JAG1 and all phenotype


JAG1<- single_snp %>% filter(locus_name=='JAG1')
qq(JAG1$P)
JAG1$label<-ifelse(JAG1$P < 10e-10,paste0(as.character(JAG1$ID)),'')


phwas.plot<- ggplot(JAG1,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,17)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=JAG1$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


####FBN1 and all phenotype


FBN1<- single_snp %>% filter(locus_name=='FBN1')
qq(FBN1$P)
FBN1$label<-ifelse(FBN1$P < 10e-10,paste0(as.character(FBN1$ID)),'')


phwas.plot<- ggplot(FBN1,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,15)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='FBN1',y='-log10(pvalue)') +
  #geom_text(aes(label=FBN1$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot



HP0012594<- single_snp %>% filter(outcome=='HP0012594')
HP0004421<- single_snp %>% filter(outcome=='HP0004421')

par(mfrow=c(2,2))
manhattan(HP0000821,chr='CHROM',bp='POS',P='P',snp="ID",
          ,main='plot hypothyroidism (HP0000821)')
qq(HP0000821$P,main='plot hypothyroidism (HP0000821)')



manhattan(HP0000002,chr='CHROM',bp='POS',P='P',snp="ID",
          main='standing/sitting height (HP0000002)')
qq(HP0000002$P,main='standing/sitting height (HP0000002)')
boxplot(x=phenotype$HP0000002)
hist(x=phenotype$HP0000002, breaks=600)


manhattan(HP0005117,chr='CHROM',bp='POS',P='P',snp="ID",
          main='plot diastolic BP (HP0005117)')
qq(HP0005117$P,main='plot diastolic BP (HP0005117)')
boxplot(x=phenotype$HP0005117)
hist(x=phenotype$HP0005117, breaks=600)


manhattan(HP0004421,chr='CHROM',bp='POS',P='P',snp="ID",
          main='plot systolic BP (HP0004421)')
qq(HP0004421$P,main='plot systolic BP (HP0005117)')
boxplot(x=phenotype$HP0004421)
hist(x=phenotype$HP0004421, breaks=600)



manhattan(HP0012594,chr='CHROM',bp='POS',P='P',snp="ID",
          main='urine creatinine level(HP0012594)')
qq(HP0012594$P,main='urine creatinine level(HP0012594)')

boxplot(x=phenotype$HP0012594)
hist(x=phenotype$HP0012594, breaks=600)


########################
### gene level plot ####

single_gene<-fread('phewas_SD.all_loci.all_phenotypes.skat.txt.meta_datav0723.txt',header = T,sep='\t')


single_gene$label<-ifelse(single_gene$p_value < 10e-8,paste0(as.character(single_gene$outcome)),'')

single_gene<- single_gene[order(single_gene$locus_name ,decreasing = F),]

single_gene$rownames<-rep(1:length(single_gene$locus_name))

## get max score for each gene
max<-single_gene %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

## change the level to have level as on the plot 

### manhattan plot with gene as color
phwas.plot<- ggplot(single_gene,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  geom_text(aes(label=single_gene$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))
qq(single_gene$p_value,main='qqplot gene results')

## get max score for each gene
max<-single_gene %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

## change the level to have level as on the plot 
single_gene<- single_gene[order(single_gene$cat_text ,decreasing = F),]

### manhattan plot with gene as color

single_gene$label<-ifelse(single_gene$p_value < 10e-5,paste0(as.character(single_gene$pheno_text)),'')

phwas.plot<- ggplot(single_gene,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=cat_text),size=2,show.legend = T) +
  geom_hline(yintercept=5,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  geom_text(aes(label=single_gene$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))


### manhattan plot for just set=all and weight CAD the all and MAF 

single_gene1<- filter(single_gene,set=='all',weight=='maf')

phwas.plot<- ggplot(single_gene1,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,11)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms (maf and all snps methods',x='gene',y='-log10(pvalue)') +
  geom_text(aes(label=single_gene1$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))
qq(single_gene1$p_value,xlim = c(0,8), ylim = c(0,10))

#### set all and weigt CAD
single_gene2<- filter(single_gene,set=='all',weight=='weight_cadd')

phwas.plot<- ggplot(single_gene2,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms (weight Cadd and all snps methods',x='gene',y='-log10(pvalue)') +
  geom_text(aes(label=single_gene2$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))
qq(single_gene2$p_value,xlim = c(0,50), ylim = c(0,51))

p2<-single_gene2[!is.na(single_gene2$p_value) ,'p_value']
observed2 <- sort(p2)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2)) 
lexp2 <- -(log10(expected2 / (length(expected2)+1)))

plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", 
     ylab="Observed (-logP)", xlim=c(0,6), ylim=c(0,15), las=1, xaxs="i", yaxs="i", bty="l",main="")
points(lexp, lobs, pch=25, cex=.4, bg="blue",col='blue') 
points(lexp2, lobs2, pch=25, cex=.4, bg="blue",col='blue') 
title(main=list('QQplot PheWas all gene with MAF and weight methods for  select HPO terms ',cex=1,col='black'))

############################################
##### plot result for snp threasure <0.5 ###
#############################################

setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/tcheandj_phewas_phewas.SD.co.ra.variants')

### single snp result
single_snp<-read.table('phewas.SD.co.ra.variants_single_variant.all_loci.all_phenotypes.txt',header = T,sep='\t')

single_snp$label<-ifelse(single_snp$P<10e-8,paste0(as.character(single_snp$ID),'*',
                                                    as.character(single_snp$outcome)),'')

single_snp<- single_snp[order(single_snp$CHROM , single_snp$POS,decreasing = F),]

single_snp$rownames<-rep(1:length(single_snp$CHROM))

## get max score for each gene
max<-single_snp %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]
## change the level to have level as on the plot 
levels(single_snp$locus_name)
single_snp$locus_name <- factor(single_snp$locus_name, 
                                levels = c('NOTCH2', 'RIT1','SOS1','RAF1','BRAF' ,'KAT6B','A2ML1','KRAS' ,'PTPN11',
                                           'SOS2' ,'FBN1' ,'MAP2K1','RRAS', 'JAG1', 'TBX1'  ,'LZTR1' ))

### manhattan plot with gene as color
phwas.plot<- ggplot(single_snp,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,89)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_snp$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))



#### qqplot snp level ####
####qq plot at gene level
p<-single_snp[!is.na(single_snp$P) ,'P']
observed <- sort(p)
lobs <- -log10(observed)

expected <- c(1:length(observed)) 
lexp <- -log10(expected / (length(expected)+1))

plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", 
     ylab="Observed (-logP)", xlim=c(0,6), ylim=c(0,89), las=1, xaxs="i", yaxs="i", bty="l",main="")
points(lexp, lobs, pch=25, cex=.4, bg="blue",col='blue') 
title(main=list('QQplot PheWas all snp in all gene with select HPO terms ',cex=1,col='black'))

### gene level plot ####

single_gene<-read.table('phewas.SD.co.ra.variants.all_loci.all_phenotypes.skat.txt.meta_data.txt',header = T,sep='\t')


single_gene$label<-ifelse(single_gene$p_value < 10e-5,paste0(as.character(single_gene$locus_name),'*',
                                                             as.character(single_gene$outcome)),'')

single_gene<- single_gene[order(single_gene$locus_name ,decreasing = F),]

single_gene$rownames<-rep(1:length(single_gene$locus_name))

## get max score for each gene
max<-single_gene %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

## change the level to have level as on the plot 

### manhattan plot with gene as color
phwas.plot<- ggplot(single_gene,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=5,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))

## get max score for each gene
max<-single_gene %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

## change the level to have level as on the plot 
single_gene<- single_gene[order(single_gene$cat_text ,decreasing = F),]

### manhattan plot with gene as color

single_gene$label<-ifelse(single_gene$p_value < 10e-5,paste0(as.character(single_gene$pheno_text)),'')

phwas.plot<- ggplot(single_gene,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=cat_text),size=2,show.legend = F) +
  geom_hline(yintercept=5,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))


### manhattan plot for just set=all and weight CAD the all and MAF 

single_gene1<- filter(single_gene,set=='all',weight=='maf')

phwas.plot<- ggplot(single_gene1,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=cat_text),size=2,show.legend = T) +
  geom_hline(yintercept=5,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms (maf and all snps methods',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene1$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))


####qq plot at gene level
p<-single_gene1[!is.na(single_gene1$p_value) ,'p_value']
observed <- sort(p)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", 
     ylab="Observed (-logP)", xlim=c(0,6), ylim=c(0,11), las=1, xaxs="i", yaxs="i", bty="l",main="")
points(lexp, lobs, pch=25, cex=.4, bg="blue",col='black') 
title(main=list('QQplot PheWas all gene ith MAF methods and  select HPO terms ',cex=1,col='black'))


#### set all and weigt CAD
single_gene2<- filter(single_gene,set=='all',weight=='weight_cadd')

phwas.plot<- ggplot(single_gene2,aes(x=rownames, y=-log10(p_value))) + geom_point(aes(col=cat_text),size=2,show.legend = F) +
  geom_hline(yintercept=5,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms (weight Cadd and all snps methods',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_gene2$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))

p2<-single_gene2[!is.na(single_gene2$p_value) ,'p_value']
observed2 <- sort(p2)
lobs2 <- -(log10(observed2))

expected2 <- c(1:length(observed2)) 
lexp2 <- -(log10(expected2 / (length(expected2)+1)))

plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", 
     ylab="Observed (-logP)", xlim=c(0,6), ylim=c(0,15), las=1, xaxs="i", yaxs="i", bty="l",main="")()
points(lexp, lobs, pch=25, cex=.4, bg="blue",col='blue') 
points(lexp2, lobs2, pch=25, cex=.4, bg="black",col='black') 
title(main=list('QQplot PheWas all gene with MAF and weight methods for  select HPO terms ',cex=1,col='black'))