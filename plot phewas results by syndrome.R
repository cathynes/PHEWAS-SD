### plot for results for phewas on syndromic diseases
library(dplyr)
library(ggplot2)
library(grid)
library(lattice)
library(qqman)
library(data.table)
library(gridExtra)
library(wesanderson)
library(RColorBrewer)
library(PheWAS)
###plot results by syndrome SNP and gene level 
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

### manhattan plot snp level

single_snp<- single_snp[order(single_snp$`#CHROM`,single_snp$POS,decreasing = F),]

single_snp$rownames<-rep(1:nrow(single_snp))

single_snp$label<-ifelse(single_snp$P < 5e-08,paste(as.character(single_snp$ID),as.character(single_snp$outcome),sep='*'),'')

#### all result swith color by syndrom
meta.data<-fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/meta.data.txt')
meta.data$outcome<-meta.data$name
single.snp.metadata<-merge(single_snp,meta.data,all=T)

table(single.snp.metadata$cat_text,exclude=NULL)

single.snp.metadata$cat_text<-ifelse(is.na(single.snp.metadata$cat_text) | single.snp.metadata$cat_text=='','Alagille syndrome, DIGEORGE syndrome, Noonan syndrome, Marfan syndrome',
                                     single.snp.metadata$cat_text)

table(single.snp.metadata$cat_text,exclude=NULL)

single.snp.metadata<- single.snp.metadata[order(as.factor(single.snp.metadata$cat_text),decreasing = F),]

single.snp.metadata$rownames<-rep(1:nrow(single.snp.metadata))
fwrite(single.snp.metadata,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/single.snp.metadata.txt')

single.snp.metadata2<- single.snp.metadata %>% mutate(SNP=ID,pval=P,phenotype=paste0(outcome,'(',pheno_text,')')) %>% 
                        select(SNP,phenotype,pval,locus_name)

fwrite(single.snp.metadata2,'/Users/catherine.t/Documents/PheWas_syndromic_diseases/tableau_results_phewas/single.snp.metadata2.txt',sep='\t')


phwas.plot<- ggplot(single.snp.metadata,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=cat_text),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + ylim(0,20)+
  labs(title='',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=single.snp.metadata$label),size=4,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name')
p1<-phwas.plot
p1

### plot Alagille Syndrome #####
AS<- single_snp %>% filter(locus_name=='JAG1' |locus_name=='NOTCH2' )

AS$label<-ifelse(AS$P < 10e-6,paste(as.character(AS$ID),as.character(AS$outcome),sep='*'),'')

phwas.plot<- ggplot(AS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,20)+
  labs(title='',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=AS$label),size=4,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(AS$locus_name))), values=c("Brown", "blue"))
p1<-phwas.plot
p1
### marfan syndrome
MF<- single_snp %>% filter(locus_name=='FBN1' )
MF$label<-ifelse(MF$P < 10e-8,paste(as.character(MF$ID),as.character(MF$outcome),sep='*'),'')


phwas.plot<- ggplot(MF,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') +  ylim(0,30)+
  labs(title='',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=MF$label),size=4,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(AS$locus_name))), values=c("darkblue"))
p2<-phwas.plot
p2
### noonan syndrome
gene<-c('PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B')

NS<- single_snp %>% filter(locus_name%in%gene )
# Count the number of colors we'll need for each bar
ncol = table(NS$locus_name)


phewas.plot<- ggplot(NS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red')  + scale_fill_manual(values=pal) +  guides(fill=FALSE) + ylim(0,60)+
  #labs(title='Association between SNPs Noonan, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
  scale_color_manual(breaks = c(levels(as.factor(NS$locus_name))), values=c("darkblue",'darkgreen',"brown4",
                                                                            'darkred','purple','orange4','green',
                                                                            'brown','aquamarine4','coral4','chocolate1','blue4',
                                                                            'blue','dark grey','blueviolet','yellow'))

p3<-phewas.plot +theme(legend.position='top',legend.title = element_text(size=20),legend.text = element_text(size=18))
p3


### Digeorge syndrome
gene<-c('PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B','FBN1',"JAG1",'NOTCH2')

DS<- single_snp %>% filter(!(locus_name%in%gene) & `#CHROM`==22)
# Count the number of colors we'll need for each bar
ncol = table(DS$locus_name)

###plot 


phewas.plot<- ggplot(DS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red')  +  ylim(0,15)+
  #labs(title='Association between SNPs Noonan, and HPO terms',x='HPO terms',y='-log10(pvalue)',fill="Gene") +
  #geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + 
  scale_color_manual(breaks = c(levels(as.factor(DS$locus_name))), values=c("darkblue",'darkgreen',"brown4",
                                                                            'darkred','purple','orange4','green'))

p4<-phewas.plot +theme(legend.position='top',legend.title = element_text(size=20),legend.text = element_text(size=18))

p4

grid.arrange(p1,p2,p3,p4)


## single snps Phewas
snp<- single_snp %>% filter(ID=='12:112883476_G_A'| ID=='14:50655357_G_C')

phwas.plot<- ggplot(snp,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=ID),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') +  ylim(0,60)+
  #labs(title='Associati, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=MF$label),size=4,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(AS$locus_name))), values=c("brown4",'darkgreen'))

phwas.plot +theme(legend.position='top',legend.title = element_text(size=20),legend.text = element_text(size=18))

snp12450655357_G_C<- single_snp %>% filter(ID=='14:50655357_G_C')

phwas.plot<- ggplot(snp12450655357_G_C,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') +  ylim(0,60)+
  labs(title='Association between 14:50655357_G_C (SOS2), and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  #geom_text(aes(label=MF$label),size=4,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(AS$locus_name))), values=c("brown4"))

phwas.plot

##### CVN results ######

setwd('/Users/catherine.t/Documents/PheWas_syndromic_diseases/results_CNV/')


#### import PRS 
results.cnv<-fread('results.target.cnv.txt',header = T,fill=TRUE)


results.cnv$label<-ifelse(results.cnv$P < 5e-9,as.character(results.cnv$Outcome),'')

phwas.plot<- ggplot(results.cnv,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=as.factor(CHROM)),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Plot pvalue association betwen whole genome CNV and HPO terms',x='chrom',y='-log10(pvalue)') +
  #geom_text(aes(label=results.cnv$label),size=3,vjust=-1.0,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

max<-results.cnv %>% group_by(CHROM) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

p1<-phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = max$CHROM)

#### plot with HPO terms
phwas.plot2<- ggplot(results.cnv,aes(x=Outcome, y=-log10(P))) + geom_point(aes(col=as.factor(CHROM)),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + 
  labs(title='Plot pvalue association betwen whole genome CNV and HPO terms',x='HPO',y='-log10(pvalue)') +
  #geom_text(aes(label=results.cnv$label),size=3,vjust=-1.0,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot2

grid.arrange(p1,phwas.plot2)

#####
results.gene<-fread('results.target.genecnv.txt',header = T,fill=TRUE)


results.gene$label<-ifelse(results.gene$P < 1e-03,as.character(results.gene$Outcome),'')

phwas.plot<- ggplot(results.cnv,aes(x=gene, y=-log10(P))) + geom_point(aes(col=as.factor(gene)),size=2,show.legend = T) +
  geom_hline(yintercept=3,colour='blue') + 
  labs(title='Plot pvalue association betwen whole genome CNV and HPO terms',x='chrom',y='-log10(pvalue)') +
  geom_text(aes(label=results.gene$label),size=3,vjust=-1.0,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

#### plot with HPO terms
phwas.plot2<- ggplot(results.gene,aes(x=Outcome, y=-log10(P))) + geom_point(aes(col=as.factor(gene)),size=2,show.legend = T) +
  geom_hline(yintercept=3,colour='red') + 
  labs(title='Plot pvalue association betwen whole genome CNV and HPO terms',x='HPO',y='-log10(pvalue)') +
  geom_text(aes(label=results.gene$label),size=3,vjust=-1.0,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot2

grid.arrange(phwas.plot,phwas.plot2)
