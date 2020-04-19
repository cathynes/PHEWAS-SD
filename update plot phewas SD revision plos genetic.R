library(lattice)
library(ggplot2)
library(data.table)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(gtable)
library(MASS)
library(gridExtra)
library(viridis)
library(grid)


################################################
#### create a plot for single SNP association
################################################
## get color for each gene
gene<-data.frame(locus_name=as.factor(c("JAG1",'NOTCH2','FBN1','PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
                                        'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B', "DGCR6" , "PRODH", "DGCR2" ,"TBX1","COMT","DGCR8" ,"LZTR1",
                                        "HIC2")),number=1:26)

### get color and select those that represewnt my gene
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,74), col=col_vector)
col2=col_vector[c(1:3,5:8,10:11,52,47:50,53,33,35,9,67:69,22:25,9)]
pie(rep(1,26), col=col2)

names(col2)<- gene$locus_name
gene$colors<-col2

### create the heatmap for the set of independent SNP
single.snp <- fread("/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/data_revision_plosgenetic/results/all.SNPs.single.level.asso.with.annotation.csv")

indep.set<-fread("/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/data_revision_plosgenetic/results/Ind.snp.set.all.pheno_with_annotation.csv")

### here I'm removing and sudset of SNPs that are not ind
indep_set<- single.snp %>% filter(ID%in%indep.set$ID, P<3.2e-07,!rsid%in%c("rs8036173",'rs9835593','rs6040076')) %>% rename(b="BETA",se="SE") %>% 
                           mutate(rsid=ifelse(ID=="7:140422661_C_G","rs1263647022",rsid),p=ifelse(P<1e-20,1e-20,P))

levels=c("JAG1",'NOTCH2','FBN1','PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
         'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B', "DGCR6" , "PRODH", "DGCR2" ,"TBX1","COMT","DGCR8")

indep_set$genes<-factor(indep_set$genes,levels=levels[indep_set$genes%in%levels])

indep_set<- indep_set %>% mutate(zscore=b/se,logp=-log10(p),snp=paste0(genes,",",rsid)) %>% 
                          mutate(zscore=ifelse(zscore>10,10,zscore))

## merge genes and ind_set to have color names
indep_set<- merge(indep_set,gene,by.x="genes",by.y="locus_name")

## i will have to ocreate a specific order that mach what i want 

color.snp.uniq<- indep_set %>% filter(!duplicated(rsid)) %>%  dplyr::select(rsid,number,colors) 

### make sure the color smatched the order of interest
colors<- color.snp.uniq$colors[order(color.snp.uniq$number,color.snp.uniq$rsid,decreasing = F)]
  


ind.plot<-ggplot(indep_set, aes(x=reorder(snp,number), y=pheno_text,color=zscore,size=logp)) + 
  geom_point(alpha=0.5) + 
  scale_color_viridis(option = "D") + 
  scale_size_continuous(breaks = c(6,8,10,12,16,18,20),labels = c(6,8,10,12,16,18,">20")) +
  theme_bw() + scale_color_gradient2(low = "darkblue",high = "darkred",mid = "grey") +
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=11,face="italic",colour = colors)) +
  theme(axis.text.y=element_text(size=12),legend.position='right',legend.title = element_text(size=8),
        legend.text = element_text(size=8),legend.direction = "vertical", 
        plot.background = element_rect(fill = "grey99"))



setwd('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/paper phewas/revision_plos_genetics/')

tiff('plot_indep_snps_from_gcta.tiff',width = 1200, height = 460)
ind.plot   
dev.off()

### import the data 
single.snp <- fread("/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/data_revision_plosgenetic/results/all.SNPs.single.level.asso.with.annotation.csv")

single.snp<- single.snp %>% filter(is.na(A1_CASE_CT) | A1_CASE_CT>=20) 


### get HPO term with significant association


pleio.hpo<- single.snp %>% group_by(pheno_text) %>%  filter(P<3.2e-07)  %>% 
  dplyr::select(pheno_text)

pheno<-unique(pleio.hpo$pheno_text)

all.pheno<-unique(single.snp$pheno_text)

### create the shape for uniq phenotype
shape.pheno<- c(0:14,16,18,rep(111,67))
all.pheno2<- c(pheno,all.pheno[!all.pheno%in%pheno])

names(shape.pheno) <- all.pheno2

### Screate a new p value for y axis in order to get discountinues axis
threshold<-30
pos <- c(0,5,10,20,30,50, 60)
# we apply the same transform to `pos` to get `new_pos`
new_pos <- ifelse(pos > threshold, 0.20 * pos, pos - threshold * 0.80)

single.snp$new_log10.p <- ifelse(-log10(single.snp$P) > threshold, 0.20 * -log10(single.snp$P), 
                                -log10(single.snp$P) - threshold * 0.80)


str(single.snp)

## create value for the X axis 
single.snp<- merge(single.snp,gene,by.x="genes", by.y="locus_name")

single.snp<- single.snp %>% arrange(number,outcome) %>% 
                           mutate(n.row=1:nrow(single.snp))

### get the median value of n.row for each genes to use as value for x axis label
breaks<- single.snp %>% group_by(genes) %>% summarise(median = median(n.row, na.rm = TRUE)) %>%
                        mutate(median=as.integer(median)) %>% arrange(median)

### create the label only for the most significant association for each gene
single.snp <- single.snp %>% group_by(genes) %>% mutate(label2=ifelse(P==min(P) & min(P)<=3.2e-07, paste(rsid,pheno_text,sep=', '),''),
                                                        median=ifelse(n.row==as.integer(median(n.row)),n.row,"")) %>% 
              mutate(label2=ifelse(ID=='11:119116748_G_A'| P<1e-20,"",label2),
                     label3=ifelse(P<1e-20 & P==min(P),paste(rsid,pheno_text,sep=', '),'')) %>% ungroup()

## now that I have the exact number of value needed, I will select them
### plot results by SNP with the label on the lower pvalue 
setwd('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/paper phewas/revision_plos_genetics/')

phwas.plot<- ggplot(single.snp,aes(x=n.row, y=new_log10.p)) + 
  geom_point(aes(col=genes,shape=pheno_text),size=2,show.legend=F) +
  geom_hline(yintercept=6.5-threshold * 0.80,colour='brown',lty = 2) + 
  geom_text(aes(label=single.snp$label2,angle=80,hjust=0,vjust=0), 
          position=position_nudge(x = 0.3, y = 0.3),size=2.7) +
  labs(title='',x='gene',y='-log10(pvalue)') +
  scale_y_continuous(breaks=new_pos,label=pos,expand = c(0,0),limits=c(-24,13)) +
  scale_x_continuous(breaks =as.vector(breaks$median), labels =as.vector(breaks$genes),expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank(), axis.text.x=element_text(angle=80,vjust=1,hjust=1,size=8,face="italic"),
        axis.text.y=element_text(size=14)) + 
  scale_color_manual(values=col2) + scale_shape_manual(values=shape.pheno)

pdf('plot.allsnp.allphenotype2.pdf',width = 12, height = 8,paper="a4r")
phwas.plot + geom_text_repel(aes(label = single.snp$label3),size=3)
dev.off()

#### get the annotation for the snp plot will plot only significan phenotype 
sign.pheno<- single.snp %>% filter(pheno_text%in%pheno)

phwas.plot2<- ggplot(sign.pheno,aes(x=n.row, y=new_log10.p)) + 
  geom_point(aes(col=genes,shape=pheno_text),size=2,show.legend=T) +
  geom_hline(yintercept=6.5-threshold * 0.80,colour='brown',lty = 2) + 
  #geom_text(aes(label=single.snp$label2,angle=80,hjust=0,vjust=0), 
    #      position=position_nudge(x = 0.3, y = 0.3),size=2.7) +
  labs(title='',x='gene',y='-log10(pvalue)') + 
  #scale_y_continuous(breaks=new_pos,label=pos,expand = c(0,0),limits=c(-24,13)) +
 # scale_x_continuous(breaks =as.vector(breaks$median), labels =as.vector(breaks$genes),expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank(), axis.text.x=element_text(angle=80,vjust=1,hjust=1,size=8,face="italic"),
        axis.text.y=element_text(size=14),legend.position = "top") + guides(colour=FALSE) +
  scale_color_manual(values=col2) + scale_shape_manual(values=shape.pheno)


pdf('plot.allsnp.allphenotype.shap.legend.pdf',width = 12, height = 8,paper="a4r")
phwas.plot2
dev.off()


##### get a 

png('plot.allsnp.allphenotype2.png',width = 1200, height = 880)
phwas.plot   
dev.off()

### create the plot for the gene level association
gene.result <- fread("/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/data_revision_plosgenetic/results/results_skat_all.with.annotation.csv")

str(gene.result)

## order by gene as in the snp result
gene<-data.frame(locus_name=as.factor(c("JAG1",'NOTCH2','FBN1','PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
                                        'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B', "DGCR6" , "PRODH", "DGCR2" ,"TBX1","COMT","DGCR8" ,"LZTR1",
                                        "HIC2"))
                 ,number=1:26)

gene.result<-merge(gene,gene.result,by.x=c("locus_name"),by.y="g")

gene.result$locus_name<- factor(gene.result$locus_name,levels=c("JAG1",'NOTCH2','FBN1','PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
                                                                'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B', "DGCR6" , "PRODH", "DGCR2" ,"TBX1","COMT","DGCR8" ,"LZTR1",
                                                                "HIC2"))

gene.result<- gene.result %>%  arrange(locus_name,outcome) %>%  
  filter(!is.na(P.weightCADD)) %>% mutate(rownames=1:1794,
                                          ax=paste(outcome,locus_name,sep='_'))  ### sort by gene to have a correct order

## create a label for the gene level association 
gene.result <- gene.result %>%  group_by(locus_name) %>% mutate(label=ifelse(P.weightCADD==min(P.weightCADD) & P.weightCADD<4e-04 | P.weightCADD<=5e-05,pheno_text,''))

### create the shape for phenotype
pheno<-unique(pleio.hpo$outcome) 

### get the set of significant phenotype
sign.pheno<- unique(gene.result$pheno_text[gene.result$P.cadd.frd<1e-02])
other.pheno<- unique(gene.result$pheno_text[gene.result$P.cadd.frd>1e-02])
shape.pheno<- c(1:15,rep(16,69))
all.pheno2<- c(sign.pheno,other.pheno)
names(shape.pheno) <- all.pheno2


### create a new y axis 

threshold <- 20

## create a new p so that the value <20 have 90 percent of the data
#### to avoid missing value and  chiti plot turn p>=1 tp 0 o.95
gene.result$P.weightCADD <- ifelse(gene.result$P.weightCADD>=1,0.95,gene.result$P.weightCADD)

gene.result$new_log10.p <- ifelse(-log10(gene.result$P.weightCADD) > threshold, 0.1 * -log10(gene.result$P.weightCADD), 
                                  -log10(gene.result$P.weightCADD) - threshold * 0.9)

pos <- c(0,4,8,12,20,40, 80)
## we apply the same transform to `pos` to get `new_pos`
new_pos <- ifelse(pos > threshold, 0.1 * pos, pos - threshold * 0.9)

### add label on the x axis
max<-gene.result %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]


#gene.result$rownames<-as.factor(gene.result$rownames)
str(gene.result)

geneplot<- ggplot(data=gene.result,aes(x=rownames, y=new_log10.p)) +
  geom_point(aes(col=locus_name,shape=pheno_text),size=3,show.legend=F) +
  labs(title='',x='gene',y='-log10(pvalue)') +
  geom_hline(aes(yintercept=3.1- threshold * 0.9), colour="#990000", linetype="dashed") +
  #geom_text(aes(label=gene.result$label,angle=45,hjust=0.3,vjust=-0.7), check_overlap = T,
  #          position=position_nudge(x = 1, y = 0.3),size=3) +
  scale_x_continuous( breaks = as.vector(max$mean), 
                      label = as.vector(max$locus_name),expand = c(0,0)) + 
  scale_y_continuous(breaks=new_pos,label=pos,expand = c(0,0),limits = c(-18,9)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank(), axis.text.x=element_text(angle=80,vjust=1,hjust=1,size=8,face="italic"),
        axis.text.y=element_text(size=14)) + 
  scale_color_manual(values=col2) +  scale_shape_manual(values=shape.pheno) +
  theme(legend.position='top',legend.title = element_text(size=10),
        legend.text = element_text(size=12),legend.direction = "vertical") + 
  guides(colour = guide_legend(nrow = 2))

#geneplot

png('plot.allgene.allphenotype.png',width = 1000, height = 680)
geneplot + geom_text_repel(aes(label = gene.result$label),size=2.6) 
dev.off()


### create a correlation plot for the significant association
sign.gene.pheno<- gene.result %>% filter(P.weightCADD<1e-03) %>% dplyr::select(outcome,locus_name)

sign.gene<- gene.result %>% filter(outcome%in%sign.gene.pheno$outcome & 
                                     locus_name%in%sign.gene.pheno$locus_name) %>%
                            mutate(logp=ifelse(P.weightCADD<1e-10,10,-log10(P.weightCADD)))


sign.gene.plot<- ggplot(sign.gene, aes(x=locus_name, y=pheno_text)) + 
  geom_tile(aes(fill =logp), colour = "grey",show.legend = T) + 
  scale_fill_gradient2(low = "grey90",mid="white",midpoint = 3,high = "darkred") +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=10,face = "italic"),
        legend.position='right',legend.title = element_text(size=10),
        legend.text = element_text(size=10)) + 
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) + 
  theme(panel.background=element_rect(fill="white", colour="black"))

sign.gene.plot

### create a correlation plot for the SNP with significat association

sign.snp<- single.snp  %>% filter(P<3.2e-07) %>% mutate(MAF=ifelse(A1_FREQ>0.5,1-A1_FREQ,A1_FREQ),
                                                        beta.maf=ifelse(A1==ALT & A1_FREQ>0.5, -1*BETA,BETA))


### correlation plot
corr.plot<- ggplot(sign.snp,aes(x=MAF, y=abs(BETA))) + 
  geom_point(aes(col=genes,shape=pheno_text),size=2,show.legend=F) +
  scale_x_continuous( breaks = c(seq(min(sign.snp$MAF),0.51,0.05)),
                      labels = c(6.7e-05,round(seq(0.05,0.5,0.05),2))) +
  labs(title='',x='MAF',y='abs(BETA)') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank(),legend.position='top') + 
  scale_color_manual(values=col2) + scale_shape_manual(values=shape.pheno)

## plot restrict to rare variant only 
corr.plot2<- ggplot(sign.snp[sign.snp$MAF<0.01,],aes(x=MAF, y=abs(BETA))) + 
  geom_point(aes(col=genes,shape=pheno_text),size=2,show.legend=F) +
  scale_x_continuous( breaks = c(seq(0.0000675,0.01,0.002)),
                      labels = c("6.7e-05","2e-03","4e-03","6e-03","8e-03")) +
  labs(title='',x='MAF',y='abs(BETA)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank()) + 
  scale_color_manual(values=col2) + scale_shape_manual(values=shape.pheno)

corr.plot3<- ggplot(sign.snp[sign.snp$MAF>=0.01,],aes(x=MAF, y=abs(BETA))) + 
  geom_point(aes(col=genes,shape=pheno_text),size=2,show.legend=F) +
  labs(title='',x='MAF',y='abs(BETA)') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank()) + 
  scale_color_manual(values=col2) + scale_shape_manual(values=shape.pheno)

pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(corr.plot2, vp = vplayout(1, 1))
print(corr.plot3, vp = vplayout(1, 2))
print(corr.plot, vp = vplayout(2, 1:2))
