### description study pop phewas ###
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)


single.snp<-fread('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/phewas_SD.final.results.txt')

### remove SNP with inf OR or p
single.snp <- single.snp %>% filter(abs(BETA)<=10,SE<10,SE>0,abs(BETA)>1e-4)

### add label and some missing info
single_snp<- single.snp %>% group_by(outcome) %>% 
  mutate(valeur=ifelse(is.na(ALT_CTRL_CT),'Continue','Categorical'),lower_ci=BETA-1.96*SE,upper_ci=BETA+1.96*SE, MAF=ALT_FREQ) %>% 
  mutate(pheno_text=ifelse(outcome=='HP0000146','ovarian cyst',pheno_text)) %>% 
  mutate(pheno_text=ifelse(outcome=='HP0000486','Strabism',pheno_text)) %>% 
  mutate(pheno_text=ifelse(outcome=='HP0000077','abnormality of the kidney',pheno_text)) %>%
  mutate(label=paste0(pheno_text,' (',outcome,')')) %>%  
  mutate(label2=ifelse(P==min(P) & min(P)<=3.2e-07, paste(rsid,outcome,sep=','),'')) %>% ungroup()

###  create a random set of color corresponding to the 
### add label on the x axis
### get pleiotropie SNPs --> defined by SNPs displaying association with at least 2 phenotypes

pleio.snp<- single_snp %>% group_by(ID) %>%  filter(P<3.2e-07) %>% ungroup() %>% filter(duplicated(ID)) %>% select(ID)

pleio.snp<- unique(pleio.snp$ID)
#### create a pleiotropie SNP data filter(ID%in%pleio.snp) 
gene<-data.frame(locus_name=as.factor(c("JAG1",'NOTCH2','FBN1','PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
                               'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B', "DGCR6" , "PRODH", "DGCR2" ,"TBX1","COMT","DGCR8" ,"LZTR1",
                               "HIC2"))
                 ,number=1:26)
pleio_snp<-merge(single_snp,gene)
pleio_snp$locus_name<- factor(pleio_snp$locus_name,levels=c("JAG1",'NOTCH2','FBN1','PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
                                                            'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B', "DGCR6" , "PRODH", "DGCR2" ,"TBX1","COMT","DGCR8" ,"LZTR1",
                                                            "HIC2"))
## create a numeric vector for ooutcome   
outcome<-data.frame(outcome=unique(pleio_snp$outcome),val=1:67)
pleio_snp<-merge(pleio_snp,outcome)
str(pleio_snp)

pleio_snp <- pleio_snp %>% arrange(number,outcome) %>% filter(!is.na(P)) %>% 
  mutate(n.row=1:nrow(single_snp))  

t<-data.frame(table(pleio_snp$locus_name,pleio_snp$outcome))
t<- t %>%  mutate(colum=row.names(t)) %>% rename(locus_name='Var1',outcome='Var2')
pleio_snp<-merge(pleio_snp,t)

str(pleio_snp)

### get color and select those that represewnt my gene
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col2=col_vector[c(1:3,5:8,10:11,47:51,33:36,66:69,22:25)]
pie(rep(1,26), col=col2)

names(col2)<- gene$locus_name

### get HPO term with significant association 
pleio.hpo<- single_snp %>% group_by(outcome) %>%  filter(P<3.2e-07)  %>% select(outcome)
pheno<-unique(pleio.hpo$outcome)
all.pheno<-unique(single_snp$outcome)
### create the shape for uniq phenotype
shape.pheno<- c(1,20,42,43,8,18,10,14,3,4,rep(111,57))
all.pheno2<- c(pheno,all.pheno[!all.pheno%in%pheno])

names(shape.pheno) <- all.pheno2

### Screate a new p value for y axis in order to get discountinues axis
threshold<-30
pos <- c(0,5,10,20,30,50, 60)
# we apply the same transform to `pos` to get `new_pos`
new_pos <- ifelse(pos > threshold, 0.20 * pos, pos - threshold * 0.80)

pleio_snp$new_log10.p <- ifelse(-log10(pleio_snp$P) > threshold, 0.20 * -log10(pleio_snp$P), 
                                        -log10(pleio_snp$P) - threshold * 0.80)


str(pleio_snp)
pleio_snp$colum<-as.factor(pleio_snp$colum)
### create a continue scale for my x axis
str(gene)
#### because I want my axis be plot in the middle of the gene band I will add 900 the count of gene
gene$ncol<- as.character(gene$number+1000)

### select the value of those number in my colum data use as x axis and because t
breaks<- pleio_snp %>%filter(colum%in%gene$ncol) %>% select(locus_name,colum) %>% filter(!duplicated(colum))

## now that I have the exact number of value needed, I will select them
### plot results by SNP with the label on the lower pvalue 
setwd('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/paper phewas/figure')

png('plot.allsnp.allphenotype.png',width = 1020, height = 880)
phwas.plot<- ggplot(pleio_snp,aes(x=reorder(colum,n.row), y=new_log10.p)) + 
    geom_point(aes(col=locus_name,shape=outcome),size=3,show.legend=F) +
     geom_hline(yintercept=6.5-threshold * 0.80,colour='brown') + 
   geom_text(aes(label=pleio_snp$label2),size=4,vjust=-1.2,hjust=0.1) +
     labs(title='',x='gene',y='-log10(pvalue)') +
   scale_y_continuous(breaks=new_pos,label=pos,expand = c(0,0),limits=c(-24,13) ) +
   scale_x_discrete(breaks = as.vector(breaks$colum), 
                    labels = as.vector(breaks$locus_name)) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
                    panel.background = element_blank(), axis.text.x=element_text(angle=-90,vjust=0,hjust=0,size=14)) + 
   scale_color_manual(values=col2) + scale_shape_manual(values=shape.pheno,guide = FALSE)
 
phwas.plot   
dev.off()



phwas.plot<- ggplot(single_snp,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=6.5,colour='red') + 
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=single_snp$label2),size=4,vjust=-1,hjust=0.5) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') 
phwas.plot +theme(legend.position='bottom',
                  legend.title = element_text(size=12),
                  label.theme=element_text(angle=90),
                  legend.text = element_text(size=12))


### table of the most significant result
top.snp<- single_snp %>% filter(label2!='') %>%  select(ID,rsid)

tabl.top.snp<- single_snp %>% filter(rsid%in%top.snp$rsid,P<5e-04) %>% select(outcome:all_traits)
## save the tablke
fwrite(tabl.top.snp,'/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/table.top.snp.phewas.txt',sep='\t')
### marfan syndrome

MF<- single_snp %>% filter(locus_name=='FBN1' )

MF <- MF %>% group_by(outcome) %>% mutate(rsid=ifelse(is.na(rsid), ID,rsid)) %>% 
  mutate(label=ifelse(P==min(P) & min(P)< 3.16e-7,rsid,"")) %>% ungroup()

png('phewasplot.MS.png',width = 920, height = 680)
phwas.plot<- ggplot(MF,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = F) +
  geom_hline(yintercept=6.5,colour='red') +  ylim(0,30)+
  labs(title='',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=MF$label),size=4,vjust=-1.2,hjust=0.1) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(MF$locus_name))), values=c("darkblue"))
phwas.plot
dev.off()


### noonan syndrome #####
gene<-c('PTPN11', 'SOS1','SOS2', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B')

NS<- single_snp %>% filter(locus_name%in%gene )

# Count the number of colors we'll need for each bar
ncol = table(NS$locus_name)
NS <- NS %>% group_by(outcome) %>% mutate(rsid=ifelse(is.na(rsid), ID,rsid)) %>% 
  mutate(label=ifelse(P==min(P) & min(P)< 3.16e-7,paste0(rsid,'(',locus_name,")"),"")) %>% ungroup()

png('phewasplot.NS.png',width = 920, height = 680)
phewas.plot<- ggplot(NS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=6.5,colour='red')  +  guides(fill=FALSE) + ylim(0,60)+
  labs(title='Association between SNPs Noonan, and HPO terms',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=NS$label),size=3,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(NS$locus_name))), values=c("darkblue",'darkgreen',"brown4",
                                                                            'darkred','purple','orange4','green',
                                                                            'brown','aquamarine4','coral4','chocolate1','blue4',
                                                                            'blue','dark grey','blueviolet','yellow'))

phewas.plot +theme(legend.position='top',legend.title = element_text(size=12),legend.text = element_text(size=12))

dev.off()

### plot Alagille Syndrome #####
AS<- single_snp %>% filter(locus_name=='JAG1' |locus_name=='NOTCH2' )

AS <- AS %>% group_by(outcome) %>% mutate(rsid=ifelse(is.na(rsid), ID,rsid)) %>% 
  mutate(label=ifelse(P==min(P) & min(P)< 3.16e-7,paste0(rsid,'(',locus_name,")"),"")) %>% ungroup()

png('phewasplot.AS.png',width = 920, height = 680)
phwas.plot<- ggplot(AS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=6.5,colour='red') + ylim(0,20)+
  labs(title='',x='HPO terms',y='-log10(pvalue)') +
  geom_text(aes(label=AS$label),size=4,vjust=-1.2,hjust=0.80) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(AS$locus_name))), values=c("Brown", "blue"))
phwas.plot +theme(legend.position='top',legend.title = element_text(size=12),legend.text = element_text(size=12))

dev.off()



### DiGeorge syndrome
gene<-c('PTPN11', 'SOS1', 'RAF1', 'KRAS', 'RIT1', 'BRAF', 'A2ML1', 'RRAS', 
        'LZTR1', 'NRAS', 'RASA2','CBL', 'SHOC2', 'MAP2K1', 'KAT6B','FBN1',"JAG1",'NOTCH2')

DS<- single_snp %>% filter(!(locus_name%in%gene) & `#CHROM`==22)
# Count the number of colors we'll need for each bar
ncol = table(DS$locus_name)

DS <- DS %>% group_by(outcome) %>% mutate(rsid=ifelse(is.na(rsid), ID,rsid)) %>% 
  mutate(label=ifelse(P==min(P) & min(P)< 3.16e-7,paste0(rsid,'(',locus_name,")"),"")) %>% ungroup()

###plot 

png('phewasplot.DS.png',width = 920, height = 680)
phewas.plot<- ggplot(DS,aes(x=outcome, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=6.5,colour='red')  +  ylim(0,15)+
  labs(title='',x='HPO terms',y='-log10(pvalue)',fill="Gene") +
  geom_text(aes(label=DS$label),size=3,vjust=-1,hjust=0) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50")) + labs(colour='Gene name') +
  scale_color_manual(breaks = c(levels(as.factor(DS$locus_name))), values=c("darkblue",'darkgreen',"brown4",
                                                                         'darkred','purple','orange4','green'))
phewas.plot +theme(legend.position='top',legend.title = element_text(size=12),legend.text = element_text(size=12))

dev.off()


### plot gene results for noonan syndrom
NS.generesult<-fread('/Users/catherine/Documents/Documents/PheWas_syndromic_diseases/result_phewas_121218_with_standa_vars/NS.generesult.txt',sep='\t')

NS.generesult.cadd<- NS.generesult %>%  arrange(locus_name,outcome) %>%  
   filter(!is.na(P.weightCADD)) %>% mutate(rownames=1:nrow(NS.generesult),ax=paste(outcome,locus_name,sep='_')) 

col=c("darkblue",'darkgreen',"brown4",'red3','purple','orange','blue','slateblue4',
      'aquamarine4','coral','black','yellowgreen',
      'yellow','dark grey','goldenrod','magenta')

## assign colors top my set
names(col) <- levels(as.factor(NS.generesult.cadd$locus_name))

shapes<-c(3:13,22:24,7:9,22:24,7:9,11,12,8,1)
names(shapes) <- levels(as.factor(NS.generesult.cadd$outcome))

## create a new threshold for the y axis 
threshold<-10
## create a new p so that the value <20 have 90 percent of the data

NS.generesult.cadd$new_log10.p <- ifelse(-log10(NS.generesult.cadd$P.weightCADD) > threshold, 0.1 * -log10(NS.generesult.cadd$P.weightCADD), 
                      -log10(NS.generesult.cadd$P.weightCADD) - threshold * 0.9)

pos <- c(min(-log10(NS.generesult.cadd$P.weightCADD)),4,8,12,20,40, 80)
## we apply the same transform to `pos` to get `new_pos`
new_pos <- ifelse(pos > threshold, 0.1 * pos, pos - threshold * 0.9)

### add label on the x axis
max<-NS.generesult.cadd %>% group_by(locus_name) %>% summarise(mean = mean(rownames))
max<- max[order(max$mean,decreasing = F),]

nsplot<- qplot(x=rownames, y=new_log10.p, data=NS.generesult.cadd,colour=locus_name, margins = T,show.legend = F,
                xlab="", ylab="-log10(p_value)")+ geom_point(aes(shape = outcome),size=3) + 
  geom_hline(aes(yintercept=(3.1-threshold * 0.9)), colour="#990000", linetype="dashed") +
  scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name)) + 
  scale_y_continuous(breaks=new_pos,label=pos) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(), panel.border = element_rect(fill = NA,colour = "black"),
        panel.background = element_blank(), axis.text.x=element_text(angle=-90,vjust=0,hjust=0,size=14)) + 
  scale_color_manual(values=col) +  scale_shape_manual(values=shapes) +theme(legend.position='bottom',legend.title = element_text(size=10),
                                                                             legend.text = element_text(size=12),legend.direction = "vertical") + 
  guides(colour = guide_legend(nrow = 2))


