#### exploration of my phenotypes variable
library(dplyr)
library(data.table)
library(ggplot2)
library(gtable)
library(grid)
library(lattice)
library(qqman)

phenotype<-fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS/phenotype.txt')

str(phenotype)

meeas<- phenotype %>% select(HP0045081:HP0000002) %>% 
     mutate(HP0000002=ifelse(HP0000002 < median(HP0000002, na.rm = T)-4*sd(HP0000002, na.rm = T) 
                             | HP0000002 > median(HP0000002, na.rm = T)+4*sd(HP0000002, na.rm = T), NA,HP0000002))
                             
mas<- phenotype %>% select(app13721, HP0045081:HP0000002) %>% 
  mutate_at(vars(-app13721),funs(ifelse(. < median(., na.rm = T)-4*sd(., na.rm = T) 
                          | . > median(., na.rm = T)+4*sd(., na.rm = T), NA,.)))
hist(mas$HP0000002)

## plot for continue phenotype
summary(phenotype$HP0001518)
hist(phenotype$HP0000002,freq=T)
summary(phenotype$HP0000002)
sd(phenotype$HP0000002, na.rm = T)

plot(phenotype$HP0000002)

par(mfrow=c(3,2))
boxplot(x=phenotype$HP0000002)
hist(x=phenotype$HP0000002, breaks=600,xlim=c(0.46,0.60))

boxplot(x=phenotype$HP0001518)
hist(x=phenotype$HP0001518, breaks=600)
median(phenotype$HP0001518, na.rm = T)
sd(phenotype$HP0001518, na.rm = T)

boxplot(x=phenotype$HP0005117)
hist(x=phenotype$HP0005117, breaks=600)


boxplot(x=phenotype$HP0004421)
hist(x=phenotype$HP0004421, breaks=600)

boxplot(x=phenotype$HP0000120)
hist(x=phenotype$HP0000120, breaks=600)

par(mfrow=c(1,2))
boxplot(x=phenotype$HP0012594)
hist(x=phenotype$HP0012594, breaks=600,xlim=c(0,200))

boxplot(x=phenotype$HP0003081)
hist(x=phenotype$HP0003081, breaks=600)

### density plot 

p1 <- ggplot(phenotype, aes(x = HP0001518, colour = as.factor(sex)))+  geom_line(stat="density") + 
  labs(title='Small gestational age',  x='(HP0001518)') 
summary(phenotype$HP0001518)

p2 <- ggplot(phenotype, aes(x = HP0000002, colour = as.factor(sex))) + geom_line(stat="density") + 
  labs(title='Sitting/standing height ratio',  x='(HP0000002)') 


p3 <- ggplot(phenotype, aes(x = HP0005117, colour = as.factor(sex))) + geom_line(stat="density") + 
  labs(title='Diastolic blood pressure',  x=' (HP0005117)')

p4 <- ggplot(phenotype, aes(x = HP0004421 , colour = as.factor(sex))) + geom_line(stat="density") + 
  labs(title='Systolic blood pressure',  x='(HP0004421)')

p5 <- ggplot(phenotype, aes(x = HP0000120, colour = as.factor(sex))) + geom_line(stat="density") + 
  labs(title='Urine microlbumin level', x='(HP0000120)')

p6 <- ggplot(phenotype, aes(x = HP0012594, colour = as.factor(sex))) + geom_line(stat="density") + 
  labs(title='Urine creatinine level', x='(HP0012594)')

p7 <- ggplot(phenotype, aes(x = HP0003081, colour = as.factor(sex))) + geom_line(stat="density") +
   labs(title='Urine potassium level', x='(HP0003081)')
# Move to a new page
grid.newpage()

# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 3)))
# A helper function to define a region on the layout
vplayout <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))
print(p3, vp = vplayout(3, 1))
print(p4, vp = vplayout(1, 2))
print(p5, vp = vplayout(2, 2))
print(p6, vp = vplayout(3, 2))
print(p7, vp = vplayout(1, 3))

### check plot result chr15 and HP000002
HP000002<- fread('/Users/catherine.t/Documents/PheWas_syndromic_diseases/select_hpoS/result.HP000002.txt')

HP000002<- HP000002 %>% 
  mutate_at(.funs=funs(as.numeric), vars(CHROM,POS,P)) %>%
  filter(!is.na(P))

manhattan(HP000002,chr='CHROM',bp='POS',P='P',snp="ID")


plot(x=HP000002$POS,y=-log10(HP000002$P))
qq(HP000002$P)
ggplot(data=HP000002, aes(x=POS, y=len, group=1)) +
  geom_line()





## change the level to have level as on the plot 
single_snp<- single_snp[order(single_snp$CHROM, single_snp$POS,decreasing = F),]
single_snp$rownames<-rep(1:length(single_snp$CHROM))

###plot for 
phwas.plot<- ggplot(single_snp,aes(x=rownames, y=-log10(P))) + geom_point(aes(col=CHROM),size=2,show.legend = F) +
  geom_hline(yintercept=8,colour='red') + ylim(0,89)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_snp$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


single_snp$label<-ifelse(single_snp$P<10e-10,paste0(as.character(single_snp$ID),'*',
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
phwas.plot<- ggplot(single_snp,aes(x=locus_name, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  # geom_text(aes(label=single_snp$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot

### change the scale of my plot 
phwas.plot + scale_x_continuous( breaks = as.vector(max$mean), label = as.vector(max$locus_name))


### plot for HP0000002 phenotype only

### manhattan plot with gene as color
phwas.plot<- ggplot(single_snp[single_snp$outcome=='HP0000002',],aes(x=locus_name, y=-log10(P))) + geom_point(aes(col=locus_name),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + ylim(0,26)+
  labs(title='Plot pvalue association betwen snp in selected gene and HPO terms',x='gene',y='-log10(pvalue)') +
  #geom_text(aes(label=single_snp$label),size=3,vjust=-1,hjust=0.70) +
  theme(axis.text.x=element_text(angle=-90,vjust=0,hjust=0),
        panel.background = element_rect(fill = "white", colour = "grey50"))
phwas.plot


### manhattan plot with gene as color
### manhattan plot with gene as color
phwas.plot<- ggplot(single_snp[single_snp$locus_name=='FBN1' & single_snp$outcome=='HP0000002',],aes(x=POS, y=-log10(P))) + geom_point(aes(col=outcome),size=2,show.legend = T) +
  geom_hline(yintercept=8,colour='red') + ylim(0,26)+
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
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", 
     ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,28), las=1, xaxs="i", yaxs="i", bty="l",main="")
points(lexp, lobs, pch=25, cex=.4, bg="blue",col='blue') 
title(main=list('QQplot PheWas all snp in all gene with select HPO terms ',cex=1,col='black'))