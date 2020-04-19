
R
library(data.table)
library(dplyr)


### binary phenotypes
path="/scratch/users/tcheandj/revision_phewas_plos_genetic/"

files<-list.files(path=path,'*.glm.linear',full.names=T)

gwas<- list.files(path=path,'*.glm.linear')

gwas.names<- as.list(sapply(strsplit(gwas,"\\."),"[",2))

for (i in 1:length(gwas.names)){

files.gwas[[i]]<- files[files%in%files[grep(gwas.names[[i]],files)]]

gwas_results<-lapply(files.gwas[[i]],fread,header=T,fill=TRUE)

gwas_results<- rbind(data.frame,gwas_results)
gwas_results<- Reduce(function(...) merge(...,all = TRUE), gwas_results)

names(gwas_results)<- c("''")

fwrite(gwas_results,paste0(path, "result.gwas.daniala.",gwas.names[[i]],".txt"),sep='\t')
 }
