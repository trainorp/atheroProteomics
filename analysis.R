########### Prereqs ###########
options(stringsAsFactors = FALSE,scipen = 600)
library(tidyverse)

setwd("~/gdrive/AthroProteomics")
peptides<-read.csv('peptide_table_20180514.csv')
proteins<-read.csv('protein_table_20180514.csv')

reg1<-regexpr("JVE.",text=names(peptides))
reg2<-regexpr(".File.Area.s.",text=names(peptides))
reg1DF<-data.frame(include=reg1>0,start=as.integer(reg1)+attr(reg1,"match.length"),
                   stop=as.integer(reg2)-1L)
names1<-substr(names(peptides),start=reg1DF$start,stop=reg1DF$stop)
names2<-str_split(names1,"\\.",simplify=TRUE)
pheno<-data.frame(oldName=names(peptides)[names1!=""],ptid=names2[names1!="",1],
                  timept=names2[names1!="",2],replicate=as.integer(names2[names1!="",3]))
pheno$newName<-paste0("rep_",1L:length(pheno$oldName))
names(peptides)[match(pheno$oldName,names(peptides))]<-pheno$newName

# Phenotype data:
groups<-read.csv(file="~/gdrive/Athro/groups_20180515.csv")
groups$ptid<-as.character(groups$ptid)
groups$Group<-"sCAD"
groups$Group[!is.na(groups$MIGroup)]<-groups$MIGroup[!is.na(groups$MIGroup)]

# Wide data:
peptidesW<-peptides %>% select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
  gather(key="rep",value="value",-Name)
peptidesW<-pheno %>% left_join(groups) %>% left_join(peptidesW,by=c("newName"="rep"))
peptidesW$GroupTime<-paste(peptidesW$Group,peptidesW$timept,sep=".")
peptidesW$GroupTime<-factor(peptidesW$GroupTime,
  levels=c("sCAD.FU","sCAD.T0","Type 1.FU","Type 1.T0","Type 2.FU","Type 2.T0",
           "Indeterminate.FU","Indeterminate.T0"))

peptidesW<-peptidesW %>% arrange(GroupTime,ptid,replicate)
temp1<-peptidesW %>% select(GroupTime,ptid,replicate) %>% unique()
temp1$uID<-as.factor(1L:nrow(temp1))
peptidesW<-temp1 %>% left_join(peptidesW)

########### Plots ###########
png(file="plots/peptideNoNorm.png",height=5,width=15,units="in",res=300)
ggplot(data=peptidesW,aes(x=uID,group=newName,color=GroupTime,y=log2(value)))+
         geom_boxplot()+theme_bw()
dev.off()

########### Normalization ###########
m0<-as.matrix(peptides[,grepl("rep",names(peptides))])

# Columnwise total intensity normalization:
cSums<-apply(m0,2,sum)
cFac<-cSums/mean(cSums)
