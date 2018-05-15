########### Prereqs ###########
options(stringsAsFactors = FALSE,scipen = 600)
library(tidyverse)
library(gridExtra)

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
makeWideFun<-function(data){
  peptidesW<-data %>% select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
    gather(key="rep",value="Intensity",-Name)
  peptidesW<-pheno %>% left_join(groups) %>% left_join(peptidesW,by=c("newName"="rep"))
  peptidesW$GroupTime<-paste(peptidesW$Group,peptidesW$timept,sep=".")
  peptidesW$GroupTime<-factor(peptidesW$GroupTime,
                              levels=c("sCAD.FU","sCAD.T0","Type 1.FU","Type 1.T0","Type 2.FU","Type 2.T0",
                                       "Indeterminate.FU","Indeterminate.T0"))
  
  peptidesW<-peptidesW %>% arrange(GroupTime,ptid,replicate)
  temp1<-peptidesW %>% select(GroupTime,ptid,replicate) %>% unique()
  temp1$uID<-as.factor(1L:nrow(temp1))
  peptidesW<-temp1 %>% left_join(peptidesW)
  return(peptidesW)
}
peptidesW<-makeWideFun(peptides)

########### Plots ###########
png(file="plots/peptideNoNorm.png",height=5,width=10,units="in",res=300)
p0<-ggplot(data=peptidesW,aes(x=uID,group=newName,color=GroupTime,y=log2(Intensity)))+
         geom_boxplot()+theme_bw()+xlab("")+theme(axis.text.x=element_blank())
show(p0)
dev.off()

########### Normalization ###########
m0<-as.matrix(peptides[,grepl("rep",names(peptides))])

# Columnwise total intensity normalization:
cSums<-apply(m0,2,sum)
cFac<-cSums/mean(cSums)
m1<-m0
for(i in 1:ncol(m1)) m1[,i]<-m1[,i]/cFac[i]
peptides1<-peptides
peptides1[,grepl("rep",names(peptides1))]<-m1
peptidesW1<-makeWideFun(peptides1)

png(file="plots/peptideColNorm.png",height=5,width=10,units="in",res=300)
p1<-ggplot(data=peptidesW1,aes(x=uID,group=newName,color=GroupTime,y=log2(Intensity)))+
  geom_boxplot()+theme_bw()+xlab("")+theme(axis.text.x=element_blank())
show(p1)
dev.off()

png(file="plots/peptideNoneVColNorm.png",height=10,width=10,units="in",res=300)
grid.arrange(p0,p1,nrow=2)
dev.off()