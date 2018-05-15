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
  peptidesL<-data %>% select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
    gather(key="rep",value="Intensity",-Name)
  peptidesL<-pheno %>% left_join(groups) %>% left_join(peptidesL,by=c("newName"="rep"))
  peptidesL$GroupTime<-paste(peptidesL$Group,peptidesL$timept,sep=".")
  peptidesL$GroupTime<-factor(peptidesL$GroupTime,
                              levels=c("sCAD.FU","sCAD.T0","Type 1.FU","Type 1.T0","Type 2.FU","Type 2.T0",
                                       "Indeterminate.FU","Indeterminate.T0"))
  
  peptidesL<-peptidesL %>% arrange(GroupTime,ptid,replicate)
  temp1<-peptidesL %>% select(GroupTime,ptid,replicate) %>% unique()
  temp1$uID<-as.factor(1L:nrow(temp1))
  peptidesL<-temp1 %>% left_join(peptidesL)
  return(peptidesL)
}
peptidesL<-makeWideFun(peptides)

########### Plots ###########
png(file="plots/peptideNoNorm.png",height=5,width=10,units="in",res=300)
p0<-ggplot(data=peptidesL,aes(x=uID,group=newName,color=GroupTime,y=log2(Intensity)))+
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
peptidesL1<-makeWideFun(peptides1)
png(file="plots/peptideColNorm.png",height=5,width=10,units="in",res=300)
p1<-ggplot(data=peptidesL1,aes(x=uID,group=newName,color=GroupTime,y=log2(Intensity)))+
  geom_boxplot()+theme_bw()+xlab("")+theme(axis.text.x=element_blank())
show(p1)
dev.off()

png(file="plots/peptideNoneVColNorm.png",height=10,width=10,units="in",res=300)
grid.arrange(p0,p1,nrow=2)
dev.off()

# Beta gal:
bGalFun<-function(data){
  bGal<-data %>% filter(grepl("P00722",ParentProtein.FullName))
  bGalL<-bGal %>% select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
    gather(key="rep",value="Intensity",-Name)
  bGalLNames<-bGalL %>% select(Name) %>% unique()
  bGalLNames$id<-as.factor(1:nrow(bGalLNames))
  bGalL<-bGalLNames %>% left_join(bGalL)
  bGalPlot<-ggplot(bGalL,aes(x=rep,y=log2(Intensity),group=Name,color=id))+
    geom_line()+theme_bw()
  return(bGalPlot)
}

# Quantile normalization:
m2<-preprocessCore::normalize.quantiles(m0)
peptides2<-peptides
peptides2[,grepl("rep",names(peptides2))]<-m2
peptidesL2<-makeWideFun(peptides2)
png(file="plots/peptideQuantNorm.png",height=5,width=10,units="in",res=300)
p2<-ggplot(data=peptidesL2,aes(x=uID,group=newName,color=GroupTime,y=log2(Intensity)))+
  geom_boxplot()+theme_bw()+xlab("")+theme(axis.text.x=element_blank())
show(p2)
dev.off()

bGalFun(peptides)
bGalFun(peptides1)
bGalFun(peptides2)
