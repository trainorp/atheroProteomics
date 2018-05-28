########### Prereqs ###########
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)
library(gridExtra)

setwd("~/gdrive/AthroProteomics/data")
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
names(proteins)[match(pheno$oldName,names(proteins))]<-pheno$newName
pheno$uSamp<-paste(pheno$ptid,pheno$timept,sep="_")

# Phenotype data:
groups<-read.csv(file="~/gdrive/Athro/groups_20180515.csv")
groups$ptid<-as.character(groups$ptid)
groups$Group<-"sCAD"
groups$Group[!is.na(groups$MIGroup)]<-groups$MIGroup[!is.na(groups$MIGroup)]

pheno<-pheno %>% left_join(groups)

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

# Export peptide sequences for protein query:
pepSeqs<-peptides$Name
pepSeqs<-gsub("(\\[.*?\\])","",pepSeqs)
pepSeqsStr<-paste(pepSeqs,collapse=",")
# write.table(pepSeqsStr,file="pepSeqsStr.txt",quote=FALSE,row.names=FALSE)

# Export other peptide annotation:
pepAnno<-peptides %>% select(Name,Quality.Score,ParentProtein.FullName,Use.For.Quant,
                             Percent.Files.With.Good.Quant)
# save(pepAnno,file="pepAnno.RData")

setwd("~/gdrive/AthroProteomics/")

########### Normalization ###########
m0<-as.matrix(peptides[,grepl("rep",names(peptides))])

# Minimum value imputation:
mins<-apply(m0,2,function(x) min(x[x>0]))
for(i in 1:ncol(m0)){
  m0[,i][m0[,i]<1e-6]<-mins[i]
}

# Log-transform:
m00<-m0
peptides00<-peptides
peptides00[,grepl("rep",names(peptides00))]<-m00
m0<-log2(m0)
peptides[,grepl("rep",names(peptides))]<-m0

peptidesL<-makeWideFun(peptides)
# png(file="plots/peptideNoNorm.png",height=5,width=10,units="in",res=300)
p0<-ggplot(data=peptidesL,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p0)
# dev.off()

# Columnwise total intensity normalization:
cSums<-apply(m00,2,sum)
cFac<-cSums/mean(cSums)
m1<-m00
for(i in 1:ncol(m1)) m1[,i]<-m1[,i]/cFac[i]
m1b<-m1
# Median normalization
cMed<-apply(m1,2,median)
cFac2<-cMed/mean(cMed)
for(i in 1:ncol(m1b)) m1b[,i]<-m1b[,i]/cFac2[i]

peptides1b<-peptides1<-peptides00
peptides1[,grepl("rep",names(peptides1))]<-log2(m1)
peptides1b[,grepl("rep",names(peptides1b))]<-log2(m1b)

peptidesL1<-makeWideFun(peptides1)
peptidesL1b<-makeWideFun(peptides1b)
# png(file="plots/peptideColNorm.png",height=5,width=10,units="in",res=300)
p1<-ggplot(data=peptidesL1,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p1)
# dev.off()

# png(file="plots/peptideColMedNorm.png",height=5,width=10,units="in",res=300)
p1b<-ggplot(data=peptidesL1b,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p1b)
# dev.off()

# png(file="plots/peptideNoneVColNorm.png",height=10,width=10,units="in",res=300)
grid.arrange(p0,p1,nrow=2)
# dev.off()

# Quantile normalization:
m2<-preprocessCore::normalize.quantiles(m0)
peptides2<-peptides
peptides2[,grepl("rep",names(peptides2))]<-m2
peptidesL2<-makeWideFun(peptides2)
# png(file="plots/peptideQuantNorm.png",height=5,width=10,units="in",res=300)
p2<-ggplot(data=peptidesL2,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p2)
# dev.off()

# Cyclic loess:
m3<-limma::normalizeCyclicLoess(m0)
peptides3<-peptides
peptides3[,grepl("rep",names(peptides3))]<-m3
peptidesL3<-makeWideFun(peptides3)
# png(file="plots/peptideFastCyclicLoess.png",height=5,width=10,units="in",res=300)
p3<-ggplot(data=peptidesL3,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p3)
# dev.off()

# Cyclic loess bGal only:
bGalWeights<-as.numeric(grepl("P00722",peptides$Parent.Protein) & 
         peptides$Percent.Files.With.Good.Quant>.99)
m4<-limma::normalizeCyclicLoess(m0,weights = bGalWeights)
peptides4<-peptides
peptides4[,grepl("rep",names(peptides4))]<-m4
peptidesL4<-makeWideFun(peptides4)
# png(file="plots/peptideFastCyclicLoessbGalOnly.png",height=5,width=10,units="in",res=300)
p4<-ggplot(data=peptidesL4,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p4)
# dev.off()

# Cyclic loess half & half bGal only:
bGalWeights2<-bGalWeights+.2
bGalWeights2<-ifelse(bGalWeights2>1,1.0,bGalWeights2)
m5<-limma::normalizeCyclicLoess(m0,weights = bGalWeights2)
peptides5<-peptides
peptides5[,grepl("rep",names(peptides5))]<-m5
peptidesL5<-makeWideFun(peptides5)
# png(file="plots/peptideFastCyclicLoessbGalHeavy.png",height=5,width=10,units="in",res=300)
p5<-ggplot(data=peptidesL5,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p5)
# dev.off()

# Cyclic loess:
m6<-limma::normalizeMedianAbsValues(m0)
peptides6<-peptides
peptides6[,grepl("rep",names(peptides6))]<-m6
peptidesL6<-makeWideFun(peptides6)
# png(file="plots/peptideMAD.png",height=5,width=10,units="in",res=300)
p6<-ggplot(data=peptidesL6,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p6)
# dev.off()

# Beta gal:
bGalFun<-function(data,method){
  # betaGal data:
  bGal<-data %>% filter(grepl("P00722",ParentProtein.FullName) & 
                          Percent.Files.With.Good.Quant>.99)
  set.seed(3)
  bGalShort<-sample(bGal$Name,8)
  bGalL<-bGal %>% select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
    gather(key="rep",value="Intensity",-Name)
  bGalLNames<-bGalL %>% select(Name) %>% unique()
  bGalLNames$id<-as.factor(1:nrow(bGalLNames))
  bGalL<-bGalLNames %>% left_join(bGalL)
  
  # CVs:
  cv<-bGalL %>% group_by(Name) %>% summarize(cv=sd(2**Intensity)/mean(2**Intensity)*100)
  names(cv)[names(cv)=="cv"]<-paste0(method,"CV")
  
  # Plots:
  fName<-paste0("plots/bGalPeptides_",method,"_Full.png")
  # png(filename=fName,height=5,width=8,units="in",res=300)
  bGalp<-ggplot(bGalL,aes(x=rep,y=Intensity,group=Name,color=id))+
    geom_line()+ylim(23,32)+theme_bw()+xlab("Replicate")+ylab("Intensity (log scale)")
  print(bGalp)
  # dev.off()
  fName2<-paste0("plots/bGalPeptides_",method,"_Samp.png")
  # png(filename=fName2,height=5,width=8,units="in",res=300)
  bGalp2<-ggplot(bGalL %>% filter(Name %in% bGalShort),
                    aes(x=rep,y=Intensity,group=Name,color=id))+
    geom_line()+ylim(26,31.75)+theme_bw()+xlab("Replicate")+ylab("Intensity (log scale)")
  print(bGalp2)
  # dev.off()
  
  return(cv)
}
bGalCVs<-cbind(
  bGalFun(data=peptides,method="none"),
  bGalFun(data=peptides1,method="columnTI")[,-1],
  bGalFun(data=peptides1b,method="columnTIMed")[,-1],
  bGalFun(peptides2,method="Quant")[,-1],
  bGalFun(peptides3,method="FCycLoess")[,-1],
  bGalFun(peptides4,method="FCycLoessbGal")[,-1],
  bGalFun(peptides5,method="FCycLoessbGalLight")[,-1],
  bGalFun(peptides6,method="MAD")[,-1]
)
# write.csv(bGalCVs,file="bGalCVs.csv",row.names=FALSE)

########### Un log-transform ###########
peptides1[,grepl("rep_",names(peptides1))]<-
  2**peptides1[,grepl("rep_",names(peptides1))]
peptides2[,grepl("rep_",names(peptides2))]<-
  2**peptides2[,grepl("rep_",names(peptides2))]
peptides3[,grepl("rep_",names(peptides3))]<-
  2**peptides3[,grepl("rep_",names(peptides3))]
peptides4[,grepl("rep_",names(peptides4))]<-
  2**peptides4[,grepl("rep_",names(peptides4))]
peptides5[,grepl("rep_",names(peptides5))]<-
  2**peptides5[,grepl("rep_",names(peptides5))]

########### My beta-gal protein normalization ###########
load("~/gdrive/AthroProteomics/data/pepAnno2.RData")
pinEcoli<-pepAnno2 %>% filter(grepl("P00722",proteins) & goodQuant>.99 
                              & pinUseQuant=="Yes")
myEColi<-pepAnno2 %>% filter(grepl("P00722",proteins) & goodQuant>.99 
                             & !grepl("(\\[.*?\\])",pepSeq))

myEColiIntensPep00<-peptides00[peptides00$Name %in% myEColi$pepSeq,
                          grepl("rep_",names(peptides00))]
myEColiIntensPep1<-peptides1[peptides1$Name %in% myEColi$pepSeq,
                               grepl("rep_",names(peptides1))]
myEColiIntensPep2<-peptides2[peptides2$Name %in% myEColi$pepSeq,
                                 grepl("rep_",names(peptides2))]
myEColiIntensPep3<-peptides3[peptides3$Name %in% myEColi$pepSeq,
                                 grepl("rep_",names(peptides3))]
myEColiIntensPep4<-peptides4[peptides4$Name %in% myEColi$pepSeq,
                                 grepl("rep_",names(peptides4))]
myEColiIntensPep5<-peptides5[peptides5$Name %in% myEColi$pepSeq,
                                 grepl("rep_",names(peptides5))]

sd(apply(myEColiIntensPep00,2,sum))/mean(apply(myEColiIntensPep00,2,sum))
sd(apply(myEColiIntensPep1,2,sum))/mean(apply(myEColiIntensPep1,2,sum))
sd(apply(myEColiIntensPep2,2,sum))/mean(apply(myEColiIntensPep2,2,sum))
sd(apply(myEColiIntensPep3,2,sum))/mean(apply(myEColiIntensPep3,2,sum))
sd(apply(myEColiIntensPep4,2,sum))/mean(apply(myEColiIntensPep4,2,sum))
sd(apply(myEColiIntensPep5,2,sum))/mean(apply(myEColiIntensPep5,2,sum))

########### Combine injections ###########
combFun<-function(Names,data){
  # Output dataframe:
  df2<-data.frame(Name=Names)
  
  # Combine injections
  for(uSamp in unique(pheno$uSamp)){
    df1<-data.frame(value=apply(data[,pheno$newName[pheno$uSamp==uSamp]],1,mean))
    names(df1)[names(df1)=="value"]<-paste0("rep_",uSamp)
    df2<-cbind(df2,df1)
  }
  
  # Median scale & log-transform
  meds<-apply(df2[,grepl("rep",names(df2))],1,median)
  for(i in 1:nrow(df2[,grepl("rep",names(df2))])){
    df2[,grepl("rep",names(df2))][i,]<-log2(df2[,grepl("rep",names(df2))][i,]/meds[i])
  }
  return(df2)
}
idk<-combFun(Names=peptides1$Name,data=peptides1)

########### Peptide difference at baseline ###########
idk2<-idk %>% gather(key="rep",value="Intensity",-Name)
idk2$ptid<-str_split(idk2$rep,"_",simplify=TRUE)[,2]
idk2$timept<-str_split(idk2$rep,"_",simplify=TRUE)[,3]


########### How peptides were aggregated into proteins ###########
pep2prot<-peptides00 %>% select(Name,Parent.Protein,Use.For.Quant,rep_1,rep_2) %>% 
  filter(Use.For.Quant=="Yes") %>% 
  group_by(Parent.Protein) %>% summarize(sum(rep_1),sum(rep_2))

