########### Prereqs ###########
# Begin always run:
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)
library(gridExtra)
library(emmeans)
library(multcomp)
library(corrplot)
# Ensembl library()

# End always run

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
# Remove 2010_T0 because there are 4 replicates!
peptides<-peptides[,!(colnames(peptides) %in% pheno$newName[pheno$uSamp=="2010_T0"])]
pheno<-pheno %>% filter(uSamp != "2010_T0")

# Wide data:
makeWideFun<-function(data){
  peptidesL<-data %>% dplyr::select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
    gather(key="rep",value="Intensity",-Name)
  peptidesL<-pheno %>% left_join(groups) %>% left_join(peptidesL,by=c("newName"="rep"))
  peptidesL$GroupTime<-paste(peptidesL$Group,peptidesL$timept,sep=".")
  peptidesL$GroupTime<-factor(peptidesL$GroupTime,
                              levels=c("sCAD.FU","sCAD.T0","Type 1.FU","Type 1.T0","Type 2.FU","Type 2.T0",
                                       "Indeterminate.FU","Indeterminate.T0"))
  
  peptidesL<-peptidesL %>% arrange(GroupTime,ptid,replicate)
  temp1<-peptidesL %>% dplyr::select(GroupTime,ptid,replicate) %>% unique()
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
pepAnno<-peptides %>% dplyr::select(Name,Quality.Score,ParentProtein.FullName,Use.For.Quant,
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
# png(file="Plots/peptideNoNorm.png",height=5,width=10,units="in",res=300)
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
# png(file="Plots/peptideColNorm.png",height=5,width=10,units="in",res=300)
p1<-ggplot(data=peptidesL1,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p1)
# dev.off()

# png(file="Plots/peptideColMedNorm.png",height=5,width=10,units="in",res=300)
p1b<-ggplot(data=peptidesL1b,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p1b)
# dev.off()

# png(file="Plots/peptideNoneVColNorm.png",height=10,width=10,units="in",res=300)
grid.arrange(p0,p1,nrow=2)
# dev.off()

# Quantile normalization:
m2<-preprocessCore::normalize.quantiles(m0)
peptides2<-peptides
peptides2[,grepl("rep",names(peptides2))]<-m2
peptidesL2<-makeWideFun(peptides2)
# png(file="Plots/peptideQuantNorm.png",height=5,width=10,units="in",res=300)
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
# png(file="Plots/peptideFastCyclicLoess.png",height=5,width=10,units="in",res=300)
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
# png(file="Plots/peptideFastCyclicLoessbGalOnly.png",height=5,width=10,units="in",res=300)
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
# png(file="Plots/peptideFastCyclicLoessbGalHeavy.png",height=5,width=10,units="in",res=300)
p5<-ggplot(data=peptidesL5,aes(x=uID,group=newName,color=GroupTime,y=Intensity))+
  geom_boxplot()+theme_bw()+xlab("")+ylab("Intensity (log scale)")+
  theme(axis.text.x=element_blank())
show(p5)
# dev.off()

# MAD:
m6<-limma::normalizeMedianAbsValues(m0)
peptides6<-peptides
peptides6[,grepl("rep",names(peptides6))]<-m6
peptidesL6<-makeWideFun(peptides6)
# png(file="Plots/peptideMAD.png",height=5,width=10,units="in",res=300)
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
  bGalL<-bGal %>% dplyr::select(-Quality.Score,-(Parent.Protein:Percent.Files.With.Good.Quant)) %>%
    gather(key="rep",value="Intensity",-Name)
  bGalLNames<-bGalL %>% dplyr::select(Name) %>% unique()
  bGalLNames$id<-as.factor(1:nrow(bGalLNames))
  bGalL<-bGalLNames %>% left_join(bGalL)
  
  # CVs:
  cv<-bGalL %>% group_by(Name) %>% summarize(cv=sd(2**Intensity)/mean(2**Intensity)*100)
  names(cv)[names(cv)=="cv"]<-paste0(method,"CV")
  
  # Plots:
  fName<-paste0("Plots/bGalPeptides_",method,"_Full.png")
  # png(filename=fName,height=5,width=8,units="in",res=300)
  bGalp<-ggplot(bGalL,aes(x=rep,y=Intensity,group=Name,color=id))+
    geom_line()+ylim(23,32)+theme_bw()+xlab("Replicate")+ylab("Intensity (log scale)")
  print(bGalp)
  # dev.off()
  fName2<-paste0("Plots/bGalPeptides_",method,"_Samp.png")
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
peptides6[,grepl("rep_",names(peptides6))]<-
  2**peptides6[,grepl("rep_",names(peptides6))]

########### Rownames ###########
rownames(peptides00)<-peptides00$Name
rownames(peptides1)<-peptides1$Name
rownames(peptides2)<-peptides2$Name
rownames(peptides3)<-peptides3$Name
rownames(peptides4)<-peptides4$Name
rownames(peptides5)<-peptides5$Name
rownames(peptides6)<-peptides6$Name

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
myEColiIntensPep6<-peptides6[peptides6$Name %in% myEColi$pepSeq,
                             grepl("rep_",names(peptides6))]

cvFun1<-function(data){
  return(sd(apply(data,2,sum))/mean(apply(data,2,sum)))
}
bGalProtCVs<-
  data.frame(technique=c("None","Column TI","Quantile","CyclicLoess","CyclicLoessBGal",
                       "CyclicLoessBGalLt","MAD"),
           cv=c(cvFun1(myEColiIntensPep00),cvFun1(myEColiIntensPep1),
                cvFun1(myEColiIntensPep2),cvFun1(myEColiIntensPep3),
                cvFun1(myEColiIntensPep4),cvFun1(myEColiIntensPep5),
                cvFun1(myEColiIntensPep6)))
# write.csv(bGalProtCVs,file="bGalProtCVs.csv",row.names=FALSE)

########### Beta gal peptide correlations ###########
myEColiIntensPep1<-rbind(myEColiIntensPep1,apply(myEColiIntensPep1,2,sum))
rownames(myEColiIntensPep1)[nrow(myEColiIntensPep1)]<-"Total"
cor1<-cor(t(myEColiIntensPep1))

# png("Plots/cor1.png",height=6,width=6,units="in",res=300)
corrplot(cor1,type="upper",tl.cex=.4,is.corr=FALSE,cl.lim=c(-.4,1))
# dev.off()

myEColiIntensPep5<-rbind(myEColiIntensPep5,apply(myEColiIntensPep5,2,sum))
rownames(myEColiIntensPep5)[nrow(myEColiIntensPep5)]<-"Total"
cor5<-cor(t(myEColiIntensPep5))

# png("Plots/cor5.png",height=6,width=6,units="in",res=300)
corrplot(cor5,type="upper",tl.cex=.4,is.corr=FALSE,cl.lim=c(-.4,1))
# dev.off()

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
pep1<-combFun(Names=peptides1$Name,data=peptides1)

# Major cleanup:
rm(bGalCVs,bGalProtCVs,cor1,cor5,m0,m00,m1,m1b,m2,m3,m4,m5,m6,
   myEColi,myEColiIntensPep00,myEColiIntensPep1,myEColiIntensPep2,
   myEColiIntensPep3,myEColiIntensPep4,myEColiIntensPep5,myEColiIntensPep6,
   names2,p0,p1,p1b,p2,p3,p4,p5,p6,peptides,peptides00,peptides1b,
   peptides2,peptides3,peptides4,peptides5,peptides6,peptidesL,peptidesL1b,
   peptidesL2,peptidesL3,peptidesL4,peptidesL5,peptidesL6,proteins,
   reg1DF,bGalWeights,bGalWeights2,cFac,cFac2,cMed,cSums,i,mins,names1,
   pepSeqs,pepSeqsStr,reg1,reg2,bGalFun,combFun,cvFun1,makeWideFun,
   pinEcoli)

########### Peptide difference at baseline ###########
pepDF<-pep1 %>% gather(key="rep",value="Intensity",-Name)
pepDF$ptid<-str_split(pepDF$rep,"_",simplify=TRUE)[,2]
pepDF$timept<-str_split(pepDF$rep,"_",simplify=TRUE)[,3]
pepDF<-pepDF %>% left_join(groups)

unqPep<-unique(pepDF$Name)
pepDFT0Res<-data.frame(unqPep=unqPep,T0_sCAD=NA,T0_Type1=NA,T0_Type2=NA,
                     T0_Anova=NA,T0_Type1_sCAD=NA,T0_Type2_sCAD=NA,T0_Type1_Type2=NA, 
                     T0_Type1_sCAD_p=NA,T0_Type2_sCAD_p=NA,T0_Type1_Type2_p=NA)
for(i in 1:length(unqPep)){
  # Linear Model
  lm1<-lm(Intensity~Group,data=pepDF %>% 
            filter(Name==unqPep[i] & Group!="Indeterminate" & timept=="T0"))
  
  # Overall T0 ANOVA:
  lm1FStat<-summary(lm1)$fstatistic
  pepDFT0Res$T0_Anova[i]<-pf(lm1FStat[1],lm1FStat[2],lm1FStat[3],lower.tail=FALSE)
  
  # T0 Means:
  lm1Emmeans<-as.data.frame(emmeans(lm1,~Group))
  pepDFT0Res$T0_sCAD[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="sCAD"]
  pepDFT0Res$T0_Type1[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 1"]
  pepDFT0Res$T0_Type2[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 2"]

  # Pairwise T0: 
  lm1Pairs<-as.data.frame(pairs(emmeans(lm1,~Group),adjust="none"))
  pepDFT0Res$T0_Type1_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 1"])
  pepDFT0Res$T0_Type2_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 2"])
  pepDFT0Res$T0_Type1_Type2[i]<-
    (lm1Pairs$estimate[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  # Pairwise T0 p-value
  pepDFT0Res$T0_Type1_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 1"])
  pepDFT0Res$T0_Type2_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 2"])
  pepDFT0Res$T0_Type1_Type2_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  print(i)
}
pepDFT0Res<-pepDFT0Res %>% left_join(pepAnno2,by=c("unqPep"="pepSeq"))
pepDFT0ResGood<-pepDFT0Res %>% 
  filter(goodQuant>.8 & T0_Type1_sCAD_p<.1 & T0_Type1_Type2_p<.1)

########### Peptide Temporal Difference analysis ###########
pepDFw<-pepDF %>% dplyr::select(-rep) %>% 
  tidyr::spread(key="timept",value="Intensity")
pepDFw$d<-pepDFw$T0-pepDFw$FU

pepDFDRes<-data.frame(unqPep=unqPep,D_sCAD=NA,D_Type1=NA,D_Type2=NA,
                       D_Anova=NA,D_Type1_sCAD=NA,D_Type2_sCAD=NA,D_Type1_Type2=NA, 
                       D_Type1_sCAD_p=NA,D_Type2_sCAD_p=NA,D_Type1_Type2_p=NA)
for(i in 1:length(unqPep)){
  # Linear Model
  lm1<-lm(d~Group,data=pepDFw %>% 
            filter(Name==unqPep[i] & Group!="Indeterminate"))
  
  # Overall T0 ANOVA:
  lm1FStat<-summary(lm1)$fstatistic
  pepDFDRes$D_Anova[i]<-pf(lm1FStat[1],lm1FStat[2],lm1FStat[3],lower.tail=FALSE)
  
  # Time-point Means:
  lm1Emmeans<-as.data.frame(emmeans(lm1,~Group))
  pepDFDRes$D_sCAD[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="sCAD"]
  pepDFDRes$D_Type1[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 1"]
  pepDFDRes$D_Type2[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 2"]
  
  # Pairwise D: 
  lm1Pairs<-as.data.frame(pairs(emmeans(lm1,~Group),adjust="none"))
  pepDFDRes$D_Type1_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 1"])
  pepDFDRes$D_Type2_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 2"])
  pepDFDRes$D_Type1_Type2[i]<-
    (lm1Pairs$estimate[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  # Pairwise T0 p-value
  pepDFDRes$D_Type1_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 1"])
  pepDFDRes$D_Type2_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 2"])
  pepDFDRes$D_Type1_Type2_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  print(i)
}
pepDFDRes<-pepDFDRes %>% left_join(pepAnno2,by=c("unqPep"="pepSeq"))
pepDFDResGood<-pepDFDRes %>% 
  filter(goodQuant>.8 & D_Type1_sCAD_p<.1 & D_Type1_Type2_p<.1)

rm(lm1,lm1Emmeans,lm1Pairs,i,lm1FStat,unqPep)

########### Protein Aggregation ###########
protList<-paste(pepAnno2$proteins,collapse=";")
protList<-unique(unlist(str_split(protList,";")))
Prot<-data.frame(prot=protList,lvl=NA,nPep=NA,peps=NA)
for(prot in Prot$prot){
  peps<-pepAnno2$pepSeq[grepl(prot,pepAnno2$proteins,fixed=TRUE) & pepAnno2$protN==1 &
             pepAnno2$goodQuant>.3]
  pepsAll<-pepAnno2$pepSeq[grepl(prot,pepAnno2$proteins,fixed=TRUE)]
  if(length(peps)>0){
    Prot$peps[Prot$prot==prot]<-paste(peps,collapse=";")
    Prot$nPep[Prot$prot==prot]<-length(peps)
  }
  if(length(pepsAll)>0){
    Prot$pepsAll[Prot$prot==prot]<-paste(pepsAll,collapse=";")
    Prot$nPepAll[Prot$prot==prot]<-length(pepsAll)
  }
}
Prot0<-Prot
ProtOut<-Prot[is.na(Prot$nPep),]
Prot<-Prot[!is.na(Prot$nPep),]

# Calculate mean to make the protein abundances
pepsInProtList<-list()
for(i in 1:nrow(Prot)){
  pepsInProt<-unlist(str_split(Prot$peps[i],";"))
  pepsInProtDF<-pep1[pep1$Name %in% pepsInProt,]
  pepsInProtDF$Name<-NULL
  pepsInProtDF<-t(t(apply(pepsInProtDF,2,mean)))
  pepsInProtDF<-as.data.frame(pepsInProtDF)
  names(pepsInProtDF)[1]<-"value"
  pepsInProtDF$prot<-Prot$prot[i]
  pepsInProtDF$rep<-rownames(pepsInProtDF)
  pepsInProtList[[i]]<-pepsInProtDF
}
prots<-do.call("rbind",pepsInProtList)

# Add sample annotation:
prots$rep<-gsub("rep_","",prots$rep)
prots<-prots %>% left_join(pheno %>% dplyr::select(uSamp,Group,ptid,timept) 
                           %>% unique(),
                           by=c("rep"="uSamp"))
unqProts<-unique(prots$prot)

########### Join prot data to peptide data ###########
# Baseline abundances:
pepDFT0Res$OtherPepGood<-pepDFT0Res$OtherPepTotal<-pepDFT0Res$otherCor<-pepDFT0Res$otherGoodCor<-NA
for(i in 1:nrow(pepDFT0Res)){
  tempProts<-pepDFT0Res$proteins[i]
  tempProts<-unlist(str_split(tempProts,";"))
  tempProts<-str_split(str_split(tempProts,"\\|",simplify=TRUE)[,2],
                       "-",simplify=TRUE)[,1]
  allPepsDF<-pepDFT0Res[grepl(paste(tempProts,collapse="|"),
                              pepDFT0Res$proteins),]
  allPeps<-unique(allPepsDF$unqPep)
  allPepsGood<-unique(allPepsDF[allPepsDF$goodQuant>.3,]$unqPep)
  pepDFT0Res$OtherPepTotal[i]<-length(allPeps)
  pepDFT0Res$OtherPepGood[i]<-length(allPepsGood)
  
  # Correlation analysis:
  mat1<-as.matrix(pep1[pep1$Name %in% allPeps,names(pep1)!="Name"])
  mat2<-as.matrix(pep1[pep1$Name %in% allPepsGood,names(pep1)!="Name"])
  corMat1<-cor(t(mat1))
  corMat2<-cor(t(mat2))
  pepDFT0Res$otherCor[i]<-mean(corMat1[rownames(corMat1)!=pepDFT0Res$unqPep[i],
          colnames(corMat1)==pepDFT0Res$unqPep[i]])
  pepDFT0Res$otherGoodCor[i]<-mean(corMat2[rownames(corMat2)!=pepDFT0Res$unqPep[i],
                                       colnames(corMat2)==pepDFT0Res$unqPep[i]])
  print(i)
}

# Change across time:
pepDFDRes$OtherPepGood<-pepDFDRes$OtherPepTotal<-pepDFDRes$otherCor<-pepDFDRes$otherGoodCor<-NA
for(i in 1:nrow(pepDFDRes)){
  tempProts<-pepDFDRes$proteins[i]
  tempProts<-unlist(str_split(tempProts,";"))
  tempProts<-str_split(str_split(tempProts,"\\|",simplify=TRUE)[,2],
                       "-",simplify=TRUE)[,1]
  allPepsDF<-pepDFDRes[grepl(paste(tempProts,collapse="|"),
                              pepDFDRes$proteins),]
  allPeps<-unique(allPepsDF$unqPep)
  allPepsGood<-unique(allPepsDF[allPepsDF$goodQuant>.3,]$unqPep)
  pepDFDRes$OtherPepTotal[i]<-length(allPeps)
  pepDFDRes$OtherPepGood[i]<-length(allPepsGood)
  
  # Correlation analysis:
  mat1<-as.matrix(pep1[pep1$Name %in% allPeps,names(pep1)!="Name"])
  mat2<-as.matrix(pep1[pep1$Name %in% allPepsGood,names(pep1)!="Name"])
  corMat1<-cor(t(mat1))
  corMat2<-cor(t(mat2))
  pepDFDRes$otherCor[i]<-median(corMat1[rownames(corMat1)!=pepDFDRes$unqPep[i],
                                       colnames(corMat1)==pepDFDRes$unqPep[i]])
  pepDFDRes$otherGoodCor[i]<-median(corMat2[rownames(corMat2)!=pepDFDRes$unqPep[i],
                                           colnames(corMat2)==pepDFDRes$unqPep[i]])
  print(i)
}

# Join:
pepDFRes<-pepDFT0Res %>% full_join(
  pepDFDRes %>% dplyr::select("unqPep","D_sCAD","D_Type1","D_Type2","D_Anova",
                              "D_Type1_sCAD","D_Type2_sCAD","D_Type1_Type2",
                              "D_Type1_sCAD_p","D_Type2_sCAD_p","D_Type1_Type2_p"),
  by=c("unqPep"))

########### Add peptide criteria ###########
# Criteria 1: Fold change across time in Type 1 (ln-scale) > 1; 
# Diff between sCAD and Type 1 fold change across time (ln-scale) > 0.75
# pepDFRes$crit1 <- abs(pepDFRes$D_Type1 > 1) & abs(pepDFRes$D_Type1_sCAD > .75)

# New Criteria 1: 
pepDFRes$crit1<-pepDFRes$D_Type1_sCAD_p<.05

# Criteria 2: If >1 other peptide from protein or family, correlation > .5?
pepDFRes$crit2 <- FALSE
pepDFRes$crit2[(pepDFRes$OtherPepTotal > 1 & pepDFRes$otherCor > .5) | 
                 (pepDFRes$OtherPepTotal ==1)] <- TRUE

# Criteria 3: No missed cleavages
pepDFRes$crit3 <- pepDFRes$miss == 1

# Criteria 4: No oxidized Met
pepDFRes$crit4 <- !grepl("M\\[Oxid\\]",pepDFRes$unqPep)

# Criteria 5: If #1-#4
pepDFRes$crit5 <- pepDFRes$crit1 & pepDFRes$crit2 & pepDFRes$crit3 & pepDFRes$crit4

# Criteria 6: If different between type 1 MI and sCAD at T0:
# pepDFRes$crit6 <- abs(pepDFRes$T0_Type1_sCAD)>.5

# New Criteria 6:
pepDFRes$crit6<-pepDFRes$T0_Type1_sCAD_p<.1

# Criteria 7: Criteria 5 and 6:
pepDFRes$crit7 <- pepDFRes$crit5 & pepDFRes$crit6

# Export peptide results:
# write.csv(pepDFRes,"pepDFResV2.csv")

########### Protein level analysis ###########
# Baseline
protDFT0Res<-data.frame(prot=unqProts,T0_sCAD=NA,T0_Type1=NA,T0_Type2=NA,
                       T0_Anova=NA,T0_Type1_sCAD=NA,T0_Type2_sCAD=NA,T0_Type1_Type2=NA, 
                       T0_Type1_sCAD_p=NA,T0_Type2_sCAD_p=NA,T0_Type1_Type2_p=NA)
for(i in 1:nrow(protDFT0Res)){
  # Linear Model
  lm1<-lm(value~Group,data=prots %>% 
            filter(prot==protDFT0Res$prot[i] & Group!="Indeterminate" & timept=="T0"))
  
  # Overall T0 ANOVA:
  lm1FStat<-summary(lm1)$fstatistic
  protDFT0Res$T0_Anova[i]<-pf(lm1FStat[1],lm1FStat[2],lm1FStat[3],lower.tail=FALSE)
  
  # T0 Means:
  lm1Emmeans<-as.data.frame(emmeans(lm1,~Group))
  protDFT0Res$T0_sCAD[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="sCAD"]
  protDFT0Res$T0_Type1[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 1"]
  protDFT0Res$T0_Type2[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 2"]
  
  # Pairwise T0: 
  lm1Pairs<-as.data.frame(pairs(emmeans(lm1,~Group),adjust="none"))
  protDFT0Res$T0_Type1_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 1"])
  protDFT0Res$T0_Type2_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 2"])
  protDFT0Res$T0_Type1_Type2[i]<-
    (lm1Pairs$estimate[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  # Pairwise T0 p-value
  protDFT0Res$T0_Type1_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 1"])
  protDFT0Res$T0_Type2_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 2"])
  protDFT0Res$T0_Type1_Type2_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  print(i)
}

############ Protein Change across time ############
protsW<-prots %>% dplyr::select(-rep) %>% 
  tidyr::spread(key="timept",value="value")
protsW$d<-protsW$T0-protsW$FU

protDFDRes<-data.frame(prot=unqProts,D_sCAD=NA,D_Type1=NA,D_Type2=NA,
                        D_Anova=NA,D_Type1_sCAD=NA,D_Type2_sCAD=NA,D_Type1_Type2=NA, 
                        D_Type1_sCAD_p=NA,D_Type2_sCAD_p=NA,D_Type1_Type2_p=NA)
for(i in 1:nrow(protDFDRes)){
  # Linear Model
  lm1<-lm(d~Group,data=protsW %>% 
            filter(prot==protDFDRes$prot[i] & Group!="Indeterminate"))
  
  # Overall T0 ANOVA:
  lm1FStat<-summary(lm1)$fstatistic
  protDFDRes$D_Anova[i]<-pf(lm1FStat[1],lm1FStat[2],lm1FStat[3],lower.tail=FALSE)
  
  # D Means:
  lm1Emmeans<-as.data.frame(emmeans(lm1,~Group))
  protDFDRes$D_sCAD[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="sCAD"]
  protDFDRes$D_Type1[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 1"]
  protDFDRes$D_Type2[i]<-lm1Emmeans$emmean[lm1Emmeans$Group=="Type 2"]
  
  # Pairwise D: 
  lm1Pairs<-as.data.frame(pairs(emmeans(lm1,~Group),adjust="none"))
  protDFDRes$D_Type1_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 1"])
  protDFDRes$D_Type2_sCAD[i]<-
    (-lm1Pairs$estimate[lm1Pairs$contrast=="sCAD - Type 2"])
  protDFDRes$D_Type1_Type2[i]<-
    (lm1Pairs$estimate[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  # Pairwise D p-value
  protDFDRes$D_Type1_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 1"])
  protDFDRes$D_Type2_sCAD_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="sCAD - Type 2"])
  protDFDRes$D_Type1_Type2_p[i]<-
    (lm1Pairs$p.value[lm1Pairs$contrast=="Type 1 - Type 2"])
  
  print(i)
}

# Save
rm(corMat1,corMat2,lm1,lm1Emmeans,lm1FStat,lm1Pairs,mat1,mat2,allPeps,
   allPepsGood,i,peps,pepsAll,pepsInProt,prot,protList,tempProts,unqProts)
save.image(file="working_20181007.RData")

########### Peptide plots ###########
load(file="working_20181007.RData")

sigPeps<-pepDFRes$unqPep[pepDFRes$crit7 & pepDFRes$length<25 & pepDFRes$goodQuant > .05]
sigPeps2<-pepDFRes$unqPep[grepl("P04275",pepDFRes$proteins)]
sigPeps<-c(sigPeps,sigPeps2)
pepsSig<-pepDFRes[pepDFRes$unqPep %in% sigPeps,]
customTitles<-matrix(c(
  "NSEEFAAAMSR","Apolipoprotein B-100",
  "YYTYLIMNK","Complement C3",
  "C[Carboxymethyl]EWETPEGC[Carboxymethyl]EQVLTGK","C4b-binding protein alpha chain",
  "LQQVLHAGSGPC[Carboxymethyl]LPHLLSR","Alpha-2-antiplasmin",
  "LC[Carboxymethyl]MGSGLNLC[Carboxymethyl]EPNNK","Serotransferrin",
  "KPVAFSDYIHPVC[Carboxymethyl]LPDR","Prothrombin",
  "YEITTIHNLFR","Heparin cofactor 2",
  "EVVADSVWVDVK","Complement C3",
  "TVMVNIENPEGIPVK","Complement C3",
  "TSSFALNLPTLPEVK","Apolipoprotein B-100",
  "NFVASHIANILNSEELDIQDLK","Apolipoprotein B-100",
  "YTIAALLSPYSYSTTAVVTNPK","Transthyretin",
  "EALQGVGDMGR","Serum amyloid A-4 protein",
  "EWFSETFQK","Apolipoprotein C-I",
  "NIQEYLSILTDPDGK","Apolipoprotein B-100",
  "SGSSTASWIQNVDTK","Apolipoprotein B-100",
  "AGPHC[Carboxymethyl]PTAQLIATLK","Platelet Factor 4",
  "AASGTTGTYQEWK","Apolipoprotein B-100",
  "GTHC[Carboxymethyl]NQVEVIATLK","Platelet basic protein",
  "LIVAMSSWLQK","Apolipoprotein B-100",
  "HIQNIDIQHLAGK","Apolipoprotein B-100",
  "EDVYVVGTVLR","C4b-binding protein alpha chain",
  "EELC[Carboxymethyl]TMFIR","Apolipoprotein B-100",
  "HITSLEVIK","Platelet factor 4",
  "SLMLHYEFLQR","Complement component C8 beta chain",
  "ALVEQGFTVPEIK","Apolipoprotein B-100",
  "AVLTIDEK","Alpha-1-antitrypsin",
  "LC[Carboxymethyl]GC[Carboxymethyl]LELLLQFDQK","RUN and FYVE domain-containing protein 4",
  "NIQSLEVIGK","Platelet basic protein",
  "ATAVVDGAFK","Peroxiredoxin-2",
  "VIVIPVGIGPHANLK","von Willebrand factor",
  "YTLFQIFSK","von Willebrand factor"
  ),
  ncol=2,byrow=TRUE)

customTitles<-as.data.frame(customTitles)
names(customTitles)<-c("unqPep","Title")
pepsSig<-pepsSig %>% left_join(customTitles)

# Plot function:
pFun1<-function(iter,ylower,yupper){
  temp1<-pepDF[pepDF$Name==pepsSig$unqPep[iter] & pepDF$Group != "Indeterminate",]
  
  ggplot(temp1,aes(timept,Intensity,color=Group,group=ptid))+
    geom_point()+geom_line()+ylim(ylower,yupper)+theme_bw()+facet_grid(~Group)+
    ggtitle(pepsSig$unqPep[iter],
            subtitle=paste0(pepsSig$Title[iter],". cor = ",round(pepsSig$otherCor[iter],3)))+
    xlab("Time-point")+ylab("Abundance")+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

}
pepsSig<-pepsSig %>% arrange(Title,unqPep)

# pFun1(1,-2,2)
# pFun1(2,-3,3)
pFun1(3,-1.5,1.75)
pFun1(4,-3,3)
pFun1(5,-3,4)
pFun1(6,-1.25,2.1)
pFun1(7,-2,2.25)
pFun1(8,-1,1.25)
pFun1(9,-1.5,3)
pFun1(10,-1.1,1.5)
pFun1(11,-3,3)
pFun1(12,-1.1,1.3)
pFun1(13,-1.2,1)
# pFun1(14,-3,3)
# pFun1(15,-3,3)
# pFun1(16,-3,3)
# pFun1(17,-3,3)
# pFun1(18,-3,3)
# pFun1(19,-3,3)
pFun1(20,-1.1,1)
# pFun1(21,-3,3)
pFun1(22,-7.5,5.25)
pFun1(23,-2,3)
pFun1(24,-6.5,5.25)
pFun1(25,-5,5)
pFun1(26,-.9,.7)
# pFun1(27,-3,3)
# pFun1(28,-3,3)
# pFun1(29,-3,3)
pFun1(30,-1.2,1)
pFun1(31,-5,3)
pFun1(32,-4.5,3.2)

########### Peptide pull protein analysis ###########
setwd("~/gdrive/AthroProteomics")
load(file="working_20180821.RData")

# Which proteins:
str_split(pepDFRes$proteins,";")