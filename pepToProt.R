########### Prereqs ###########
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)

setwd("~/gdrive/AthroProteomics/data/")

# Import data:
pepToProt<-read.csv("pepToProt_20180524.csv",header=TRUE)
pepToProtwIso<-read.csv("pepToProt_20180603.csv",header=TRUE)
pepToProt<-rbind(pepToProt,pepToProtwIso)
pepToProt<-pepToProt[!duplicated(pepToProt),]

srma<-read.csv("HumanSRMAtlasPeptidesFinalAnnotated.csv",header=TRUE)
# Names & human only:
names(pepToProt)<-c("pepSeq","proteins")
# Find e coli beta gal / otherwise filter for Human only:
pepToProt<-pepToProt[grepl("HUMAN",pepToProt$proteins)|
                       grepl("P00722",pepToProt$protein),]
pepToProt<-pepToProt[!duplicated(pepToProt),]
# Other peptide annotation data:
load("pepAnno.RData")
nrow(pepAnno) #3,544

########### Collapse proteins ###########
pepAnno$pepSeq<-gsub("(\\[.*?\\])","",pepAnno$Name)
pepAnno<-pepAnno %>% left_join(pepToProt)

uniquePep<-unique(pepAnno$Name)
# Collapse function
protCollapse<-function(pep){
  temp<-pepAnno[pepAnno$Name==pep,]
  temp1<-temp$proteins
  if(length(temp1)>1){
    temp1b<-str_split(temp1,"\\|",simplify=TRUE)[,2]
    temp1c<-substr(temp1b,1,6)
    protN<-length(unique(temp1c))
  }
  else{
    protN<-1L
  }
  temp2<-paste(temp$proteins,collapse=";")
  return(cbind(temp[1,names(temp)!="proteins"],
               data.frame(proteins=temp2,protN=protN,protNIso=length(temp1))))
}
protCollapseDF<-lapply(uniquePep,protCollapse)
pepAnno<-do.call("rbind",protCollapseDF)
nrow(pepAnno) #3,544

########### Counts ###########
pepAnno2<-pepAnno %>% left_join(srma,by=c("pepSeq"="sequence"))
pepAnno2<-pepAnno2 %>% select(pepSeq=Name,pepSeq2=pepSeq,quality=Quality.Score,
                    goodQuant=Percent.Files.With.Good.Quant,
                    proteins,protN,protNIso,pinProtName=ParentProtein.FullName,
                    pinUseQuant=Use.For.Quant,
                    atlasProt=Prot_acc,SSR,length,type,PA_Acc,mw,e_charge,
                    n_map_core,n_map_all)

########### Digestion & missed cleavage ###########
# digest<-read.table("uniprot_sprot_wvar_digested_Mass400to6000.txt",header=TRUE)
# Wrote that part in python
# New import:
digest<-read.table("digest2.txt")
digest<-digest[,-c(3:5)]
names(digest)<-c("protein","seq","ind")
digest$miss<-str_split(digest$ind,"\\.",simplify=TRUE)[,2]
digest<-digest %>% dplyr::select(seq,miss) %>% unique()
pepAnno2<-pepAnno2 %>% left_join(digest,by=c("pepSeq2"="seq"))

########### Export ###########
write.csv(pepAnno2,file="pepAnno2.csv",row.names=FALSE)
save(pepAnno2,file="pepAnno2.RData")
