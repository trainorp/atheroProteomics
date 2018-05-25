########### Prereqs ###########
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)

setwd("~/gdrive/AthroProteomics/data/")

# Import data:
pepToProt<-read.csv("pepToProt_20180524.csv",header=TRUE)
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
  temp2<-paste(temp$proteins,collapse=";")
  return(cbind(temp[1,names(temp)!="proteins"],
               data.frame(proteins=temp2,protN=length(temp$proteins))))
}
protCollapseDF<-lapply(uniquePep,protCollapse)
pepAnno<-do.call("rbind",protCollapseDF)
nrow(pepAnno) #3,544

########### Counts ###########
pepAnno2<-pepAnno %>% left_join(srma,by=c("pepSeq"="sequence"))
save(pepAnno2,file="pepAnno2.RData")
