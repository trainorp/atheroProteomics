########### Prereqs ###########
options(stringsAsFactors=FALSE,scipen=600)
library(tidyverse)

setwd("~/gdrive/AthroProteomics/data/")

# Import data:
pepToProt<-read.csv("pepToProt_20180524.csv",header=FALSE,skip=3)
srma<-read.csv("HumanSRMAtlasPeptidesFinalAnnotated.csv",header=TRUE)
# Names & human only:
names(pepToProt)<-c("pepSeq","proteins")
pepToProt<-pepToProt[grepl("HUMAN",pepToProt$proteins),]
# Other peptide annotation data:
load("pepAnno.RData")

########### Counts ###########
pepAnno<-pepAnno %>% left_join(pepToProt,by=c("Name"="pepSeq"))
idk<-pepToProt %>% group_by(proteins) %>% mutate(n())

pepToProt<-pepToProt %>% left_join(srma,by=c("pepSeq"="sequence"))
