library(readxl)
library(tibble)
library(dplyr)
library(plyr)
library(tidyverse)
library(vegan)
######Ba-All-dominant-rare######
#all
asv_table <- read_excel("Data/001_Ba_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
RA_asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))
all_Ba_RA <- RA_asv_table
all_T <- data.frame(SampleID=colnames(all_Ba_RA),t(all_Ba_RA))

#dominant and rare
asv_sample <- read_excel("Data/001_Ba_asv.xlsx",4)
RA_control <- subset(RA_asv_table, select = c(subset(asv_sample,Gradient=="0")$SampleID))

RA_asv_table$Average <- rowMeans(RA_control)
RA_asv_table$DR <- ifelse(RA_asv_table$Average >= 0.05/100, 'Dominant','Rare')
RA_asv_table <- subset(RA_asv_table, select = -c(Average))

DR <- data.frame(asv_id=row.names(RA_asv_table),DR=RA_asv_table$DR)
write.csv(DR,'DR.csv')

dominant_Ba_RA <- subset(RA_asv_table,DR=='Dominant')
rare_Ba_RA <- subset(RA_asv_table,DR=='Rare')

dominant_Ba_RA_T <- data.frame(SampleID=colnames(dominant_Ba_RA),t(dominant_Ba_RA))
rare_Ba_RA_T <- data.frame(SampleID=colnames(rare_Ba_RA),t(rare_Ba_RA))

######Combine with design######
asv_sample <- read_excel("Data/001_Ba_asv.xlsx",4)
NW <- subset(asv_sample, Type =='NW')
NW$Gradient_s <- (sqrt(NW$Gradient)-min(sqrt(NW$Gradient)))/(max(sqrt(NW$Gradient))-min(sqrt(NW$Gradient)))

Acid <- subset(asv_sample, Type =='Acid')
Acid$Gradient_s <- (Acid$Gradient-min(Acid$Gradient))/(max(Acid$Gradient)-min(Acid$Gradient))
asv_sample <- data.frame(rbind(NW,Acid))

#all,dominant,rare-design
all_Ba_RA <- merge(asv_sample,all_T)
dominant_Ba_RA <- merge(asv_sample,dominant_Ba_RA_T)
rare_Ba_RA <- merge(asv_sample,rare_Ba_RA_T)

#write.csv(all_Ba_RA,'all_Ba_RA.csv')
#write.csv(dominant_Ba_RA,'dominant_Ba_RA.csv')
#write.csv(rare_Ba_RA,'rare_Ba_RA.csv')

######Fu-All-dominant-rare######
#all
asv_table <- read_excel("Data/001_Fu_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
RA_asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))
all_Fu_RA <- RA_asv_table
all_T <- data.frame(SampleID=colnames(all_Fu_RA),t(all_Fu_RA))

#dominant and rare
asv_sample <- read_excel("Data/001_Fu_asv.xlsx",4)
RA_control <- subset(RA_asv_table, select = c(subset(asv_sample,Gradient=="0")$SampleID))

RA_asv_table$Average <- rowMeans(RA_control)
RA_asv_table$DR <- ifelse(RA_asv_table$Average >= 0.05/100, 'Dominant','Rare')
RA_asv_table <- subset(RA_asv_table, select = -c(Average))

DR <- data.frame(asv_id=row.names(RA_asv_table),DR=RA_asv_table$DR)
write.csv(DR,'DR.csv')

dominant_Fu_RA <- subset(RA_asv_table,DR=='Dominant')
rare_Fu_RA <- subset(RA_asv_table,DR=='Rare')

dominant_Fu_RA_T <- data.frame(SampleID=colnames(dominant_Fu_RA),t(dominant_Fu_RA))
rare_Fu_RA_T <- data.frame(SampleID=colnames(rare_Fu_RA),t(rare_Fu_RA))

######Combine with design######
asv_sample <- read_excel("Data/001_Fu_asv.xlsx",4)
NW <- subset(asv_sample, Type =='NW')
NW$Gradient_s <- (sqrt(NW$Gradient)-min(sqrt(NW$Gradient)))/(max(sqrt(NW$Gradient))-min(sqrt(NW$Gradient)))

Acid <- subset(asv_sample, Type =='Acid')
Acid$Gradient_s <- (Acid$Gradient-min(Acid$Gradient))/(max(Acid$Gradient)-min(Acid$Gradient))
asv_sample <- data.frame(rbind(NW,Acid))

#all,dominant,rare-design
all_Fu_RA <- merge(asv_sample,all_T)
dominant_Fu_RA <- merge(asv_sample,dominant_Fu_RA_T)
rare_Fu_RA <- merge(asv_sample,rare_Fu_RA_T)

#write.csv(all_Fu_RA,'all_Fu_RA.csv')
#write.csv(dominant_Fu_RA,'dominant_Fu_RA.csv')
#write.csv(rare_Fu_RA,'rare_Fu_RA.csv')