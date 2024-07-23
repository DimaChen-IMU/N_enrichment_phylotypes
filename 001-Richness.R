library(readxl)
library(tibble)
library(dplyr)
library(plyr)
library(tidyverse)
library(vegan)
######dominant-or-rare_Ba######
asv_table <- read_excel("Data/001_Ba_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
RA_asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))

asv_sample <- read_excel("Data/001_Ba_asv.xlsx",4)
RA_control <- subset(RA_asv_table, select = c(subset(asv_sample,Gradient=="0")$SampleID))

RA_asv_table$Average <- rowMeans(RA_control)
RA_asv_table$DR <- ifelse(RA_asv_table$Average >= 0.05/100, 'Dominant','Rare')
DR <- RA_asv_table[,61:62]

######richness_Ba######
asv_table <- read_excel("Data/001_Ba_asv.xlsx",2)
asv_id <- asv_table$ASV_ID;asv_table <- asv_table[,-1]
row.names(asv_table) <- asv_id
asv_flatten <-  as.data.frame(t(rrarefy(t(asv_table), min(colSums(asv_table)))))
asv_flatten <-merge(DR,asv_flatten,by=0);row.names(asv_flatten) <- asv_flatten$Row.names

Dominant <- subset(subset(asv_flatten,DR =='Dominant'),select =-c(Row.names,Average,DR))
Dominant_T <- data.frame(SampleID=colnames(Dominant),t(Dominant))

Rare <- subset(subset(asv_flatten,DR =='Rare'),select =-c(Row.names,Average,DR))
Rare_T <- data.frame(SampleID=colnames(Rare),t(Rare))

All <- subset(asv_flatten,select =-c(Row.names,Average,DR))
All_T <- data.frame(SampleID=colnames(All),t(All))

#Richness
Ba_all_rich <- ddply(All_T,~SampleID,function(x) {data.frame(Ba_all_rich=sum(x[-1]>0))})
Ba_D_rich <- ddply(Dominant_T,~SampleID,function(x) {data.frame(Ba_D_rich=sum(x[-1]>0))})
Ba_R_rich <- ddply(Rare_T,~SampleID,function(x) {data.frame(Ba_R_rich=sum(x[-1]>0))})

######dominant-or-rare_Fu######
asv_table <- read_excel("Data/001_Fu_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
RA_asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))

asv_sample <- read_excel("Data/001_Fu_asv.xlsx",4)
RA_control <- subset(RA_asv_table, select = c(subset(asv_sample,Gradient=="0")$SampleID))

RA_asv_table$Average <- rowMeans(RA_control)
RA_asv_table$DR <- ifelse(RA_asv_table$Average >= 0.05/100, 'Dominant','Rare')
DR <- RA_asv_table[,61:62]

######richness_Fu######
asv_table <- read_excel("Data/001_Fu_asv.xlsx",2)
asv_id <- asv_table$ASV_ID;asv_table <- asv_table[,-1]
row.names(asv_table) <- asv_id
asv_flatten <-  as.data.frame(t(rrarefy(t(asv_table), min(colSums(asv_table)))))
asv_flatten <-merge(DR,asv_flatten,by=0);row.names(asv_flatten) <- asv_flatten$Row.names

Dominant <- subset(subset(asv_flatten,DR =='Dominant'),select =-c(Row.names,Average,DR))
Dominant_T <- data.frame(SampleID=colnames(Dominant),t(Dominant))

Rare <- subset(subset(asv_flatten,DR =='Rare'),select =-c(Row.names,Average,DR))
Rare_T <- data.frame(SampleID=colnames(Rare),t(Rare))

All <- subset(asv_flatten,select =-c(Row.names,Average,DR))
All_T <- data.frame(SampleID=colnames(All),t(All))

#Richness
Fu_all_rich <- ddply(All_T,~SampleID,function(x) {data.frame(Fu_all_rich=sum(x[-1]>0))})
Fu_D_rich <- ddply(Dominant_T,~SampleID,function(x) {data.frame(Fu_D_rich=sum(x[-1]>0))})
Fu_R_rich <- ddply(Rare_T,~SampleID,function(x) {data.frame(Fu_R_rich=sum(x[-1]>0))})

#######Richness-Design-Combine######
asv_sample <- read_excel("Data/001_Ba_asv.xlsx",4)
NW <- subset(asv_sample, Type =='NW')
NW$Gradient_s <- (sqrt(NW$Gradient)-min(sqrt(NW$Gradient)))/(max(sqrt(NW$Gradient))-min(sqrt(NW$Gradient)))

Acid <- subset(asv_sample, Type =='Acid')
Acid$Gradient_s <- (Acid$Gradient-min(Acid$Gradient))/(max(Acid$Gradient)-min(Acid$Gradient))
asv_sample <- data.frame(rbind(NW,Acid))

df_list <- list(asv_sample, Ba_all_rich, Ba_D_rich,Ba_R_rich,Fu_all_rich,Fu_D_rich,Fu_R_rich)
Richness <- df_list %>% reduce(full_join, by='SampleID')
Richness_cal <- Richness
write.table(Richness,"Richness.txt")

#######Plot######
Richness <- Richness %>% gather(key='Group', value='Richness', Ba_all_rich:Fu_R_rich) 
Richness$Group <- factor(Richness$Group,levels = c('Ba_all_rich','Ba_D_rich','Ba_R_rich','Fu_all_rich','Fu_D_rich','Fu_R_rich'))
Richness$Type <- factor(Richness$Type,levels = c('NW','Acid'))

ggplot(Richness, aes(Gradient_s, Richness, colour = Type)) +
  geom_point(size=1) +
  geom_smooth(se = T, linewidth=0.8,method = lm,alpha=0.1,aes(fill=Type))+
  scale_color_manual(values=c("#339933","#FF9900"))+
  scale_fill_manual(values=c("#339933","#FF9900"))+
  facet_wrap(~Group,ncol=3, scales = "free_y",strip.position = "left",
             labeller = as_labeller(c(Ba_all_rich = "All bacterial phylotypes", Ba_D_rich = "Dominant bacterial phylotypes",
                                      Ba_R_rich ="Rare bacterial phylotypes",Fu_all_rich = "All fungal phylotypes",
                                      Fu_D_rich = "Dominant fungal phylotypes",Fu_R_rich = "Rare fungal phylotypes")))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  ylab(NULL)+ xlab(NULL)+theme(legend.position="none")

#export::graph2ppt(file="richness2.pptx", width=7, height=5.3)
######Statistics######
library(Rmisc)
library(multcompView)
library(lsmeans)
library(lme4)
library(MuMIn)
library(lmerTest)
#options(timeout=6000)
#packageurl <- "http://R-Forge.R-project.org/bin/windows/contrib/4.3/Matrix_1.7-0.zip"
#install.packages(packageurl, repos=NULL, type="source")
#install.packages("lme4", type = "source")

###NW###
Richness_NW <- subset(Richness_cal,Type=='NW')
#coefficents:
coefficient <- NULL
R2 <- NULL
for (i in colnames(Richness_NW[c(7:12)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Richness_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}
coefficent;R2

###Acid###
Richness_Acid <- subset(Richness_cal,Type=='Acid')
#coefficents:
coefficient <- NULL
R2 <- NULL
for (i in colnames(Richness_Acid[c(7:12)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Richness_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

coefficent;R2


