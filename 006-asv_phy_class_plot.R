library(readxl)
library(tibble)
library(data.table)
library(dplyr)
library(plyr)
library(tidyverse)
library(vegan)
library(Rmisc)
library(multcompView)
library(lsmeans)
library(lme4)
library(MuMIn)
library(lmerTest)
######Richness-Bacterial-Phylum######
###Ba_phy_rich
asv_table <- read_excel("Data/004_asv_phy_class.xlsx",1)
phy_name <- unique(asv_table$SampleID)

phy_richness <- list()

for (i in seq(length(phy_name))){ 
  phy <- subset(asv_table,SampleID == phy_name[i])
  phy_t <- data.frame(t(phy))
  colnames(phy_t) <- phy_t[1,]
  phy_t <- data.frame(SampleID=row.names(phy_t),phy_t)
  phy_t <- phy_t[-1,]
  rich <- ddply(phy_t,~SampleID,function(x) {data.frame(richness = sum(x[-1]>0))})
  colnames(rich) <- c(paste0("ID_",i),phy_name[i])
  phy_richness[[i]] <- rich
}

phy_rich_mess <-  as.data.frame(do.call(cbind, phy_richness))
phy_richness <- data.frame(SampleID = phy_rich_mess$ID_1, phy_rich_mess[,seq(0, length(phy_rich_mess), 2)])

asv_sample <- read_excel("Data/004_asv_phy_class.xlsx",3)

phy_richness <- merge(asv_sample,phy_richness,by='SampleID')

###Ba_phy_rich_statistic
#Nitrogen
Ba_phy_NW <- subset(phy_richness,Type=='NW')
count_save <- data.frame(ifelse(Ba_phy_NW[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Ba_phy_NW <- data.frame(Ba_phy_NW[,c(1:6)],subset(Ba_phy_NW, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Ba_phy_NW[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Ba_phy_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_N <- cbind(coefficient,R2)
combine_N <- as.data.frame(combine_N[,c(1,5,6,7)])
combine_N <- combine_N %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_N$R <- ifelse(combine_N$Estimate >0, sqrt(combine_N$R2m),-sqrt(combine_N$R2m))
combine_N$PF_N <- ifelse(combine_N$`Pr(>|t|)` > 0.05, 'Nnon', ifelse(combine_N$Estimate>0,'Npos', 'Nneg'))
colnames(combine_N) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_n','PF_N')
combine_N <- combine_N[,c(3,5,6)]

#Acid
Ba_phy_Acid <- subset(phy_richness,Type=='Acid')
count_save <- data.frame(ifelse(Ba_phy_Acid[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Ba_phy_Acid <- data.frame(Ba_phy_Acid[,c(1:6)],subset(Ba_phy_Acid, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Ba_phy_Acid[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Ba_phy_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_Acid <- cbind(coefficient,R2)
combine_Acid <- as.data.frame(combine_Acid[,c(1,5,6,7)])
combine_Acid <- combine_Acid %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_Acid$R <- ifelse(combine_Acid$Estimate >0, sqrt(combine_Acid$R2m),-sqrt(combine_Acid$R2m))
combine_Acid$PF_A <- ifelse(combine_Acid$`Pr(>|t|)` > 0.05, 'Anon', ifelse(combine_Acid$Estimate>0,'Apos', 'Aneg'))
colnames(combine_Acid) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_a','PF_A')
combine_Acid <- combine_Acid[,c(3,5,6)]

Phylum_Ba_rich <- merge(combine_N,combine_Acid,by='asv_id')
Phylum_Ba_rich$Type <- 'Ba_rich'




######Richness-Fungal-Class######
###Fu_class_rich
asv_table <- read_excel("Data/004_asv_phy_class.xlsx",2)
class_name <- unique(asv_table$SampleID)

class_richness <- list()

for (i in seq(length(class_name))){ 
  class <- subset(asv_table,SampleID == class_name[i])
  class_t <- data.frame(t(class))
  colnames(class_t) <- class_t[1,]
  class_t <- data.frame(SampleID=row.names(class_t),class_t)
  class_t <- class_t[-1,]
  rich <- ddply(class_t,~SampleID,function(x) {data.frame(richness = sum(x[-1]>0))})
  colnames(rich) <- c(paste0("ID_",i),class_name[i])
  class_richness[[i]] <- rich
}

class_rich_mess <-  as.data.frame(do.call(cbind, class_richness))
class_richness <- data.frame(SampleID = class_rich_mess$ID_1, class_rich_mess[,seq(0, length(class_rich_mess), 2)])

asv_sample <- read_excel("Data/004_asv_phy_class.xlsx",3)

class_richness <- merge(asv_sample,class_richness,by='SampleID')

###Fu_class_rich_statistic
#Nitrogen
Fu_class_NW <- subset(class_richness,Type=='NW')
count_save <- data.frame(ifelse(Fu_class_NW[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Fu_class_NW <- data.frame(Fu_class_NW[,c(1:6)],subset(Fu_class_NW, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Fu_class_NW[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Fu_class_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_N <- cbind(coefficient,R2)
combine_N <- as.data.frame(combine_N[,c(1,5,6,7)])
combine_N <- combine_N %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_N$R <- ifelse(combine_N$Estimate >0, sqrt(combine_N$R2m),-sqrt(combine_N$R2m))
combine_N$PF_N <- ifelse(combine_N$`Pr(>|t|)` > 0.05, 'Nnon', ifelse(combine_N$Estimate>0,'Npos', 'Nneg'))
colnames(combine_N) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_n','PF_N')
combine_N <- combine_N[,c(3,5,6)]

#Acid
Fu_class_Acid <- subset(class_richness,Type=='Acid')
count_save <- data.frame(ifelse(Fu_class_Acid[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Fu_class_Acid <- data.frame(Fu_class_Acid[,c(1:6)],subset(Fu_class_Acid, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Fu_class_Acid[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Fu_class_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_Acid <- cbind(coefficient,R2)
combine_Acid <- as.data.frame(combine_Acid[,c(1,5,6,7)])
combine_Acid <- combine_Acid %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_Acid$R <- ifelse(combine_Acid$Estimate >0, sqrt(combine_Acid$R2m),-sqrt(combine_Acid$R2m))
combine_Acid$PF_A <- ifelse(combine_Acid$`Pr(>|t|)` > 0.05, 'Anon', ifelse(combine_Acid$Estimate>0,'Apos', 'Aneg'))
colnames(combine_Acid) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_a','PF_A')
combine_Acid <- combine_Acid[,c(3,5,6)]

Class_Fu_rich <- merge(combine_N,combine_Acid,by='asv_id')
Class_Fu_rich$Type <- 'Fu_rich'


######Abundance-Bacteria######
###Bacterial phylum
asv_table <- read_excel("Data/001_Ba_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))

Ba_phylum <- data.frame(read_excel("Data/001_Ba_asv.xlsx",3)[,2])
Ba_phylum_id <- data.frame(read_excel("Data/001_Ba_asv.xlsx",3)[,1])
rownames(Ba_phylum) <- Ba_phylum_id$ASV_ID

Ba_asv <- data.table(merge(Ba_phylum,asv_table,by = 0)[,-1])
Ba_phy <- Ba_asv[, lapply(.SD, sum, na.rm=TRUE), by=Phylum ]

Ba_phy <- Ba_phy %>% tibble::column_to_rownames("Phylum")
Ba_phy_T <- data.frame(SampleID=colnames(Ba_phy),t(Ba_phy))
###SampleID combine with Ba phuylum
asv_sample <- read_excel("Data/001_Ba_asv.xlsx",4)
NW <- subset(asv_sample, Type =='NW')
NW$Gradient_s <- (sqrt(NW$Gradient)-min(sqrt(NW$Gradient)))/(max(sqrt(NW$Gradient))-min(sqrt(NW$Gradient)))

Acid <- subset(asv_sample, Type =='Acid')
Acid$Gradient_s <- (Acid$Gradient-min(Acid$Gradient))/(max(Acid$Gradient)-min(Acid$Gradient))
asv_sample <- data.frame(rbind(NW,Acid))

df_list <- list(asv_sample, Ba_phy_T)
Ba_phy <- df_list %>% reduce(full_join, by='SampleID')

###Statistic-Ba

#Nitrogen
Ba_phy_NW <- subset(Ba_phy,Type=='NW')
count_save <- data.frame(ifelse(Ba_phy_NW[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Ba_phy_NW <- data.frame(Ba_phy_NW[,c(1:6)],subset(Ba_phy_NW, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Ba_phy_NW[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Ba_phy_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_N <- cbind(coefficient,R2)
combine_N <- as.data.frame(combine_N[,c(1,5,6,7)])
combine_N <- combine_N %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_N$R <- ifelse(combine_N$Estimate >0, sqrt(combine_N$R2m),-sqrt(combine_N$R2m))
combine_N$PF_N <- ifelse(combine_N$`Pr(>|t|)` > 0.05, 'Nnon', ifelse(combine_N$Estimate>0,'Npos', 'Nneg'))
colnames(combine_N) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_n','PF_N')
combine_N <- combine_N[,c(3,5,6)]

#Acid
Ba_phy_Acid <- subset(Ba_phy,Type=='Acid')
count_save <- data.frame(ifelse(Ba_phy_Acid[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Ba_phy_Acid <- data.frame(Ba_phy_Acid[,c(1:6)],subset(Ba_phy_Acid, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Ba_phy_Acid[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Ba_phy_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_Acid <- cbind(coefficient,R2)
combine_Acid <- as.data.frame(combine_Acid[,c(1,5,6,7)])
combine_Acid <- combine_Acid %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_Acid$R <- ifelse(combine_Acid$Estimate >0, sqrt(combine_Acid$R2m),-sqrt(combine_Acid$R2m))
combine_Acid$PF_A <- ifelse(combine_Acid$`Pr(>|t|)` > 0.05, 'Anon', ifelse(combine_Acid$Estimate>0,'Apos', 'Aneg'))
colnames(combine_Acid) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_a','PF_A')
combine_Acid <- combine_Acid[,c(3,5,6)]

Phylum_Ba_abun <- merge(combine_N,combine_Acid,by='asv_id')
Phylum_Ba_abun$Type <- 'Ba_abun'

######Abundance-Fungi######
###Fungal class
asv_table <- read_excel("Data/001_Fu_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))

Fu_class <- data.frame(read_excel("Data/001_Fu_asv.xlsx",3)[,3])
Fu_class_id <- data.frame(read_excel("Data/001_Fu_asv.xlsx",3)[,1])
rownames(Fu_class) <- Fu_class_id$ASV_ID

Fu_asv <- data.table(merge(Fu_class,asv_table,by = 0)[,-1])
Fu_class <- Fu_asv[, lapply(.SD, sum, na.rm=TRUE), by=Class ]

Fu_class <- Fu_class %>% tibble::column_to_rownames("Class")
Fu_class_T <- data.frame(SampleID=colnames(Fu_class),t(Fu_class))

###SampleID combine with Fu class
asv_sample <- read_excel("Data/001_Fu_asv.xlsx",4)
NW <- subset(asv_sample, Type =='NW')
NW$Gradient_s <- (sqrt(NW$Gradient)-min(sqrt(NW$Gradient)))/(max(sqrt(NW$Gradient))-min(sqrt(NW$Gradient)))

Acid <- subset(asv_sample, Type =='Acid')
Acid$Gradient_s <- (Acid$Gradient-min(Acid$Gradient))/(max(Acid$Gradient)-min(Acid$Gradient))
asv_sample <- data.frame(rbind(NW,Acid))

df_list <- list(asv_sample, Fu_class_T)
Fu_class <- df_list %>% reduce(full_join, by='SampleID')

###Statistic-Fu
library(Rmisc)
library(multcompView)
library(lsmeans)
library(lme4)
library(MuMIn)
library(lmerTest)
#Nitrogen
Fu_class_NW <- subset(Fu_class,Type=='NW')
count_save <- data.frame(ifelse(Fu_class_NW[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Fu_class_NW <- data.frame(Fu_class_NW[,c(1:6)],subset(Fu_class_NW, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Fu_class_NW[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Fu_class_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_N <- cbind(coefficient,R2)
combine_N <- as.data.frame(combine_N[,c(1,5,6,7)])
combine_N <- combine_N %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_N$R <- ifelse(combine_N$Estimate >0, sqrt(combine_N$R2m),-sqrt(combine_N$R2m))
combine_N$PF_N <- ifelse(combine_N$`Pr(>|t|)` > 0.05, 'Nnon', ifelse(combine_N$Estimate>0,'Npos', 'Nneg'))
colnames(combine_N) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_n','PF_N')
combine_N <- combine_N[,c(3,5,6)]

#Acid
Fu_class_Acid <- subset(Fu_class,Type=='Acid')
count_save <- data.frame(ifelse(Fu_class_Acid[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
Fu_class_Acid <- data.frame(Fu_class_Acid[,c(1:6)],subset(Fu_class_Acid, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(Fu_class_Acid[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = Fu_class_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_Acid <- cbind(coefficient,R2)
combine_Acid <- as.data.frame(combine_Acid[,c(1,5,6,7)])
combine_Acid <- combine_Acid %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE) %>% unnest('R2m', .drop = FALSE)
combine_Acid$R <- ifelse(combine_Acid$Estimate >0, sqrt(combine_Acid$R2m),-sqrt(combine_Acid$R2m))
combine_Acid$PF_A <- ifelse(combine_Acid$`Pr(>|t|)` > 0.05, 'Anon', ifelse(combine_Acid$Estimate>0,'Apos', 'Aneg'))
colnames(combine_Acid) <- c('Estimate','Pr(>|t|)','asv_id','R2m','R_a','PF_A')
combine_Acid <- combine_Acid[,c(3,5,6)]

Class_Fu_abun <- merge(combine_N,combine_Acid,by='asv_id')
Class_Fu_abun$Type <- 'Fu_abun'

###plot_data
phylum_class <- rbind(Phylum_Ba_rich,Class_Fu_rich,Phylum_Ba_abun,Class_Fu_abun)
phylum_class$label <- gsub("^.{0,3}", "", phylum_class$asv_id)
phylum_class$label <- toupper(substr(phylum_class$label, 1, 3))

phylum_class <- phylum_class %>% 
  mutate(PF =case_when(PF_N=='Npos'& PF_A=='Apos' ~ 'EP1', PF_N=='Npos'& PF_A=='Anon' ~ 'EP2', PF_N=='Npos'& PF_A=='Aneg' ~ 'EP3', 
                       PF_N=='Nnon'& PF_A=='Apos' ~ 'EP4', PF_N=='Nnon'& PF_A=='Anon' ~ 'EP5', PF_N=='Nnon'& PF_A=='Aneg' ~ 'EP6', 
                       PF_N=='Nneg'& PF_A=='Apos' ~ 'EP7', PF_N=='Nneg'& PF_A=='Anon' ~ 'EP8', PF_N=='Nneg'& PF_A=='Aneg' ~ 'EP9' ))


phylum_class$Type <- factor(phylum_class$Type, levels = c("Ba_rich","Fu_rich","Ba_abun","Fu_abun"),labels =c("Ba_rich","Fu_rich","Ba_abun","Fu_abun"))

rich_Ba <- subset(phylum_class,Type=='Ba_rich')
rich_Fu <- subset(phylum_class,Type=='Fu_rich')
abun_Ba <- subset(phylum_class,Type=='Ba_abun')
abun_Fu <- subset(phylum_class,Type=='Fu_abun')

p1 <- ggplot(rich_Ba, aes(x= R_n, y= R_a, colour=PF, label=label))+
  geom_point() +geom_text(hjust=0, vjust=0,size=2.5)+
  scale_color_manual(values =alpha(c('#FF9900','#1313FF','#653200','#339933','#612B00','#0000FF','#FF9900'),1))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  scale_y_continuous(limits = c(-0.8, 0.8))+scale_x_continuous(limits = c(-0.9, 0.6))+
  geom_hline(yintercept=0.32, linetype="dashed")+geom_hline(yintercept=-0.31, linetype="dashed")+
  geom_vline(xintercept=0.35, linetype="dashed")+geom_vline(xintercept=-0.36, linetype="dashed")+
  theme(legend.position = "none")+xlab('r of ASV richness in N addition')+ylab('r of ASV richness in A addition')


p2 <- ggplot(rich_Fu, aes(x= R_n, y= R_a, colour=PF, label=label))+
  geom_point() +geom_text(hjust=0, vjust=0,size=2.5)+
  scale_color_manual(values =alpha(c('#1313FF','#339933','#612B00','#0000FF','#FF9900'),1))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  scale_y_continuous(limits = c(-0.65, 0.5))+scale_x_continuous(limits = c(-0.8, 0.6))+
  geom_hline(yintercept=0.34, linetype="dashed")+geom_hline(yintercept=-0.34, linetype="dashed")+
  geom_vline(xintercept=0.34, linetype="dashed")+geom_vline(xintercept=-0.34, linetype="dashed")+
  theme(legend.position = "none")+xlab('r of ASV richness in N addition')+ylab('r of ASV richness in A addition')

p3 <- ggplot(abun_Ba, aes(x= R_n, y= R_a, colour=PF, label=label))+
  geom_point() +geom_text(hjust=0, vjust=0,size=2.5)+
  scale_color_manual(values =alpha(c('#FF9900','#1313FF','#653200','#339933','#612B00','#0000FF','#FF9900'),1))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  scale_y_continuous(limits = c(-0.8, 0.8))+scale_x_continuous(limits = c(-1, 0.8))+
  geom_hline(yintercept=0.35, linetype="dashed")+geom_hline(yintercept=-0.35, linetype="dashed")+
  geom_vline(xintercept=0.35, linetype="dashed")+geom_vline(xintercept=-0.35, linetype="dashed")+
  theme(legend.position = "none")+xlab('r of abundance of N addition')+ylab('r of abundance of A addition')

p4 <- ggplot(abun_Fu, aes(x= R_n, y= R_a, colour=PF, label=label))+
  geom_point() +geom_text(hjust=0, vjust=0,size=2.5)+
  scale_color_manual(values =alpha(c('#FF9900','#1313FF','#339933','#612B00','#0000FF','#FF9900'),1))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  scale_y_continuous(limits = c(-0.8, 0.8))+scale_x_continuous(limits = c(-0.8, 1))+
  geom_hline(yintercept=0.35, linetype="dashed")+geom_hline(yintercept=-0.35, linetype="dashed")+
  geom_vline(xintercept=0.35, linetype="dashed")+geom_vline(xintercept=-0.35, linetype="dashed")+
  theme(legend.position = "none")+xlab('r of abundance of N addition')+ylab('r of abundance of A addition')


gridExtra::grid.arrange(p1, p2,p3,p4, nrow=2)
#export::graph2ppt(file="phyclas.pptx", width=7, height=6.5)

ggplot(abun_Fu, aes(x= R_n, y= R_a, colour=PF, label=asv_id))+
  geom_point() +geom_text(hjust=0, vjust=0,size=2.5)+
  scale_color_manual(values =alpha(c('#FF9900','#1313FF','#339933','#612B00','#0000FF','#FF9900'),1))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  scale_y_continuous(limits = c(-0.8, 0.8))+scale_x_continuous(limits = c(-0.8, 1))+
  geom_hline(yintercept=0.35, linetype="dashed")+geom_hline(yintercept=-0.35, linetype="dashed")+
  geom_vline(xintercept=0.35, linetype="dashed")+geom_vline(xintercept=-0.35, linetype="dashed")+
  theme(legend.position = "none")+xlab('r of abundance of N addition')+ylab('r of abundance of A addition')
