library(readxl)
library(tibble)
library(dplyr)
library(plyr)
library(tidyverse)
library(vegan)

env <- read_excel("Data/005_env.xlsx",1)
env_cal <- env

env <- env %>% gather(key='Group', value='env', pH:PSR) 

env$Group <- factor(env$Group,levels = c('pH','TIN','ANPP','SOC','TSN','PSR'))
env$Type <- factor(env$Type,levels = c('NW','Acid'))

ggplot(env, aes(Gradient_s, env, colour = Type)) +
  geom_point(size=1) +
  geom_smooth(se = T, linewidth=0.8,method = lm,alpha=0.1,aes(fill=Type))+
  scale_color_manual(values=c("#339933","#FF9900"))+
  scale_fill_manual(values=c("#339933","#FF9900"))+
  facet_wrap(~Group,ncol=3, scales = "free_y",strip.position = "left",
             labeller = as_labeller(c(pH = "Soil pH", TIN = "Soil inorganic nitrogen (mg kg-1)",
                                      ANPP ="ANPP (g m-2)",SOC = "Soil organic C (g kg-1)",
                                      TSN = "Total soil N (g kg-1)",PSR = "Plant species richness")))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  ylab(NULL)+ xlab(NULL)+theme(legend.position="none")

#export::graph2ppt(file="env.pptx", width=7, height=5.3)

######Statistics######
library(Rmisc)
library(multcompView)
library(lsmeans)
library(lme4)
library(MuMIn)
library(lmerTest)

###NW###
env_NW <- subset(env_cal,Type=='NW')
#coefficents:
coefficient <- NULL
R2 <- NULL
for (i in colnames(env_NW[c(7:12)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = env_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}
coefficient;R2

###Acid###
env_Acid <- subset(env_cal,Type=='Acid')
#coefficents:
coefficient <- NULL
R2 <- NULL
for (i in colnames(env_Acid[c(7:12)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = env_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

coefficient;R2
