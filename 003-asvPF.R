library(Rmisc)
library(multcompView)
library(lsmeans)
library(lme4)
library(MuMIn)
library(lmerTest)

######Bacteria-statistics######

###Nitrogen
all_Ba_RA <- read_excel("Data/002_asv_RA.xlsx",1)
all_Ba_RA_NW <- subset(all_Ba_RA,Type=='NW')
count_save <- data.frame(ifelse(all_Ba_RA_NW[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
all_Ba_RA_NW <- data.frame(all_Ba_RA_NW[,c(1:6)],subset(all_Ba_RA_NW, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(all_Ba_RA_NW[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = all_Ba_RA_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_N <- cbind(coefficient,R2)
combine_N <- as.data.frame(combine_N[,c(1,5,6)])
combine_N$PF_N <- ifelse(combine_N$`Pr(>|t|)` > 0.05, 'Nnon', ifelse(combine_N$Estimate>0,'Npos', 'Nneg'))
combine_N <- combine_N %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE)

###Acid
all_Ba_RA_Acid <- subset(all_Ba_RA,Type=='Acid')
count_save <- data.frame(ifelse(all_Ba_RA_Acid[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
all_Ba_RA_Acid <- data.frame(all_Ba_RA_Acid[,c(1:6)],subset(all_Ba_RA_Acid, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(all_Ba_RA_Acid[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = all_Ba_RA_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}
combine_Acid <- cbind(coefficient,R2)
combine_Acid <- as.data.frame(combine_Acid[,c(1,5,6)])
combine_Acid$PF_A <- ifelse(combine_Acid$`Pr(>|t|)` > 0.05, 'Anon', ifelse(combine_Acid$Estimate>0,'Apos', 'Aneg'))
combine_Acid <- combine_Acid %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE)

###combine-Ba-Nitrogen-Acid

asv_table <- read_excel("Data/002_asv_RA.xlsx",7)

asv_table <- data.frame(asv_id=asv_table$asv_id,DR=asv_table$DR,RA=rowSums(asv_table[,-c(1,2)]))

df_list <- list(combine_N,combine_Acid,asv_table )
Class_Ba <- df_list %>% reduce(full_join, by='asv_id')
Class_Ba <- Class_Ba[,-c(1,2,5,6)]
Class_Ba$PF_N[is.na(Class_Ba$PF_N)] <- 'Nnon'
Class_Ba$PF_A[is.na(Class_Ba$PF_A)] <- 'Anon'

Class_Ba <- Class_Ba %>% mutate(Class = case_when(PF_N=='Npos'& PF_A=='Apos' ~ 'EP1', PF_N=='Npos'& PF_A=='Anon' ~ 'EP2', PF_N=='Npos'& PF_A=='Aneg' ~ 'EP3', 
                                        PF_N=='Nnon'& PF_A=='Apos' ~ 'EP4', PF_N=='Nnon'& PF_A=='Anon' ~ 'EP5', PF_N=='Nnon'& PF_A=='Aneg' ~ 'EP6', 
                                        PF_N=='Nneg'& PF_A=='Apos' ~ 'EP7', PF_N=='Nneg'& PF_A=='Anon' ~ 'EP8', PF_N=='Nneg'& PF_A=='Aneg' ~ 'EP9' ))

Class_Ba$taxon <- 'Ba'

write.csv(Class_Ba,'Class_Ba.csv')

######Fungi-statistics######

###Nitrogen
all_Fu_RA <- read_excel("Data/002_asv_RA.xlsx",4)
#Nitrogen
all_Fu_RA_NW <- subset(all_Fu_RA,Type=='NW')
count_save <- data.frame(ifelse(all_Fu_RA_NW[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
all_Fu_RA_NW <- data.frame(all_Fu_RA_NW[,c(1:6)],subset(all_Fu_RA_NW, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(all_Fu_RA_NW[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = all_Fu_RA_NW)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}

combine_N <- cbind(coefficient,R2)
combine_N <- as.data.frame(combine_N[,c(1,5,6)])
combine_N$PF_N <- ifelse(combine_N$`Pr(>|t|)` > 0.05, 'Nnon', ifelse(combine_N$Estimate>0,'Npos', 'Nneg'))
combine_N <- combine_N %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE)

#Acid
all_Fu_RA_Acid <- subset(all_Fu_RA,Type=='Acid')
count_save <- data.frame(ifelse(all_Fu_RA_Acid[,-c(1:6)] == 0,0,1))
save_id <- colnames(count_save)[which(colSums(count_save) >= 1)]
all_Fu_RA_Acid <- data.frame(all_Fu_RA_Acid[,c(1:6)],subset(all_Fu_RA_Acid, select = save_id))

#coefficents:
coefficient <- NULL
R2 <- NULL

for (i in colnames(all_Fu_RA_Acid[-c(1:6)])){ 
  y <- i
  f <- as.formula(paste(y, paste("Gradient_s + (1|Block)", collapse=" * "), sep=" ~ "))
  mod <- lmer(f, data = all_Fu_RA_Acid)
  sum <- summary(mod)
  coef <- sum$coefficients[2,]
  coef$asv_id <- c(y)
  coefficient <- rbind(coefficient,coef)
  Rsq <- r.squaredGLMM(mod)
  R2 <- rbind(R2,Rsq)
}
combine_Acid <- cbind(coefficient,R2)
combine_Acid <- as.data.frame(combine_Acid[,c(1,5,6)])
combine_Acid$PF_A <- ifelse(combine_Acid$`Pr(>|t|)` > 0.05, 'Anon', ifelse(combine_Acid$Estimate>0,'Apos', 'Aneg'))
combine_Acid <- combine_Acid %>% unnest(Estimate, .drop = FALSE) %>% unnest('Pr(>|t|)', .drop = FALSE) %>% unnest('asv_id', .drop = FALSE)

###combine-Fu-Nitrogen-Acid

asv_table <- read_excel("Data/002_asv_RA.xlsx",8)

asv_table <- data.frame(asv_id=asv_table$asv_id,DR=asv_table$DR,RA=rowSums(asv_table[,-c(1,2)]))

df_list <- list(combine_N,combine_Acid,asv_table )
Class_Fu <- df_list %>% reduce(full_join, by='asv_id')
Class_Fu <- Class_Fu[,-c(1,2,5,6)]
Class_Fu$PF_N[is.na(Class_Fu$PF_N)] <- 'Nnon'
Class_Fu$PF_A[is.na(Class_Fu$PF_A)] <- 'Anon'

Class_Fu <- Class_Fu %>% mutate(Class = case_when(PF_N=='Npos'& PF_A=='Apos' ~ 'EP1', PF_N=='Npos'& PF_A=='Anon' ~ 'EP2', PF_N=='Npos'& PF_A=='Aneg' ~ 'EP3', 
                                                  PF_N=='Nnon'& PF_A=='Apos' ~ 'EP4', PF_N=='Nnon'& PF_A=='Anon' ~ 'EP5', PF_N=='Nnon'& PF_A=='Aneg' ~ 'EP6', 
                                                  PF_N=='Nneg'& PF_A=='Apos' ~ 'EP7', PF_N=='Nneg'& PF_A=='Anon' ~ 'EP8', PF_N=='Nneg'& PF_A=='Aneg' ~ 'EP9' ))

Class_Fu$taxon <- 'Fu'

write.csv(Class_Fu,'Class_Fu.csv')
