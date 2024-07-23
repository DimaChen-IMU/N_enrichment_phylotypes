library(readxl)
library(tibble)
library(dplyr)
library(plyr)
library(tidyverse)
library(vegan)
library(grid)
library(gridExtra)
library(ggforce)

Ba_pf <- read_excel("Data/003_asv_pf.xlsx",1)
Fu_pf <- read_excel("Data/003_asv_pf.xlsx",2)
pf_id <- read_excel("Data/003_asv_pf.xlsx",3)

PF <- data.frame(rbind(Ba_pf,Fu_pf))

PF_sum <-PF %>% dplyr::group_by(PF_N,PF_A,DR,Class,taxon) %>% dplyr::summarise(RA = sum(RA),count = count(DR)) 
PF_sum <- data.frame(PF_N=PF_sum$PF_N,PF_A=PF_sum$PF_A,DR=PF_sum$DR,Class=PF_sum$Class, taxon=PF_sum$taxon, RA=PF_sum$RA, count=PF_sum$count$freq)

PF_sum <- merge(pf_id,PF_sum, by=c("PF_N","PF_A","DR","Class","taxon"), all = TRUE)
PF_sum[is.na(PF_sum)] <- 0

PF_sum$DR_PF_N = paste(PF_sum$DR, PF_sum$PF_N, sep="_")
PF_sum <- PF_sum %>% dplyr::group_by(DR_PF_N,taxon) %>% dplyr::mutate(count_sum = sum(count),RA_sum=sum(RA))
PF_sum$count100 <- PF_sum$count/PF_sum$count_sum
PF_sum$RA100 <- PF_sum$RA/PF_sum$RA_sum

PF_sum <- gather(PF_sum, key="Type", value="per", count100:RA100)

PF_sum$DR_PF_N <- factor(PF_sum$DR_PF_N, levels = c("All_Npos","All_Nnon","All_Nneg",
                                                    "Dominant_Npos","Dominant_Nnon","Dominant_Nneg",
                                                    "Rare_Npos","Rare_Nnon","Rare_Nneg"),
                         labels =c("All_Npos","All_Nnon","All_Nneg",
                                   "Dominant_Npos","Dominant_Nnon","Dominant_Nneg",
                                   "Rare_Npos","Rare_Nnon","Rare_Nneg"))

PF_sum$PF_A <- factor(PF_sum$PF_A, levels = c("Aneg","Anon","Apos"),labels =c("Aneg","Anon","Apos"))

ggplot(PF_sum, aes(fill=PF_A, y=per, x=DR_PF_N)) + 
  geom_bar(position="fill", stat="identity", width=0.5)+
  facet_wrap(Type~taxon,ncol=2, scales = "free_y")+
  scale_fill_manual(values =alpha(c('#7AA706','#FFC926','#CF37B7'),1))+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  scale_y_continuous(limits = c(0, 1),expand = c(0, 0))
#export::graph2ppt(file="RA.pptx", width=7, height=7)
###
Ba_pf <- read_excel("Data/003_asv_pf.xlsx",1)
Fu_pf <- read_excel("Data/003_asv_pf.xlsx",2)
pf_id <- read_excel("Data/003_asv_pf.xlsx",3)

PF <- data.frame(rbind(Ba_pf,Fu_pf))

PF_sum <-PF %>% dplyr::group_by(PF_N,PF_A,DR,Class,taxon) %>% dplyr::summarise(RA = sum(RA),count = count(DR)) 
PF_sum <- data.frame(PF_N=PF_sum$PF_N,PF_A=PF_sum$PF_A,DR=PF_sum$DR,Class=PF_sum$Class, taxon=PF_sum$taxon, RA=PF_sum$RA, count=PF_sum$count$freq)

PF_sum <- merge(pf_id,PF_sum, by=c("PF_N","PF_A","DR","Class","taxon"), all = TRUE)
PF_sum[is.na(PF_sum)] <- 0

#count_number per DR&N
PF_sum_count <- PF_sum %>% dplyr::group_by(DR,taxon,PF_N) %>% dplyr::mutate(RA_DR_count = sum(count))
PF_sum_count <- unique(data.frame(taxon =PF_sum_count$taxon, PF_N = PF_sum_count$PF_N,DR =PF_sum_count$DR,RA_DR_count=PF_sum_count$RA_DR_count))

#RA per DR&N
PF_sum_RA <- PF_sum %>% dplyr::group_by(DR,taxon) %>% dplyr::mutate(RA_DR_sum = sum(RA))
PF_sum_RA <- PF_sum_RA %>% dplyr::group_by(DR,taxon,PF_N) %>% dplyr::mutate(RA_DR_PF_N_sum = sum(RA))
PF_sum_RA$proportion <- PF_sum_RA$RA_DR_sum/60

PF_sum_RA <- unique(data.frame(taxon =PF_sum_RA$taxon, PF_N = PF_sum_RA$PF_N,DR =PF_sum_RA$DR,
                               RA_DR_sum=PF_sum_RA$RA_DR_sum,RA_DR_PF_N_sum=PF_sum_RA$RA_DR_PF_N_sum,
                               proportion =PF_sum_RA$proportion))
PF_sum_RA$per <- PF_sum_RA$proportion* PF_sum_RA$RA_DR_PF_N_sum/PF_sum_RA$RA_DR_sum



