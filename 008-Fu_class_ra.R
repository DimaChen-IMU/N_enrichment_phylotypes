library(readxl)
library(tibble)
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
###Dominant or Rare
asv_table <- read_excel("Data/001_Fu_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
RA_asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))
all_Fu_RA <- RA_asv_table
all_T <- data.frame(SampleID=colnames(all_Fu_RA),t(all_Fu_RA))

asv_sample <- read_excel("Data/001_Fu_asv.xlsx",4)
RA_control <- subset(RA_asv_table, select = c(subset(asv_sample,Gradient=="0")$SampleID))

RA_asv_table$Average <- rowMeans(RA_control)
RA_asv_table$DR <- ifelse(RA_asv_table$Average >= 0.05/100, 'Dominant','Rare')
RA_asv_table <- subset(RA_asv_table, select = -c(Average))

asv_table <- data.frame(RA=rowSums(RA_asv_table[,-61]),DR=RA_asv_table$DR)

###class
Fu_class <- data.frame(read_excel("Data/001_Fu_asv.xlsx",3)[,3])
Fu_class_id <- data.frame(read_excel("Data/001_Fu_asv.xlsx",3)[,1])
rownames(Fu_class) <- Fu_class_id$ASV_ID

Fu_class <- data.frame(merge(Fu_class,asv_table,by = 0)[,-1])
Fu_class$name <- gsub("^.{0,3}", "", Fu_class$Class)

Fu_class_sum <- aggregate(RA ~ name + DR, Fu_class, sum)
Fu_class_sum <- spread(Fu_class_sum,DR,RA)
Fu_class_sum[is.na(Fu_class_sum)] <- 0
Fu_class_sum <- gather(Fu_class_sum,DR,RA,Dominant:Rare)

Fu_class_sum$DR <- factor(Fu_class_sum$DR, levels = c("Rare","Dominant"),labels =c("Rare","Dominant"))

a1 <- ggplot(Fu_class_sum, aes(fill=DR, y=RA/60, x=reorder(name, -RA))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  scale_fill_manual(values =alpha(c('#A9D18E','#E99A77'),1))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  xlab('')+ylab('Relative abundance (%)')

a2 <- ggplot(Fu_class_sum, aes(fill=DR, y=RA, x=reorder(name, -RA))) + 
  geom_bar(position="fill", stat="identity")+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  scale_fill_manual(values =alpha(c('#A9D18E','#E99A77'),1))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  xlab('')+ylab('Relative abundance (%)')

gridExtra::grid.arrange(a1, a2, nrow=2)
#export::graph2ppt(file="class.pptx", width=7, height=8)