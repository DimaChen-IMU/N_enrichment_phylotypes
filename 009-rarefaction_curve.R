# ##Usearch10.exe
#usearch10 -otutab_norm otu_table.txt -output otu_table_norm.txt -sample_size 3000 
#usearch10 -alpha_div_rare otu_table_norm.txt -output alpha_rare.txt -method without_replacement

library(ggplot2)
library(readxl)
library(reshape2)
library(ggpubr)
library(gridExtra) 
library(ggalt)
library(export)

rare <- read_excel("Data/001_Ba_asv.xlsx",6)
rare <- rare[,-1]

design <- read_excel("Data/001_Ba_asv.xlsx",4)
rare$x <- rownames(rare)
rare_melt <- melt(rare, id.vars = c('x'))
rare_melt$x <- factor(rare_melt$x, levels=1:100)
colnames(rare_melt) <- c('x','SampleID','value')

new <- merge(design,rare_melt)

p1 <- ggplot(new,aes(x=x,y=value,group=SampleID))+
  geom_line(aes(linetype=Type, color=Gradient))+xlab('Rarefraction Percentage')+
  ylab('Richness (Observed OTUs)')+theme_bw()+
  scale_color_gradient(low="blue", high="red")+
  scale_x_discrete(breaks = c(1:10)*10, labels=c(1:10)*10)+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))
###
rare <- read_excel("Data/001_Fu_asv.xlsx",6)
rare <- rare[,-1]

design <- read_excel("Data/001_Fu_asv.xlsx",4)
rare$x <- rownames(rare)
rare_melt <- melt(rare, id.vars = c('x'))
rare_melt$x <- factor(rare_melt$x, levels=1:100)
colnames(rare_melt) <- c('x','SampleID','value')

new <- merge(design,rare_melt)

p2 <- ggplot(new,aes(x=x,y=value,group=SampleID))+
  geom_line(aes(linetype=Type, color=Gradient))+xlab('Rarefraction Percentage')+
  ylab('Richness (Observed OTUs)')+theme_bw()+
  scale_color_gradient(low="blue", high="red")+
  scale_x_discrete(breaks = c(1:10)*10, labels=c(1:10)*10)+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))

ggarrange(p1, p2,ncol = 2,widths = c(1, 1), nrow = 1)
graph2ppt(file="rarefraction.pptx", width=8, height=3)


unique(design$Gradient)
