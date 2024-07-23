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
library(data.table)
######DR in each phylum or class######
###Dominant or Rare
asv_table <- read_excel("Data/001_Ba_asv.xlsx",2)
asv_table <- asv_table %>% tibble::column_to_rownames("ASV_ID")
RA_asv_table <- as.data.frame(apply(asv_table, 2, function(x){x/sum(x)}))
all_Ba_RA <- RA_asv_table
all_T <- data.frame(SampleID=colnames(all_Ba_RA),t(all_Ba_RA))

asv_sample <- read_excel("Data/001_Ba_asv.xlsx",4)
RA_control <- subset(RA_asv_table, select = c(subset(asv_sample,Gradient=="0")$SampleID))

RA_asv_table$Average <- rowMeans(RA_control)
RA_asv_table$DR <- ifelse(RA_asv_table$Average >= 0.05/100, 'Dominant','Rare')
RA_asv_table <- subset(RA_asv_table, select = -c(Average))

asv_table <- data.frame(RA=rowSums(RA_asv_table[,-61]),DR=RA_asv_table$DR)

###phylum
Ba_phylum <- data.frame(read_excel("Data/001_Ba_asv.xlsx",3)[,2])
Ba_phylum_id <- data.frame(read_excel("Data/001_Ba_asv.xlsx",3)[,1])
rownames(Ba_phylum) <- Ba_phylum_id$ASV_ID

Ba_phy <- data.frame(merge(Ba_phylum,asv_table,by = 0)[,-1])
Ba_phy$name <- gsub("^.{0,3}", "", Ba_phy$Phylum)

Ba_phy_sum <- aggregate(RA ~ name + DR, Ba_phy, sum)
Ba_phy_sum <- spread(Ba_phy_sum,DR,RA)
Ba_phy_sum[is.na(Ba_phy_sum)] <- 0
Ba_phy_sum <- gather(Ba_phy_sum,DR,RA,Dominant:Rare)

Ba_phy_sum$DR <- factor(Ba_phy_sum$DR, levels = c("Rare","Dominant"),labels =c("Rare","Dominant"))

a1 <- ggplot(Ba_phy_sum, aes(fill=DR, y=RA/60, x=reorder(name, -RA))) + 
  geom_bar(position="stack", stat="identity")+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  scale_fill_manual(values =alpha(c('#A9D18E','#E99A77'),1))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  xlab('')+ylab('Relative abundance (%)')

a2 <- ggplot(Ba_phy_sum, aes(fill=DR, y=RA, x=reorder(name, -RA))) + 
  geom_bar(position="fill", stat="identity")+
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  scale_fill_manual(values =alpha(c('#A9D18E','#E99A77'),1))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y=element_text(angle=90,vjust=0.5,hjust=0.5))+
  xlab('')+ylab('Relative abundance (%)')

gridExtra::grid.arrange(a1, a2, nrow=2)
#export::graph2ppt(file="phy.pptx", width=7, height=8)
######pf in each phylum or class######
###Ba
Ba_phy <- read_excel("Data/004_asv_phy_class.xlsx",4)
Ba_phy$name <- gsub("^.{0,3}", "", Ba_phy$SampleID)

Ba_phy_sum <- aggregate(RA ~ name + DR + Class, Ba_phy, sum)
Ba_phy_sum <- spread(Ba_phy_sum,Class,RA)
Ba_phy_sum[is.na(Ba_phy_sum)] <- 0

Ba_phy_sum <- gather(Ba_phy_sum,PF,RA,EP1:EP9)
Ba_phy_sum <- subset(Ba_phy_sum,name=='Proteobacteria'|name=='Actinobacteriota'|name=='Acidobacteriota'|
                     name=='Patescibacteria'|name=='Bacteroidota'|name=='Chloroflexi')

Ba_phy_sum <- data.table(Ba_phy_sum)
weight <- Ba_phy_sum[, .(weight1=sum(RA)), by=.(name,DR)]
weight_sum <- Ba_phy_sum[, .(weight_sum=sum(RA)), by=.(name)]
per1 <- Ba_phy_sum[, .(per1=sum(RA)), by=.(name,DR,PF)]

Ba_phy_sum <- merge(merge(merge(Ba_phy_sum, weight), weight_sum,by ='name'), per1,by=c("name","DR","PF"))

Ba_phy_sum$weight <- Ba_phy_sum$weight1/Ba_phy_sum$weight_sum
Ba_phy_sum$per <- Ba_phy_sum$per1/Ba_phy_sum$weight1
Ba_phy_sum$weight_per <- Ba_phy_sum$weight*Ba_phy_sum$per

Ba_phy_sum$taxon <- 'Ba'

###Fu
Fu_class <- read_excel("Data/004_asv_phy_class.xlsx",5)
Fu_class$name <- gsub("^.{0,3}", "", Fu_class$SampleID)

Fu_class_sum <- aggregate(RA ~ name + DR + Class, Fu_class, sum)
Fu_class_sum <- spread(Fu_class_sum,Class,RA)
Fu_class_sum[is.na(Fu_class_sum)] <- 0

Fu_class_sum <- gather(Fu_class_sum,PF,RA,EP1:EP9)

Fu_class_sum <- subset(Fu_class_sum,name=='Eurotiomycetes'|name=='Agaricomycetes'|name=='Mortierellomycetes'|
                         name=='Sordariomycetes'|name=='Dothideomycetes'|name=='Pezizomycetes')

Fu_class_sum <- data.table(Fu_class_sum)
weight <- Fu_class_sum[, .(weight1=sum(RA)), by=.(name,DR)]
weight_sum <- Fu_class_sum[, .(weight_sum=sum(RA)), by=.(name)]
per1 <- Fu_class_sum[, .(per1=sum(RA)), by=.(name,DR,PF)]

Fu_class_sum <- merge(merge(merge(Fu_class_sum, weight), weight_sum,by ='name'), per1,by=c("name","DR","PF"))

Fu_class_sum$weight <- Fu_class_sum$weight1/Fu_class_sum$weight_sum
Fu_class_sum$per <- Fu_class_sum$per1/Fu_class_sum$weight1
Fu_class_sum$weight_per <- Fu_class_sum$weight*Fu_class_sum$per

Fu_class_sum$taxon <- 'Fu'

#taxa_sum <- rbind(Ba_phy_sum,Fu_class_sum)
#write.table(taxa_sum,'taxa_sum.txt')

Ba_phy_sum$name <- factor(Ba_phy_sum$name,levels = c('Proteobacteria','Actinobacteriota','Acidobacteriota','Patescibacteria','Bacteroidota','Chloroflexi'))

p1 <- ggplot(data = Ba_phy_sum, aes(x = name, fill = PF, alpha= DR)) +
  scale_alpha_manual(values = c(1, 1)) +
  geom_bar(data=Ba_phy_sum[Ba_phy_sum$DR=="Dominant",], aes(y=weight_per),stat="identity",width=0.5) +
  geom_bar(data=Ba_phy_sum[Ba_phy_sum$DR=="Rare",], aes(y=-weight_per), stat="identity", width=0.5) +
  scale_fill_manual(values = c("#A9D18E","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#E99A77"))+
  scale_y_continuous(breaks=seq(-1,1,0.2),labels=abs(seq(-1,1,0.2)),position = "right") + 
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(legend.position="right") +geom_hline(yintercept=0, linetype="dashed")+
  geom_rect(data = data.frame(xmin = -Inf, xmax = Inf,ymin = c(-Inf, 0),ymax = c(0, Inf),fill = c("Dominant","Rare")),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = .2,inherit.aes = FALSE)+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=270,vjust=0.5,hjust=0.5))


Fu_class_sum$name <- factor(Fu_class_sum$name,levels = c('Eurotiomycetes','Agaricomycetes','Mortierellomycetes','Sordariomycetes','Dothideomycetes','Pezizomycetes'))

p2 <- ggplot(data = Fu_class_sum, aes(x = name, fill = PF, alpha= DR)) +
  scale_alpha_manual(values = c(1, 1)) +
  geom_bar(data=Fu_class_sum[Fu_class_sum$DR=="Dominant",], aes(y=weight_per),stat="identity",width=0.5) +
  geom_bar(data=Fu_class_sum[Fu_class_sum$DR=="Rare",], aes(y=-weight_per), stat="identity", width=0.5) +
  scale_fill_manual(values = c("#A9D18E","#a6cee3", "#1f78b4", "#b2df8a", "#33a02c","#fb9a99","#e31a1c","#ff7f00","#cab2d6","#E99A77"))+
  scale_y_continuous(breaks=seq(-1,1,0.2),labels=abs(seq(-1,1,0.2)),position = "right") + 
  theme_bw()+theme(strip.background = element_blank(), strip.placement = "outside",panel.grid = element_blank())+
  theme(legend.position="right") +geom_hline(yintercept=0, linetype="dashed")+
  geom_rect(data = data.frame(xmin = -Inf, xmax = Inf,ymin = c(-Inf, 0),ymax = c(0, Inf),fill = c("Dominant","Rare")),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = .2,inherit.aes = FALSE)+
  theme(axis.text.x=element_text(angle=0,vjust=0.5,hjust=0.5),axis.text.y=element_text(angle=270,vjust=0.5,hjust=0.5))


gridExtra::grid.arrange(p1, p2, nrow=1)


#export::graph2ppt(file="phy.pptx", width=6.38, height=4.61)