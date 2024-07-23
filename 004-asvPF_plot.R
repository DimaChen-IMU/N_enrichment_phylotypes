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

PF_sum_percentage <- PF_sum %>% dplyr::group_by(DR,taxon) %>% dplyr::mutate(count100 = count / sum(count)) %>% dplyr::mutate(RA100 = RA / sum(RA)) 
write.table(PF_sum_percentage,'PF_sum_percentage.txt')
#EP1
EP1 <- gather(subset(PF_sum_percentage,Class =='EP1'), key="Type", value="per",  count100:RA100)

ep1.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP1", gp=gpar(fontsize=28)))

ep1.2 <- ggplot(subset(EP1, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))
  

ep1.3 <- ggplot(subset(EP1, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep1.4 <- ggplot(subset(EP1, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

#layout <- rbind(c(1, 1, 2, 2),c(1, 1, 2, 2),c(3, 3, 2, 2),c(3, 3, 4, 4))

layout <- rbind(c(1, 2), c(3, 4))

ep1 <- grid.arrange(ep1.1, ep1.2, ep1.3, ep1.4, layout_matrix=layout)

#EP2
EP2 <- gather(subset(PF_sum_percentage,Class =='EP2'), key="Type", value="per",  count100:RA100)

ep2.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP2", gp=gpar(fontsize=28)))

ep2.2 <- ggplot(subset(EP2, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep2.3 <- ggplot(subset(EP2, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep2.4 <- ggplot(subset(EP2, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep2 <- grid.arrange(ep2.1, ep2.2, ep2.3, ep2.4, layout_matrix=layout)

#EP3
EP3 <- gather(subset(PF_sum_percentage,Class =='EP3'), key="Type", value="per",  count100:RA100)

ep3.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP3", gp=gpar(fontsize=28)))

ep3.2 <- ggplot(subset(EP3, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep3.3 <- ggplot(subset(EP3, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep3.4 <- ggplot(subset(EP3, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep3 <- grid.arrange(ep3.1, ep3.2, ep3.3, ep3.4, layout_matrix=layout)

#EP4
EP4 <- gather(subset(PF_sum_percentage,Class =='EP4'), key="Type", value="per",  count100:RA100)

ep4.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP4", gp=gpar(fontsize=28)))

ep4.2 <- ggplot(subset(EP4, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep4.3 <- ggplot(subset(EP4, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep4.4 <- ggplot(subset(EP4, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep4 <- grid.arrange(ep4.1, ep4.2, ep4.3, ep4.4, layout_matrix=layout)

#EP5
EP5 <- gather(subset(PF_sum_percentage,Class =='EP5'), key="Type", value="per",  count100:RA100)

ep5.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP5", gp=gpar(fontsize=28)))

ep5.2 <- ggplot(subset(EP5, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 1),expand = expansion(mult = c(0, .3))) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep5.3 <- ggplot(subset(EP5, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 1),expand = expansion(mult = c(0, .2))) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep5.4 <- ggplot(subset(EP5, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 1),expand = expansion(mult = c(0, .2))) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep5 <- grid.arrange(ep5.1, ep5.2, ep5.3, ep5.4, layout_matrix=layout)

#EP6
EP6 <- gather(subset(PF_sum_percentage,Class =='EP6'), key="Type", value="per",  count100:RA100)

ep6.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP6", gp=gpar(fontsize=28)))

ep6.2 <- ggplot(subset(EP6, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep6.3 <- ggplot(subset(EP6, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep6.4 <- ggplot(subset(EP6, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep6 <- grid.arrange(ep6.1, ep6.2, ep6.3, ep6.4, layout_matrix=layout)


#EP7
EP7 <- gather(subset(PF_sum_percentage,Class =='EP7'), key="Type", value="per",  count100:RA100)

ep7.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP7", gp=gpar(fontsize=28)))

ep7.2 <- ggplot(subset(EP7, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep7.3 <- ggplot(subset(EP7, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep7.4 <- ggplot(subset(EP7, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep7 <- grid.arrange(ep7.1, ep7.2, ep7.3, ep7.4, layout_matrix=layout)

#EP8
EP8 <- gather(subset(PF_sum_percentage,Class =='EP8'), key="Type", value="per",  count100:RA100)

ep8.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP8", gp=gpar(fontsize=28)))

ep8.2 <- ggplot(subset(EP8, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep8.3 <- ggplot(subset(EP8, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep8.4 <- ggplot(subset(EP8, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep8 <- grid.arrange(ep8.1, ep8.2, ep8.3, ep8.4, layout_matrix=layout)


#EP9
EP9 <- gather(subset(PF_sum_percentage,Class =='EP9'), key="Type", value="per",  count100:RA100)

ep9.1 <- grobTree(rectGrob(gp=gpar(fill=NA,color=NA)), textGrob("EP9", gp=gpar(fontsize=28)))

ep9.2 <- ggplot(subset(EP9, DR=='All'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep9.3 <- ggplot(subset(EP9, DR=='Dominant'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 1),expand = expansion(mult = c(0, .2))) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

ep9.4 <- ggplot(subset(EP9, DR=='Rare'), aes(x = interaction( Type,taxon), y=per, fill=interaction( Type,taxon)))+
  geom_bar(stat="identity",width = 0.5)+scale_fill_manual(values = c("#EA4C46", "#F07470", "#088CFF", "#64BAFE")) +
  geom_text(aes(label=sprintf("%0.2f", round(per*100, digits = 2))), position=position_dodge(width=0.9), vjust=0.25,hjust=0,angle = 90,size=2.2)+
  scale_y_continuous(limits = c(0, 0.5),expand = c(0, 0)) +theme_no_axes()+
  theme(legend.position = "none",panel.border = element_blank(),plot.margin = unit(c(0,0,0,0), "cm"))

layout <- rbind(c(1, 2), c(3, 4))

ep9 <- grid.arrange(ep9.1, ep9.2, ep9.3, ep9.4, layout_matrix=layout)

grid.arrange(ep1, ep4, ep7, ep2,ep5, ep8, ep3, ep6,ep9, nrow=3)

#export::graph2ppt(file="BBB.pptx", width=4.5, height=3)