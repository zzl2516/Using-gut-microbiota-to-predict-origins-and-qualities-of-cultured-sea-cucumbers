library(ggplot2)
library(reshape2)
library(dplyr)
library(vegan)

abundance <- read.table("genus.txt",header = TRUE,sep = "\t")
importance <- read.table("importance.txt",header = TRUE,sep = "\t")

abundance <- abundance[abundance$.OTU.ID %in% importance$feature,]
abundance <- abundance[-8,]

abundance[,2:9] <- decostand(abundance[,2:9],"standardize",1)

abundance <- melt(abundance)

abundance$variable <- factor(abundance$variable,levels = c("BCD","WFD","JZ","SZ","DY","YT","QD","XP"),
                             labels = c("DLE (n=32)","DLW (n=37)","JZ (n=78)","QHD (n=18)",
                                        "DY (n=30)","YT (n=31)","QD (n=30)","XP (n=33)"))
abundance$.OTU.ID <- factor(abundance$.OTU.ID,
                            levels = c("Ruegeria","Plesiocystis","Amaricoccus",
                                       "Luteolibacter","Paulinella","Bacillus",
                                       "Dinoroseobacter","Persicirhabdus",
                                       "Synechococcus","Shimia","Formosa",
                                       "Haloferula","Rubritalea","Pseudovibrio",
                                       "Roseobacter","Cupriavidus","Aureitalea",
                                       "Lutimonas","Desulfococcus","Actibacter"))

importance$feature <- factor(importance$feature,
                            levels = c("Actibacter","Desulfococcus","Lutimonas",
                                       "Aureitalea","Cupriavidus","Roseobacter",
                                       "Pseudovibrio","Rubritalea","Haloferula",
                                       "Formosa","Shimia","Synechococcus",
                                       "Persicirhabdus","Dinoroseobacter",
                                       "Bacillus","Paulinella","Luteolibacter",
                                       "Amaricoccus","Plesiocystis","Ruegeria"))


p1 <- ggplot(abundance) +
    geom_tile(aes(variable,.OTU.ID,fill = value)) +
    theme_classic() +
    theme(axis.ticks = element_blank(),axis.line = element_blank()) +
    xlab("") + ylab("") +
    scale_fill_gradient2('Abundance',low = 'white', high = 'red') +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12,face = "bold",angle = 45,
                                     vjust = 1,hjust = 1,colour = "black"),
          axis.text.y = element_text(size = 14,face = "bold.italic",
                                     colour = "black"))

p2<- ggplot(importance,aes(feature,importance)) +
    geom_bar(stat = "identity",width = 0.7) + coord_flip() +
    theme_bw()+ theme(panel.grid=element_blank()) + 
    theme(panel.border = element_blank()) +
    theme(panel.background=element_rect(fill='transparent', color='black'),
          plot.margin = unit(c(3,5,1,1),"mm")) +
    scale_y_continuous(limits = c(0,0.05),expand = c(0,0)) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(size = 12,face = "bold",colour = "black"),
          axis.title.x = element_text(size = 16,face = "bold",vjust = 10,
                                      colour = "black")) +
    xlab("") +ylab("Importance")

library(patchwork)
p3 <- p1 + p2 

png(filename="importance.png",res=600,height=4800,width=8000)
p3
dev.off()

