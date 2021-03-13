class <- read.table("alpha.txt",header = TRUE,sep = "\t",row.names = 1)

library(dplyr)
library(ggplot2)
library(reshape2)

class <- melt(class)
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")
class$Group <- factor(class$Group,levels = c("BCD","WFD","JZ","SZ","DY","YT","QD","XP"),
                      labels = c("DLE (n=32)","DLW (n=37)","JZ (n=78)","QHD (n=18)",
                                 "DY (n=30)","YT (n=31)","QD (n=30)","XP (n=33)"))
class$variable <- factor(class$variable, levels = c("Chao1","Observed_species",
                                                    "Shannon","Simpson",
                                                    "Pielou_e","Faith_pd"))

yf <- class
yd <- yf %>% group_by(variable,Group,add =TRUE) %>% summarise(Max = max(value))
yl <- yd %>% group_by(variable) %>% summarise(Max = max(Max))
ya <- data.frame()
for (i in 1:nrow(yl)) {
    ya <- rbind(ya,yl[rep(i,8),])
}
yd$Max <- yd$Max + ya$Max*0.05

test <- read.table("test.txt",header = TRUE,sep = "\t")
test <- melt(test,id.vars = "Group")
test$Group <- factor(test$Group,levels = c("BCD","WFD","JZ","SZ","DY","YT","QD","XP"),
                     labels = c("DLE (n=32)","DLW (n=37)","JZ (n=78)","QHD (n=18)",
                                "DY (n=30)","YT (n=31)","QD (n=30)","XP (n=33)"))
test$variable <- as.factor(test$variable)
test <- merge(test,yd)

png("alpha_anova.png",width = 5000,height = 6000,res = 600)
ggplot(class,aes(x = Group,y = value,color = Group)) +
    geom_boxplot(color = "black",outlier.colour = NA) +
    geom_jitter(position = position_jitter(0.2)) +
    facet_wrap(.~variable,scales = "free_y",ncol = 2) +
    scale_color_manual(values=cbbPalette,guide = FALSE) +
    geom_text(data = test,aes(x = Group,y = Max,label = value),
              size = 5,color = "black",fontface = "bold") +
    ylab("The alpha indices of gut microbiota among sea cucumbers from different origins") +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour='black', size=16,face = "bold",vjust = 3),
          axis.text.y=element_text(colour='black',size=12),
          axis.text.x=element_text(colour = "black",size = 14,face = "bold",
                                   angle = 45,hjust = 1,vjust = 1),
          plot.margin = margin(t = 5,r = 5,b = 5, l = 20, unit = "pt"),
          text = element_text(colour = "black",size = 20,face = "bold"),
          legend.position = "none")
dev.off()

