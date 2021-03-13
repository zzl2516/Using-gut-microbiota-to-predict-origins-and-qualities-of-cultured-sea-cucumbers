class <- read.table("nutrition.txt",header = TRUE,sep = "\t")

colnames(class) <- c("Group","Protein (%)","Fat (%)","Sugar (%)",
                   "Saponin (g/kg)","Collagen (g/kg)","VA (μg/kg)","VE (mg/kg)",
                   "Taurine (g/kg)","CS (mg/kg)")

library(dplyr)
library(ggplot2)
library(reshape2)

bb <- class
bb1 <- melt(class)
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")
bb$Group <- factor(bb$Group,levels = c("BCD","WFD","JZ","SZ","DY","YT","QD","XP"),
                   labels = c("DLE","DLW","JZ","QHD","DY","YT","QD","XP"))
bb1$Group <- factor(bb1$Group,levels = c("BCD","WFD","JZ","SZ","DY","YT","QD","XP"),
                    labels = c("DLE","DLW","JZ","QHD","DY","YT","QD","XP"))

library(multcomp)
bb.sample <- colnames(bb)[2:ncol(bb)]
test.b <- c()
for (i in bb.sample) {
    fit1 <- aov(as.formula(sprintf("`%s` ~ Group",i)),
                data = bb)
    tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
    res1 <- cld(tuk1,alpah=0.05)
    test.b <- cbind(test.b,res1$mcletters$Letters)
}
colnames(test.b) <- colnames(bb)[2:ncol(bb)]
test.b <- melt(test.b)
colnames(test.b) <- c("Group","variable","value")

library(tidyverse)
test.b1 <- bb %>% gather(variable,value,-Group) %>% group_by(variable,Group) %>% 
    summarise(Max = max(value))
test.b11 <- dcast(test.b1,Group~variable)
for (i in 2:ncol(test.b11)) {
    test.b11[,i] <- test.b11[,i] + max(test.b11[,i])*0.05
}
test.b11 <- melt(test.b11)
test.b1 <- merge(test.b1,test.b11,by = c("variable","Group"))
test.b2 <- merge(test.b,test.b1,by = c("variable","Group"))


png("nutrition_anova.png",width = 6400,height = 5000,res = 600)
ggplot(bb1,aes(x = Group,y = value,color = Group)) +
    geom_boxplot(color = "black",outlier.colour = NA) +
    geom_jitter(position = position_jitter(0.2)) +
    facet_wrap(.~variable,scales = "free_y",ncol = 3) +
    scale_color_manual(values=cbbPalette,guide = FALSE) +
    geom_text(data = test.b2,aes(x = Group,y = value.y,label = value.x),
              size = 5,color = "black",fontface = "bold") +
    ylab("The values of nutrients in body wall of sea cucumbers") +
    theme_bw()+
    theme(axis.ticks.length = unit(0.4,"lines"), 
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"), 
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour='black', size=20,face = "bold",vjust = 3),
          axis.text.y=element_text(colour='black',size=12),
          axis.text.x=element_text(colour = "black",size = 14,face = "bold",
                                   angle = 45,hjust = 1,vjust = 1),
          plot.margin = margin(t = 5,r = 5,b = 5, l = 20, unit = "pt"),
          text = element_text(colour = "black",size = 20,face = "bold"),
          legend.position = "none")
dev.off()

