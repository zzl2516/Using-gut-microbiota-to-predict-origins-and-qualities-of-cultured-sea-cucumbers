class <- read.table("alpha.txt",header = TRUE,sep = "\t",row.names = 1)
class$Group <- factor(class$Group,levels = c("XP","SZ","YT","QD","WFD","BCD","DY","JZ"),
                      labels = c("XP (n=33)","SZ (n=18)","YT (n=31)","QD (n=30)",
                                 "WFD (n=37)","BCD (n=32)","DY (n=30)","JZ (n=78)"))

fit <- aov(Faith_pd~Group,data = class)
Tukey <- TukeyHSD(fit)

library(multcomp)
par(mar=c(4,6,12,2))  
tuk<-glht(fit,linfct=mcp(Group="Tukey"))
res <- cld(tuk,alpah=0.05)
res