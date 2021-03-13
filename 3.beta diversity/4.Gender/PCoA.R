library(vegan)
library(ggplot2)
library(ape)

data <- read.csv("ASV.txt", head=TRUE,sep="\t",row.names = 1)

groups <- read.table("group.txt",sep = "\t",header = F,colClasses = c("character"))
data <- data[,which(names(data)%in%groups$V1)]

data <- t(data)
data[is.na(data)] <- 0
data <- vegdist(data)

cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#99999",
                "#ADD1E5")
pcoa<- pcoa(data, correction = "none", rn = NULL)
PC1 = pcoa$vectors[,1]
PC2 = pcoa$vectors[,2]
plotdata <- data.frame(rownames(pcoa$vectors),PC1,PC2)
colnames(plotdata) <-c("sample","PC1","PC2")
colnames(groups) <- c("sample","Group")
plotdata <- merge(plotdata,groups)
pc1 <-floor(pcoa$values$Relative_eig[1]*100)
pc2 <-floor(pcoa$values$Relative_eig[2]*100)
plotdata$Group <- factor(plotdata$Group,levels = c("M","F"),
                         labels = c("Male","Female"))

otu.adonis=adonis(data~Group,data = plotdata,distance = "bray")

p2<-ggplot(plotdata, aes(PC1, PC2)) +
  geom_point(aes(fill=Group),size=8,pch = 21)+
  scale_fill_manual(values=cbbPalette,name = "Gender")+
  xlab(paste("PC1 ( ",pc1,"%"," )",sep="")) + 
  ylab(paste("PC2 ( ",pc2,"%"," )",sep=""))+
  theme(text=element_text(size=30))+
  geom_vline(aes(xintercept = 0),linetype="dotted")+
  geom_hline(aes(yintercept = 0),linetype="dotted")+
  geom_text(aes(x = -0.4,y = 0.38,label = paste("PERMANOVA:\ndf = ",otu.adonis$aov.tab$Df[1],
                                               "\nR2 = 0.003",
                                               "\np-value = 0.832 ",sep = "")),
            size = 6) +
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=34),
        axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=34),
        axis.title.y=element_text(colour='black', size=34),
        axis.text=element_text(colour='black',size=28),
        legend.title=element_text(size = 28,face = "bold"),
        legend.text=element_text(size=24),
        legend.key=element_blank(),legend.position = c(0.08,0.92),
        legend.background = element_rect(colour = "black"),
        legend.key.height=unit(1,"cm")) +
  guides(fill = guide_legend(ncol = 1))

pdf("PCoA12.pdf",height=12,width=15)
p2
png(filename="PCoA12.png",res=600,height=7000,width=9000)
p2
dev.off()