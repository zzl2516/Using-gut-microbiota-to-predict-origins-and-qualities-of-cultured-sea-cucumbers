library(ggplot2)

data <- read.table("important.txt",header = TRUE,sep = "\t")

library(reshape)
data <- melt(data)

data$X <- factor(data$X,levels = c("Proteobacteria","Bacteroidetes","Actinobacteria",
                                   "Verrucomicrobia","Cyanobacteria","Firmicutes",
                                   "Chloroflexi","Acidobacteria","Planctomycetes",
                                   "TM7","TM6","Tenericutes","Chlamydiae",
                                   "Gemmatimonadetes","WPS-2","SBR1093",
                                   "Spirochaetes","WS6","GN02","BRC1","WWE1","WS2"))


p1 <- ggplot(data,aes(variable, weight=value, fill=X)) +
    geom_bar(position = "stack",width = .7,color = "black") +
    theme_bw() +
    theme(panel.grid=element_blank()) +
    labs(x="",y = "Number of orders") +
    theme(legend.title = element_blank()) +
    theme(axis.text.x=element_text(colour="black",size=15,angle = 45,hjust = 1)) +
    theme(axis.text.y=element_text(colour = "black",size = 15)) +
    theme(axis.line.x=element_line(colour = "black")) +
    theme(axis.line.y=element_line(colour = "black")) +
    theme(axis.title.y= element_text(size = 20)) +
    theme(legend.text = element_text(colour = "black",size = 15)) +
    scale_y_continuous(limits = c(0,15),expand = c(0,0)) +
    theme(legend.position = "bottom") +
    theme(panel.border = element_blank())

png(filename = "important.png",width = 5600,height = 4400,res = 600)
p1
dev.off()

