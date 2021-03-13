library(ggtree)
library(ggplot2)
data <- read.table("weighted_unifrac_distance.txt",header = TRUE,
                   row.names = 1,sep = "\t")
group <- read.table("group.txt",sep = "\t")
colnames(group) <- c("ID","Group")
data <- as.dist(data)
hc <- hclust(data, method = "average")

cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#000000")
group$Group <- factor(group$Group,levels = c("DLE","DLW","JZ","QHD","DY","YT","QD","XP"),
                      labels = c("DLE (n=32)","DLW (n=37)","JZ (n=78)","QHD (n=18)",
                                 "DY (n=30)","YT (n=31)","QD (n=30)","XP (n=33)"))


p <- ggtree(hc,layout = "circular") %<+% group +
    geom_tiplab(size = 3, aes(angle = angle,color = Group,x = x*1.02)) +
    geom_tippoint(aes(color = Group),size = 3) +
    theme(legend.position = c(0.09,0.86),
          legend.title = element_text(face = "bold",size = 15),
          legend.text = element_text(size = 12)) +
    scale_color_manual(values=cbbPalette)

png(filename="UPGMA.png",res=600,height=7000,width=9000)
p
dev.off()

rownames(group) <- group$ID
group <- as.data.frame(group[,-1])
gheatmap(p2,group,offset = 15,width = .3) +
    scale_fill_viridis_d(option = "D",values = cbbPalette)


    

