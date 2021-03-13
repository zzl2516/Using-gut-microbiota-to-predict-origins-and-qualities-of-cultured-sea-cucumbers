library(ape)
seq <- read.FASTA("align.fas")
mat <- dist.dna(seq,as.matrix = TRUE)
rownames(mat) <- gsub("-","_",rownames(mat))
colnames(mat) <- gsub("-","_",rownames(mat))

tr <- bionj(mat)

group <- read.table("group.txt",sep = "\t")
colnames(group) <- c("ID","Group")
group$Group <- factor(group$Group,levels = c("DLE","DLW","JZ","QHD","DY","YT","QD","XP"),
                      labels = c("DLE (n=32)","DLW (n=37)","JZ (n=78)","QHD (n=18)",
                                 "DY (n=30)","YT (n=31)","QD (n=30)","XP (n=33)"))
cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
                "#D55E00","#CC79A7","#000000")

library(ggtree)
library(ggplot2)
p <- ggtree(tr,layout = "fan",branch.length = "none") %<+% group +
    geom_tiplab(size = 3, aes(angle = angle,color = Group,x = x*1.02)) +
    geom_tippoint(aes(color = Group),size = 3) +
    theme(legend.position = c(0.09,0.86),
          legend.title = element_text(face = "bold",size = 15),
          legend.text = element_text(size = 12)) +
    scale_color_manual(values=cbbPalette)

png(filename="tree.png",res=600,height=7000,width=9000)
p
dev.off()