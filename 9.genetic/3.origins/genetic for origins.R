library(ape)
library(vegan)
library(reshape2)
library(ggplot2)
library(phangorn)
seq <- read.FASTA("align.fas")
names(seq) <- gsub("-","_",names(seq))
sample <- read.table("sample.txt",header = FALSE,sep = "\t")
comm <- read.table("ASV.txt",header = TRUE,row.names = 1,sep = "\t")
comm <- comm[,as.vector(sample$V1)]
group <- read.table("group.txt",header = FALSE,sep = "\t")
group <- merge(group,sample,by = "V1")
comm <- t(comm)

## XP
g <- group$V1[group$V2 == "DY"]
comm.1 <- comm[g,]
comm.1 <- comm.1[,colSums(comm.1)>0]
comm.1[is.na(comm.1)] <- 0
dis <- vegdist(comm.1)

seq.1 <- seq[g]
mat = dist.dna(seq.1,as.matrix = TRUE)
bc <- hclust(dis,method = "average")

mantel(mat,dis)
mantel(mat,dis,method = "spearman")

diag(mat) <- 0
mat[upper.tri(mat)] <- NA
mat.r <- melt(mat)
mat.r <- mat.r[-which(mat.r$Var1 == mat.r$Var2),]
mat.r <- na.omit(mat.r)

dis <- as.matrix(dis)
diag(dis) <- 0
dis[upper.tri(dis)] <- NA
dis.r <- melt(dis)
dis.r <- dis.r[-which(dis.r$Var1 == dis.r$Var2),]
dis.r <- na.omit(dis.r)

comm_dis <- cbind(mat.r,dis.r)
colnames(comm_dis) <- c("site1x","site1y","gene_dis","site2x","site2y","comm_dis")

fit1 <- lm(comm_dis~gene_dis,data = comm_dis)
summary(fit1)

lm_comm_gene <- ggplot(comm_dis,aes(gene_dis,comm_dis)) +
    geom_point() +
    geom_smooth(method = "lm",se = FALSE) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(colour = "black",fill = "transparent"),
          plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15)) +
    scale_y_continuous(limits = c(0,1)) +
    labs(y = "Bray-Curtis dissimilarity",x = "Genetic distance") +

png(filename="dis_comm_gene_DY.png",res=600,height=3000,width=3600)
lm_comm_gene
dev.off()

tr <- bionj(mat)
bc<- as.phylo(bc)

RF.dist(bc,tr,normalize = TRUE)

results <- replicate(999,RF.dist(rtree(nrow(comm.1),tip.label = bc$tip.label),tr,
                                 normalize = TRUE))
length(results[results > RF.dist(bc,tr,normalize = TRUE)])/1000


