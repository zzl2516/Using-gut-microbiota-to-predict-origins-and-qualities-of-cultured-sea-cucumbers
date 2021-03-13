library(ape)
library(vegan)
library(reshape2)
library(ggplot2)
seq <- read.FASTA("align.fas")

mat = dist.dna(seq,as.matrix = TRUE)
rownames(mat) <- gsub("-","_",rownames(mat))
colnames(mat) <- gsub("-","_",rownames(mat))

sample <- read.table("sample.txt",header = FALSE,sep = "\t")

comm <- read.table("ASV.txt",header = TRUE,row.names = 1,sep = "\t")
comm <- comm[,as.vector(sample$V1)]

comm <- t(comm)
rownames(comm) <- gsub("-","_",rownames(comm))
comm[is.na(comm)] <- 0
dis <- vegdist(comm)

bc <- hclust(dis,method = "average")

mantel.p <- mantel(mat,dis)
mantel.s <- mantel(mat,dis,method = "spearman")

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

label_fit <- data.frame(
    formula = sprintf('italic(Y) == %.7f*italic(X) + %.3f', coefficients(fit1)[2], 
                      coefficients(fit1)[1]),
    R2 = sprintf('R^2 == %.8f', 0.00164),
    p_value = sprintf('italic(P-value) == %.3f', 0.369))

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
    geom_text(x = 0.017, y = 0.2, aes(label = R2), data = label_fit, parse = TRUE,
              hjust = 0, color = 'black', show.legend = FALSE,fontface = "bold",
              size = 6) +
    geom_text(x = 0.017, y = 0.1, aes(label = p_value), data = label_fit, parse = TRUE,
              hjust = 0, color = 'black', show.legend = FALSE,fontface = "bold",
              size = 6)

png(filename="dis_comm_gene.png",res=600,height=3000,width=3600)
lm_comm_gene
dev.off()

tr <- bionj(mat)
bc<- as.phylo(bc)

library(phangorn)
RF.dist(bc,tr,normalize = TRUE)
## [1] 0.9957806

results <- replicate(999,RF.dist(rtree(240,tip.label = bc$tip.label),tr,
                                 normalize = TRUE))
p.value <- length(results[results > RF.dist(bc,tr,normalize = TRUE)])/1000
## 0.957

group <- read.table("group.txt",header = FALSE,sep = "\t")
group <- merge(group,sample,by = "V1")


mat.r$Site1 <- group$V2[match(mat.r$Var1,group$V1)]
mat.r$Site2 <- group$V2[match(mat.r$Var2,group$V1)]

dis.r$Site1 <- group$V2[match(dis.r$Var1,group$V1)]
dis.r$Site2 <- group$V2[match(dis.r$Var2,group$V1)]

library(tidyverse)
mat.s <- mat.r %>%
    group_by(Site1,Site2) %>%
    summarise(Mean = mean(value))
mat.s <- mat.s[-which(mat.s$Site1 == mat.s$Site2),]

dis.s <- dis.r %>%
    group_by(Site1,Site2) %>%
    summarise(Mean = mean(value))
dis.s <- dis.s[-which(dis.s$Site1 == dis.s$Site2),]

comm_dis_s <- cbind(mat.s,dis.s)
colnames(comm_dis_s) <- c("site1x","site1y","gene_dis","site2x","site2y","comm_dis")

fit2 <- lm(comm_dis~gene_dis,data = comm_dis_s)
summary(fit2)

label_fit <- data.frame(
    formula = sprintf('italic(Y) == %.7f*italic(X) + %.3f', coefficients(fit2)[2], 
                      coefficients(fit2)[1]),
    R2 = sprintf('R^2 == %.10f', 0.029),
    p_value = sprintf('italic(P-value) == %.3f', 0.383))

lm_comm_gene <- ggplot(comm_dis_s,aes(gene_dis,comm_dis)) +
    geom_point() +
    geom_smooth(method = "lm",se = FALSE) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(colour = "black",fill = "transparent"),
          plot.title = element_text(hjust = 0.5,face = "bold",size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15)) +
    scale_y_continuous(limits = c(0.6,1)) +
    labs(y = "Bray-Curtis dissimilarity",x = "Genetic distance") +
    geom_text(x = 0.0063, y = 0.63, aes(label = R2), data = label_fit, parse = TRUE,
              hjust = 0, color = 'black', show.legend = FALSE,fontface = "bold",
              size = 6) +
    geom_text(x = 0.0063, y = 0.6, aes(label = p_value), data = label_fit, parse = TRUE,
              hjust = 0, color = 'black', show.legend = FALSE,fontface = "bold",
              size = 6)

png(filename="dis_comm_gene_s.png",res=600,height=3000,width=3600)
lm_comm_gene
dev.off()


mat.s <- mat.r %>%
    group_by(Site1,Site2) %>%
    summarise(Mean = mean(value))
mat.s.m <- acast(mat.s,Site1~Site2)

mat.s.m <- ifelse(is.na(mat.s.m),t(mat.s.m),mat.s.m)
diag(mat.s.m) <- 0
mat.s.mt <- as.dist(mat.s.m)
mat.s.t <- hclust(mat.s.mt,method = "average")
mat.s.mt <- as.phylo(mat.s.t)

dis.s <- dis.r %>%
    group_by(Site1,Site2) %>%
    summarise(Mean = mean(value))
dis.s.m <- acast(dis.s,Site1~Site2)

dis.s.m <- ifelse(is.na(dis.s.m),t(dis.s.m),dis.s.m)
diag(dis.s.m) <- 0
dis.s.mt <- as.dist(dis.s.m)
dis.s.t <- hclust(dis.s.mt,method = "average")
dis.s.mt <- as.phylo(dis.s.t)

RF.dist(dis.s.mt,mat.s.mt,normalize = TRUE)
## 1
dist.topo(dis.s.mt,mat.s.mt,method = "score")
## 0.0684765

results <- replicate(999,RF.dist(rtree(8,tip.label = dis.s.mt$tip.label),mat.s.mt,
                                 normalize = TRUE))
p.value <- length(results[results > RF.dist(dis.s.mt,mat.s.mt,normalize = TRUE)])/1000
## 0.999

results <- replicate(999,dist.topo(rtree(8,tip.label = dis.s.mt$tip.label),mat.s.mt,
                                 method = "score"))
p.value <- length(results[results > dist.topo(dis.s.mt,mat.s.mt,method = "score")])/1000
