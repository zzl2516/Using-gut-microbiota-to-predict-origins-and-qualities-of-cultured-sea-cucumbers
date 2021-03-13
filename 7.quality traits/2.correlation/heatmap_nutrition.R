otu <- read.table("otu_table_genus.txt",header = TRUE,row.names = 1,sep = "\t")
group <- read.table("Group.txt",sep = "\t",header = FALSE)
otu <- t(otu)
otu <- otu[group$V1,]
otu <- cbind(group$V2,otu)
otu <- as.data.frame(otu)
otu$V1 <- as.factor(otu$V1)
otu <- aggregate(otu[,2:ncol(otu)],list(otu$V1),mean)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
colnames(otu) <- gsub(".*g__","",colnames(otu))

marker <- c("Actibacter","Desulfococcus","Lutimonas",
            "Aureitalea","Cupriavidus","Roseobacter",
            "Pseudovibrio","Rubritalea","Haloferula",
            "Formosa","Shimia","Synechococcus",
            "Persicirhabdus","Dinoroseobacter",
            "Bacillus","Paulinella","Luteolibacter",
            "Amaricoccus","Plesiocystis","Ruegeria")

genus <- otu[,marker]

nut <- read.table("nutrition.txt",header = TRUE,row.names = 1,sep = "\t")

library(psych)
library(pheatmap)
res=corr.test(genus,nut,method = "spearman",adjust = "holm")
write.table(res$r,"correlation.xls",sep="\t",quote=F,col.names=NA)
write.table(res$p,"pvalue.xls",sep="\t",quote=F,col.names=NA)

png(filename="correlation.png",res=600,height=7200,width=5400)
pheatmap(res$r,display_numbers = matrix(ifelse(res$p <= 0.01, "++", ifelse(res$p <= 0.05 ,"+"," ")), nrow(res$p)), 
         fontsize_number=30,fontsize = 20,cluster_rows = FALSE,cluster_cols = FALSE,fontface = "bold",angle_col = 90,
         number_color = "white")
dev.off()

