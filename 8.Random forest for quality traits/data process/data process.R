otu <- read.table("otu_table_family.txt",header = TRUE,row.names = 1,sep = "\t")
group <- read.table("Group.txt",sep = "\t",header = FALSE)
otu <- t(otu)
otu <- otu[group$V1,]
otu <- cbind(group$V2,otu)
otu <- as.data.frame(otu)
otu$V1 <- as.factor(otu$V1)
otu <- aggregate(otu[,2:ncol(otu)],list(otu$V1),mean)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
colnames(otu) <- gsub(".*f__","",colnames(otu))
write.table(otu,file = "family.txt",sep = "\t")