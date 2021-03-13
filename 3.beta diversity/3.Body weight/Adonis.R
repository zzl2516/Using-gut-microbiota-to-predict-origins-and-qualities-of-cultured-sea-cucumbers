library(vegan)
library(permute)     
library(lattice)

data.sp <- read.csv("zoo.sp.txt", head=TRUE,sep="\t")
data.su <- read.csv("zoo.su.txt", head=TRUE,sep="\t")
data.au <- read.csv("zoo.au.txt", head=TRUE,sep="\t")
data1 <- merge(data.sp,data.su,all = TRUE)
data <- merge(data1,data.au,all = TRUE)
row.names(data) <- data[,1]
data <- data[,-1]
data[is.na(data)] <- 0
data <- t(data)
group=read.table("group.txt",sep = "\t",row.names = 1)

otu.adonis=adonis(data~V2,data = group,distance = "bray")              
otu.adonis
