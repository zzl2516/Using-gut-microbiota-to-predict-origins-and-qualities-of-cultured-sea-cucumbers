library(ape)
tree <- read.tree("rep_set.tre")
otu <- read.table("ASV.txt",header = TRUE,row.names = 1,sep = "\t")
XP <- otu[,1:33]
XP <- XP[which(rowSums(XP) > 0),]
XP <- XP[1:100,]
XP <- t(XP)

library(picante)
dis <- prune.sample(XP,tree)
dis <- cophenetic(dis)
XP <- XP[,colnames(dis)]

mntd <- ses.mntd(XP,dis,abundance.weighted = TRUE,null.model = "taxa.labels",
                 runs = 2)

write.csv(mntd,file = "NTI.csv")

mntd.obs <- comdistnt(XP, dis, abundance.weighted = TRUE)
mntd.obs <- t(as.vector(t(mntd.obs)))
mntd.rand <- replicate(2,as.vector(t(comdistnt(XP,taxaShuffle(dis),
                                                abundance.weighted = TRUE))))
mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 1, FUN = mean, 
                        na.rm = TRUE)
mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 1, FUN = sd, 
                      na.rm = TRUE)
mntd.obs.z <- as.vector((mntd.obs - mntd.rand.mean)/mntd.rand.sd)
beta.mntd <- data.frame(as.vector(mntd.obs), mntd.rand.mean, mntd.rand.sd,mntd.obs.z)


a <- rownames(XP)
X <- c()
for (i in 2:length(a)) {
X <- c(X,a[i:length(a)])    
}

Y <- c()
for (i in 1:length(a)) {
Y <- c(Y,rep(a[i],length(a)-i))
}

beta.mntd <- cbind(X,Y,beta.mntd)

write.csv(beta.mntd,file = "betaNTI.csv")