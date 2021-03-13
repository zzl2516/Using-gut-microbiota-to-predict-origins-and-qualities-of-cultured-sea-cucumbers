tax <- read.table("tax.txt",header = TRUE,row.names = 1,sep = "\t")
tax$Phylum <- gsub(".*p__","",tax$taxonomy)
tax$Phylum <- gsub("\\; c__.*","",tax$Phylum)
tax$Class <- gsub(".*c__","",tax$taxonomy)
tax$Class <- gsub("\\; o__.*","",tax$Class)
tax$Order <- gsub(".*o__","",tax$taxonomy)
tax$Order <- gsub("\\; f__.*","",tax$Order)
tax$Family <- gsub(".*f__","",tax$taxonomy)
tax$Family <- gsub("\\; g__.*","",tax$Family)
tax$Genus <- gsub(".*g__","",tax$taxonomy)
tax$Genus <- gsub("\\; s__.*","",tax$Genus)
tax$Species <- gsub(".*s__","",tax$taxonomy)
tax <- subset(tax, select = -taxonomy)

otu <- read.delim("ASV.txt",row.names = 1)
map <- read.delim("map.txt",row.names = 1)
library(ape)
library(tidyverse)
library(phyloseq)
tree <- read.tree("rep_set.tre")

inputMicro = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group"){
    
    if (is.null(otu)&is.null(tax)&is.null(map)) {
        ps = ps
        map = as.data.frame(sample_data(ps))
        map = map[, group]
        colnames(map) = "Group"
        map$Group = as.factor(map$Group)
        sample_data(ps) = map
        map = NULL
    }
    
    if (is.null(ps) ) {
        
        if (!is.null(otu)) {
            head(otu)
            otu = as.matrix(otu)
            str(otu)
            
            ps <- phyloseq(otu_table(otu, taxa_are_rows=TRUE))
            
        }  
        
        if (!is.null(tax) ) {
            head(tax)
            tax = as.matrix(tax)
            # taxa_names(tax)
            x = tax_table(tax)
            ps = merge_phyloseq(ps,x)
            ps
        }
        
        
        if (!is.null(map) ){
            
            map = map[group]
            
            map[,group] = as.factor(map[,group] )
            map$Group 
            z  = sample_data(map)
            ps = merge_phyloseq(ps,z)
            ps
        }
        if (!is.null(tree) ) {
            # #导入进化树
            h = phy_tree(tree)
            ps = merge_phyloseq(ps,h)
            ps
        }
        
        
    }
    return(ps)
    
}

ps <- inputMicro(otu,tax,map,tree,group = "Treat1")

##中性理论模型

library(tidyverse)
library(vegan)
library(picante)
library(minpack.lm)
library(FSA)
library(eulerr)
library(ggplot2)
library(grid)
require(Hmisc)
require(stats4)
library(ggpubr)

neutralModel = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group"){
    
    # 抽平，默认使用最小序列抽平
    ps = inputMicro(otu,tax,map,tree,ps,group  = group)
    ps
    set.seed(72)  #设置随机种子，保证结果可重复
    psrare = rarefy_even_depth(ps)
    
    # 标准化
    ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
    
    
    
    #------------------------------------------开始计算中性模型----------------------------------------------------------
    map = as.data.frame(sample_data(psrare))
    aa = levels(map$Group)
    aa
    map$ID = row.names(map)
    # i = 17
    # i = 16
    # i = 4
    # aa[i]
    plots = list()
    i =1
    for (i in 1:length(aa)) {
        
        
        maps<- dplyr::filter(map,Group %in%aa[i])
        row.names(maps) = maps$ID
        ps_sub = psrare
        sample_data( ps_sub ) = maps;ps_sub
        
        # 提取OTU表格
        OTU.table = t(otu_table(ps_sub))
        head(OTU.table )
        # 将整个群落看做一个整体，计算每个样本的序列数，并求取均值Calculate the number of individuals in the meta community (Average read depth)
        N <- mean(apply(OTU.table, 1, sum))
        
        #计算每个OTU的的平均序列数 Calculate the average relative abundance of each taxa across communities
        p.m <- apply(OTU.table, 2, mean)
        #去除OTU序列数为0的OTU
        p.m <- p.m[p.m != 0]
        p <- p.m/N
        p.df = data.frame(p) %>%
            rownames_to_column(var="OTU")
        
        # Calculate the occurrence frequency of each taxa
        OTU.table.bi <- 1*(OTU.table>0)
        freq.table <- apply(OTU.table.bi, 2, mean)
        freq.table <- freq.table[freq.table != 0]
        freq.df = data.frame(OTU=names(freq.table), freq=freq.table)
        
        #Combine
        C <- inner_join(p.df,freq.df, by="OTU") %>%
            arrange(p)
        # Remove rows with any zero (absent in either source pool or local communities). You already did this, but just to make sure we will do it again.
        C.no0 <- C %>%
            filter(freq != 0, p != 0)
        
        #Calculate the limit of detection
        d <- 1/N
        
        ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
        p.list <- C.no0$p
        freq.list <- C.no0$freq
        m.fit <- nlsLM(freq.list ~ pbeta(d, N*m*p.list, N*m*(1-p.list), lower.tail=FALSE), start=list(m=0.1))
        m.ci <- confint(m.fit, 'm', level=0.95)
        m.sum <- summary(m.fit)
        m.coef = coef(m.fit)
        
        freq.pred <- pbeta(d, N*coef(m.fit)*p.list, N*coef(m.fit)*(1-p.list), lower.tail=FALSE)
        Rsqr <- 1 - (sum((freq.list - freq.pred)^2))/(sum((freq.list - mean(freq.list))^2))
        
        # Get table of model fit stats
        fitstats <- data.frame(m=m.coef, m.low.ci=m.ci[1], m.up.ci=m.ci[2],
                               Rsqr=Rsqr, p.value=m.sum$parameters[4], N=N,
                               Samples=nrow(OTU.table), Richness=length(p.list),
                               Detect=d)
        
        # Get confidence interval for predictions
        freq.pred.ci <- binconf(freq.pred*nrow(OTU.table), nrow(OTU.table), alpha=0.05, method="wilson", return.df=TRUE)
        
        # Get table of predictions
        pred.df <- data.frame(metacomm_RA=p.list, frequency=freq.pred,
                              frequency_lowerCI=freq.pred.ci[,2],
                              frequency_upperCI=freq.pred.ci[,3]) %>%
            unique()
        
        # Get table of observed occupancy and abundance
        obs.df = C.no0 %>%
            rename(metacomm_RA = p, frequency=freq)
        
        head(obs.df)
        
        
        
        p = ggplot(data=obs.df) +
            geom_point(data=obs.df, aes(x=log10(metacomm_RA), y=frequency),
                       alpha=.2, size=2, color="grey30") +
            geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency), color="black") +
            geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_lowerCI), linetype=2, color="black") +
            geom_line(data=pred.df, aes(x=log10(metacomm_RA), y=frequency_upperCI), linetype=2, color="black") +
            # geom_text(data=fitstats, aes(label = paste("R^2 == ", round(Rsqr, 3))),
            #           x=1, y=0.75, size=4, parse=TRUE) +
            # geom_text(data=fitstats, aes(label = paste("italic(m) ==", round(m, 3))),
            #           x=-1, y=0.85, size=4, parse=TRUE) +
            labs(x="Log10 abundance in\nmetacommunity", y="Frequency detected",
                 title = paste(aa[i],paste("R2 = ",round(fitstats$Rsqr, 3)),
                               paste("m = ",round(fitstats$m, 3)))) +
            theme_bw() +
            theme(axis.line = element_line(color="black"),
                  legend.position = "none",
                  axis.title = element_text(size=14),
                  axis.text = element_text(size=12))
        
        p
        
        
        
        plots[[aa[i]]] = p
        
    }
    
    
    # plots$ABCD
    # library(ggpubr)
    # nrow=2,,ncol=4
    p  = ggpubr::ggarrange(plotlist = plots,  common.legend = TRUE, legend="right",ncol = 2,nrow = 4)
    p
    
    return(list(p,plots))
    
}

results = neutralModel(ps = ps,group = "Treat1")

p <- results[[1]]

png("NM.png",width = 5400,height = 8400,res = 600)
p
dev.off()

##计算RCbray
RCbary = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",num = 99,thread = 1){
    ps_sub <- ps
    #----------------整理map文件
    map = as.data.frame(sample_data(ps_sub))
    map$ID = row.names(map)
    sample_data(ps) = map
    #-------------------准备OTU表格
    #-----------------抽平-不设置抽平条数，默认按照最小序列数数目抽平
    set.seed(72)  # setting seed for reproducibility
    psrare = rarefy_even_depth(ps_sub )
    #检查序列数量
    sample_sums(psrare)
    # 标准化数据
    ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
    
    
    
    
    #--------------两个函数
    # 对模拟群落计算距离
    RCbray_null_func <- function(i, freq.abd.df, alpha1, alpha2, N){
        # Get simulated communities and distance
        ## initally select OTUs weighted by their frequency. The number of OTUs selected should equal the richness of the samples.
        simcom1 = data.frame(table(sample(freq.abd.df$OTU, size=alpha1, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
        colnames(simcom1) = c("OTU","simcom1")
        simcom1$OTU = as.character(simcom1$OTU)
        simcom1 = inner_join(simcom1, freq.abd.df, by="OTU")
        simcom2 = data.frame(table(sample(freq.abd.df$OTU, size=alpha2, replace=F, prob=freq.abd.df$freq)), stringsAsFactors = F)
        colnames(simcom2) = c("OTU","simcom2")
        simcom2$OTU = as.character(simcom2$OTU)
        simcom2 = inner_join(simcom2, freq.abd.df, by="OTU")
        
        ## Now recruit OTUs based on their abundance in the metacommunity
        simcom1.abd = data.frame(table(sample(simcom1$OTU, size=N-alpha1, replace=T, prob=simcom1$p)), stringsAsFactors = F)
        colnames(simcom1.abd) = c("OTU","simcom1.abd")
        simcom1.abd$OTU = as.character(simcom1.abd$OTU)
        simcom1 = full_join(simcom1, simcom1.abd, by="OTU") %>%
            mutate(simcom1.abd = ifelse(is.na(simcom1.abd), 1, simcom1.abd)) %>%
            select(OTU, simcom1.abd)
        
        simcom2.abd = data.frame(table(sample(simcom2$OTU, size=N-alpha2, replace=T, prob=simcom2$p)), stringsAsFactors = F)
        colnames(simcom2.abd) = c("OTU","simcom2.abd")
        simcom2.abd$OTU = as.character(simcom2.abd$OTU)
        simcom2 = full_join(simcom2, simcom2.abd, by="OTU") %>%
            mutate(simcom2.abd = ifelse(is.na(simcom2.abd), 1, simcom2.abd)) %>%
            select(OTU, simcom2.abd)
        
        
        simcom = full_join(simcom1, simcom2, by="OTU")
        simcom[is.na(simcom)] = 0
        rownames(simcom) = simcom$OTU
        simcom$OTU = NULL
        
        null.dist = vegdist(t(simcom), method="bray")[1]
        return(null.dist)
    }
    
    # 计算RCbray的主功能
    Calc_RCbray <- function(physeq, reps, nproc){
        # Get OTU table from phyloseq object
        otu.table = otu_table(physeq)
        
        # Get alpha diversity for each sample
        otu.PA.table = otu.table
        otu.PA.table[otu.PA.table > 0] = 1
        alpha.df = data.frame(Sample_ID = colnames(otu.PA.table), OTU.n = colSums(otu.PA.table), stringsAsFactors = F)
        
        # Get beta diversity matrix
        beta.table = as.matrix(vegdist(t(otu.PA.table), method="bray", diag=TRUE, upper=TRUE))
        
        ## Get metacommunity
        # Calculate the number of individuals in the meta community (Average read depth)
        N <- mean(apply(t(otu.table), 1, sum))
        
        # Calculate the average relative abundance of each taxa across communities
        p.m <- apply(t(otu.table), 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
        
        # Calculate the occurrence frequency of each taxa across communities
        otu.table.bi <- 1*(t(otu.table)>0)
        freq <- apply(otu.table.bi, 2, mean)
        freq <- freq[freq != 0]
        
        # Combine
        freq.abd.df = data.frame(p=p, freq=freq) %>%
            tibble::rownames_to_column(var="OTU") %>%
            filter(p != 0, freq != 0) %>%
            arrange(p)
        
        # For each pair of samples run the RCbray analysis
        comps = combn(alpha.df$Sample_ID, m=2, simplify = F)
        RCb.df = data.frame(Site1 = character(), Site2 = character(), RCb = numeric(), stringsAsFactors = F)
        for (j in seq(1, length(comps))){
            sam = comps[[j]]
            alpha1 = alpha.df[alpha.df$Sample_ID == sam[1],]$OTU.n
            alpha2 = alpha.df[alpha.df$Sample_ID == sam[2],]$OTU.n
            # Permute "reps" many times
            rep.list = seq(1, reps)
            null.list = mclapply(rep.list, RCbray_null_func, freq.abd.df=freq.abd.df, alpha1=alpha1, alpha2=alpha2, N=N, mc.cores=nproc)
            
            RCb = (length(null.list[null.list > beta.table[sam[1], sam[2]]]) + (0.5*length(null.list[null.list == beta.table[sam[1], sam[2]]])))/reps
            RCb = (RCb - 0.5)*2
            
            RCb.df = rbind(RCb.df, data.frame(Site1=sam[1], Site2=sam[2], RCb=RCb, stringsAsFactors = F))
        }
        
        RCb.df
        return(RCb.df)
    }
    
    
    # 运行RCbray的计算，这个运算再5个小时左右999重复
    RCb = Calc_RCbray(psrare, num, thread)
    
    head(RCb)
    
    return(list(RCb))
}

result = RCbary(ps = ps ,group  = "Treat",num = 3,thread = 1)
