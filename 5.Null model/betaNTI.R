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

bNTICul = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",num = 99,thread = 1){
    ps = inputMicro(otu,tax,map,tree,ps,group  = group)
    ps
    
    ps_sub <- ps
    
    # tree = phy_tree(ps)
    # tree
    #-------------调整map文件-----------------------------------------------------------------
    #添加一个ID列
    map = as.data.frame(sample_data(ps_sub))
    map$ID = row.names(map)
    sample_data(ps) = map
    
    
    #-----------准备OTU表格---------------------抽平-不设置抽平条数，默认按照最小序列数数目抽平
    set.seed(72)  # setting seed for reproducibility
    psrare = rarefy_even_depth(ps_sub)
    #检查序列数量
    sample_sums(psrare)
    # 标准化数据
    ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))
    
    
    # ## -0--------------对样本序列数量可视化
    # read_count = data.frame("count" = colSums(otu_table(ps))) %>%
    #   rownames_to_column(var="ID") %>%
    #   inner_join(data.frame(sample_data(ps)), by="ID") %>%
    #   arrange(-count) %>%
    #   mutate(ID=factor(ID, levels=ID))
    #
    # # Now plot read count for each sample. The horizontal line represents a 2000 read threshold
    # ggplot(data=read_count, aes(x=ID, y=log10(count), fill=Group)) +
    #   geom_bar(stat="identity") +
    #   labs(x="Sample", y="Log10(Read count)") +
    #   geom_hline(yintercept=log10(10000)) +
    #   theme(text = element_text(size=16),
    #         axis.text.x = element_blank())
    # # Everything seems to be at or above 10000 total reads
    #
    # ps
    
    
    
    
    
    # 计算βMNTD对每个随机零模型群落
    bMNTD_null_func <- function(i, OTU.table, tree){
        tree$tip.label = sample(tree$tip.label)
        bMNTD_s = comdistnt(OTU.table[,colnames(dis)], cophenetic(prune.sample(OTU.table,tree)), abundance.weighted = TRUE)
        A <- attr(bMNTD_s, "Size")
        B <- if (is.null(attr(bMNTD_s, "Labels"))) sequence(A) else attr(bMNTD_s, "Labels")
        if (isTRUE(attr(bMNTD_s, "Diag"))) attr(bMNTD_s, "Diag") <- FALSE
        if (isTRUE(attr(bMNTD_s, "Upper"))) attr(bMNTD_s, "Upper") <- FALSE
        bMNTD_s.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                                Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                                bMNTD = as.vector(bMNTD_s),
                                rep=i)
        return(bMNTD_s.df)
    }
    # 计算βNTI
    Phylo_turnover <- function(physeq, reps, nproc){
        # Extract OTU table
        OTU.table = t(otu_table(physeq))
        # Extract phylogenetic tree
        tree = phy_tree(physeq)
        dis = prune.sample(OTU.table,tree)
        # Get βMNTD between all communities
        bMNTD_o = comdistnt(OTU.table[,colnames(dis)], cophenetic(dis), abundance.weighted = TRUE)
        A <- attr(bMNTD_o, "Size")
        B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
        if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
        if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
        bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                                Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                                bMNTD = as.vector(bMNTD_o))
        
        # Get βMNTD for randomized null communities
        rep.list = seq(1, reps)
        bMNTD_s.df.list = mclapply(rep.list, bMNTD_null_func, OTU.table=OTU.table, tree=tree, mc.cores=nproc)
        
        # Combine all data together and calculate βNTI for each sample pair
        bMNTD_s.df <- do.call("rbind", bMNTD_s.df.list)
        bMNTD_s.means.df = bMNTD_s.df %>%
            group_by(Sample_1, Sample_2) %>%
            dplyr::summarize(mean_bMNTD = mean(bMNTD),
                             sd_bMNTD = sd(bMNTD))
        
        bMNTD_o.df = inner_join(bMNTD_o.df, bMNTD_s.means.df, by=c("Sample_1", "Sample_2")) %>%
            mutate(bNTI = (bMNTD - mean_bMNTD)/sd_bMNTD)
        
        return(bMNTD_o.df)
    }
    
    
    #========这里一把单核就真实数据而言需要超过10个小时，跑999次，所以需要多核
    # 计算bnti，这里可以设置线程数量，是第三个参数，我们在linux下面可以设置，30个线程
    # 第二个参数设置迭代数量，这里文献一般999嘛。
    bNTI = Phylo_turnover(psrare, num, thread)
    
    return(list(bNTI))
}

set.seed(2516)

result = bNTICul(ps = ps ,group  = "Treat1",num = 1,thread = 1)
#
bNTI = result[[1]]
head(bNTI)
#
# path = "./Result/bNTI"
# dir.create(path)
# filename = paste(path,"/bNTI.csv",sep = "")
# write.csv(bNTI, filename)


