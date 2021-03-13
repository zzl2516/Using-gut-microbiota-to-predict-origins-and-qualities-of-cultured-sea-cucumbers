phylum <- read.table("phylum.txt",header = TRUE, sep = "\t",row.names = 1)
phylum <- as.matrix(phylum)
library(RColorBrewer)

png(filename = "phylum.png",width = 9600,height = 2400,res = 600)
par(mar = c(8,1,3,3))
par(xpd = TRUE)
barplot(phylum,horiz = TRUE,axes = F,names.arg = "",ylim = c(0,100),
        col = brewer.pal(6,"Set1"),width = 50)
text(x = 35.13,y = 33,labels = "Proteobacteria",font = 2,cex = 3.5)
segments(x0 = 0,y0 = -0.5,x1 = 0,y1 = 70,col = "black",lwd = 4)
axis(side = 1,at = c(0,20,40,60,80,100),labels = c(0,20,40,60,80,100),
     font.axis = 2,lwd = 4,lwd.ticks = 3,cex.axis = 2.5,padj = 0.5)
mtext("The relative abundance (%) of dominant phyla in gut microbiota of sea cucumbers",
      side = 1,font = 2,line = 4.2,cex = 2.3,at = 50)
legend(1,100,legend = rownames(phylum)[2:6],fill = brewer.pal(6,"Set1")[2:6],
       bty = "n",ncol = 5,cex = 1.8,text.font = 2)
dev.off()

