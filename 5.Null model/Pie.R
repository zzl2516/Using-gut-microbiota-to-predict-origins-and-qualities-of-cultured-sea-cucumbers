data <- read.table("rate.txt",header = TRUE,sep = "\t")

## XP
pct <- round(data$XP/sum(data$XP)*100,2)  
lab <- paste(pct,"%",sep="")
XP <- data.frame(data$Site,data$XP,lab)

## SZ
pct <- round(data$SZ/sum(data$SZ)*100,2)  
lab <- paste(pct,"%",sep="")
SZ <- data.frame(data$Site,data$SZ,lab)

## YT
pct <- round(data$YT/sum(data$YT)*100,2)  
lab <- paste(pct,"%",sep="")
YT <- data.frame(data$Site,data$YT,lab)

## QD
pct <- round(data$QD/sum(data$QD)*100,2)  
lab <- paste(pct,"%",sep="")
QD <- data.frame(data$Site,data$QD,lab)

## WFD
pct <- round(data$WFD/sum(data$WFD)*100,2)  
lab <- paste(pct,"%",sep="")
WFD <- data.frame(data$Site,data$WFD,lab)

## BCD
pct <- round(data$BCD/sum(data$BCD)*100,2)  
lab <- paste(pct,"%",sep="")
BCD <- data.frame(data$Site,data$BCD,lab)

## DY
pct <- round(data$DY/sum(data$DY)*100,2)  
lab <- paste(pct,"%",sep="")
DY <- data.frame(data$Site,data$DY,lab)

## JZ
pct <- round(data$JZ/sum(data$JZ)*100,2)  
lab <- paste(pct,"%",sep="")
JZ <- data.frame(data$Site,data$JZ,lab)

cbbPalette <- c("#B2182B","#E69F00","#56B4E9","#009E73","#F0E442")

par(mfrow = c(2,4))
pie(XP$data.XP,labels = XP$lab,col = cbbPalette,main = "XP (n = 33)")
pie(SZ$data.SZ,labels = SZ$lab,col = cbbPalette,main = "QHD (n = 18)")
pie(YT$data.YT,labels = YT$lab,col = cbbPalette,main = "YT (n = 31)")
pie(QD$data.QD,labels = QD$lab,col = cbbPalette,main = "QD (n = 30)")
pie(WFD$data.WFD,labels = WFD$lab,col = cbbPalette,main = "DLW (n = 37)")
pie(BCD$data.BCD,labels = BCD$lab,col = cbbPalette,main = "DLE (n = 32)")
pie(DY$data.DY,labels = DY$lab,col = cbbPalette,main = "DY (n = 30)")
pie(JZ$data.JZ,labels = JZ$lab,col = cbbPalette,main = "JZ (n = 78)")

## 图例
pie(XP$data.XP,labels = XP$lab,col = cbbPalette,main = "XP (n = 33)")
legend("right",legend = XP$data.Site,bty = "n",inset = c(-0.4,0),xpd = TRUE,
       fill = cbbPalette)
