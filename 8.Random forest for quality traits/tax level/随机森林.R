library(randomForest)

data <- read.table("order.txt",header = TRUE,row.names = 1,sep = "\t")
Saponin.rf <- randomForest(CS~., data = data,importance = TRUE, 
                          proximity = TRUE,ntree = 1000)

print(Saponin.rf)

png(filename = "CS.top15.png",width = 5400,height = 3600,res = 600)
varImpPlot(Saponin.rf,n.var = min(15,nrow(Saponin.rf$importance))) 
dev.off()


library(ggplot2)

lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*", "~~italic(r)^2~"="~r2,
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

plotdata <- data.frame(x = data$CS,y = Saponin.rf$predicted)

p <- ggplot(plotdata,aes(x,y)) +
    geom_point(shape = 21,color = "black",size = 6) + 
    geom_smooth(method = 'lm', formula = y ~ x) +
    geom_text(x= min(plotdata$x), y = max(plotdata$y)*1.1,
              label = lm_eqn(plotdata),parse = TRUE,size = 5,hjust = 0) +
    theme_bw() + 
    labs(x = "Actural values",y = "Predicted values",
         title = "CS (mg/kg)") +
    theme(axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(face = "bold",size = 18))

png(filename = "CS.png",width = 2800,height = 2600,res = 600)
p
dev.off()
