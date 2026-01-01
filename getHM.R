args<-commandArgs(TRUE)

library(ggplot2)
library(biganalytics)
library(tidyr)
library(pheatmap)

freqs <- read.table(args[1],header = T)
freqs_1 <- freqs[,c(2,3,6)]
freqs_1w <- data.frame(pivot_wider(freqs_1,names_from=CLST,values_from=MAF))
rownames(freqs_1w) <- freqs_1w[,1]
freqs_1w <- freqs_1w[,-c(1)]

pheatmap(scale(freqs_1w), kmeans=7,treeheight_row = 0, treeheight_col = 0,filename=paste0(args[2],"_HM.png"),fontsize=8)

#kres <- bigkmeans(as.matrix(freqs_1w),7,1000)
#cent <- data.frame(kres$centers)
#colnames(cent) <- colnames(freqs_1w)
#cent$cluster <- as.numeric(rownames(cent))
#cent_l <- pivot_longer(cent,cols=!cluster,names_to="BREED",values_to="MAF")

#p<- ggplot(cent_l,aes(x=cluster,y=BREED,fill=MAF))+geom_tile() 

#ggsave(paste0(args[2],"_HM.png"))
