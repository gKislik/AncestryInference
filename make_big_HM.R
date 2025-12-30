library(ggplot2)
#library(biganalytics)
library(tidyr)
library(pheatmap)

set.seed(7)

freqs <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.frq.strat",header = T)
freqs_1 <- freqs[,c(2,3,6)]
freqs_1w <- data.frame(pivot_wider(freqs_1,names_from=CLST,values_from=MAF))
rownames(freqs_1w) <- freqs_1w[,1]
freqs_1w <- freqs_1w[,-c(1)]

pheatmap(scale(freqs_1w), kmeans=65, filename="QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350-HM.png",fontsize=8)

pheatmap(freqs_1w, kmeans=65, filename="QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350-HM-US.png",fontsize=8)

###### not in use ######

#pheatmap(scale(freqs_1w), kmeans=50,treeheight_row = 0, treeheight_col = 0,filename="filt-ref-SRR8614076-0p03fst-HM.png",fontsize=8)

#kres <- bigkmeans(as.matrix(freqs_1w),50,1000)
#cent <- data.frame(kres$centers)
#colnames(cent) <- colnames(freqs_1w)
#cent$cluster <- as.numeric(rownames(cent))
#cent_l <- pivot_longer(cent,cols=!cluster,names_to="BREED",values_to="MAF")

#p<- ggplot(cent_l,aes(x=cluster,y=BREED,fill=MAF))+geom_tile()+scale_x_continuous(expand=c(0,0))

#ggsave("filt-ref-SRR8614076-0p03fst-HM.png",width=12, height=9, units='in')

#get cluster sizes
#clust_sizes <- cbind(seq(1,50),kres$size)
#colnames(clust_sizes) <- c("cluster", "size")
#write.table(clust_sizes, "filt-ref-SRR8614076-0p03fst-HM-seed7-clust-sizes.txt",row.names=F)
