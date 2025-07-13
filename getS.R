args<-commandArgs(TRUE)

library(ggplot2)

SNPs <- read.table(args[1],header=F)

ggplot(SNPs,aes(x=V5))+geom_histogram()

ggsave(paste0(args[2],"_Dist.png"))
