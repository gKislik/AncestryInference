library(dplyr)
library(tidyverse)
library(matrixStats)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(ggrepel)
library(data.table)
library(ggpubr)
library(vctrs)
library(plyr)
library(ggrepel)
library(Metrics)
library(irr)
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\bootstrapping')

all_iters <- data.frame(sample=as.character(),Breed=as.character(),Proportion=as.numeric(),iter=as.numeric())
for(i in seq(1,100)){
  Qhat <- read.table(paste0('Qhat_rep',i,'.txt'))
  freq <- read.table(paste0('freq_rep',i,'.txt'),header=T)
  rownames(Qhat) <- freq$CLST
  fam <- read.table(paste0('tempinp_rep',i,".fam"))
  clust <- read.table(paste0('trimclust_rep',i,".txt"))
  twocolclust <- clust[,c(1,3)]
  colnames(twocolclust) <- c('sample','breed')
  colnames(Qhat) <- fam$V1
  tr_Qhat <- t(Qhat)
  tr_Qhat <- as.data.frame(t(Qhat))
  tr_Qhat$sample <- rownames(tr_Qhat)
  long_tr_Qhat <- pivot_longer(data.frame(tr_Qhat), names_to='Breed',values_to='Proportion',cols=!sample)
  long_tr_Qhat$iter <- i
  print(head(long_tr_Qhat))
  all_iters <- rbind(all_iters, long_tr_Qhat)
  
}

setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')
clust <- read.table("modsixclust.txt")
twocolclust <- clust[,c(1,3)]
colnames(twocolclust) <- c('sample','Breed')

icc_res <- data.frame(sample=as.character(), icc=as.numeric())
for(samp in unique(all_iters$sample)){
  sdf <- all_iters[all_iters$sample==samp,]
  sdf_w <- pivot_wider(sdf ,values_from = Proportion,names_from = iter)
  sdf_w <- sdf_w[,-c(1,2)]
  intracc <- icc(sdf_w)
  print(intracc$value)
  re <- data.frame(sample=samp, icc=intracc$value)
  icc_res <- rbind(icc_res, re)
  
}

icc_res <- merge(twocolclust, icc_res, by='sample')
write.csv(icc_res,'icc_res.csv')
