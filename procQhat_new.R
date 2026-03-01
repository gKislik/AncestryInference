### new Qhat Processing

library(tidyr)
library(dplyr)
library(ggplot2)
library(ape)
library(ggtree)
library(uwot)
library(RcppHNSW)
library(Metrics)
library(ggrepel)
library(irr)

#beginning of commands to actually generate figures 
#after filtering for breed-specific SNPs, FST and removing uninformative samples
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF')
Qhat <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350Qhat.txt")
freq <- read.table("hf.txt",header=T)
rownames(Qhat) <- freq$CLST
fam <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.fam")
clust <- read.table("EXP-G-corrected-clust.txt")
twocolclust <- clust[,c(1,3)]
colnames(twocolclust) <- c('sample','breed')

colnames(Qhat) <- fam$V1

firstbreed <- !duplicated(clust$V3)
newlabs <- ifelse(firstbreed==T,clust$V3,"")

truth <- data.frame(table(clust$V2,clust$V3))
tr_Qhat <- t(Qhat)
tr_Qhat <- as.data.frame(t(Qhat))
tr_Qhat$sample <- rownames(tr_Qhat)
long_tr_Qhat <- pivot_longer(data.frame(tr_Qhat), names_to='Breed',values_to='Proportion',cols=!sample)
l_tr_tb_Qhat <- merge(long_tr_Qhat, twocolclust, by = 'sample')
colnames(l_tr_tb_Qhat) <- c('sample','Breed','Proportion','TrueBreed')
l_tr_tb_Qhat$Breed[l_tr_tb_Qhat$Breed == l_tr_tb_Qhat$TrueBreed] <- 'Correct'
l_correct <- l_tr_tb_Qhat[l_tr_tb_Qhat$Breed == 'Correct',]
l_correct <- l_correct[order(l_correct$TrueBreed),]
correct_plot <- ggplot(l_correct, aes(x=sample, y = Proportion))+
  geom_bar(stat='identity') + facet_wrap(vars(TrueBreed), scales='free_x')+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

correct_plot


#### add SD
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\bootstrapping')

all_iters <- data.frame(sample=as.character(),Breed=as.character(),Proportion=as.numeric(),iter=as.numeric())
for(i in seq(1,100)){
  Qhat <- read.table(paste0('Qhat_rep',i,'.txt'))
  freq <- read.table(paste0('freq_rep',i,'.txt'),header=T)
  rownames(Qhat) <- freq$CLST
  fam <- read.table(paste0('tempinp_rep',i,".fam"))
  clusta <- read.table(paste0('trimclust_rep',i,".txt"))
  twocolclusta <- clusta[,c(1,3)]
  colnames(twocolclusta) <- c('sample','breed')
  colnames(Qhat) <- fam$V1
  tr_Qhat <- t(Qhat)
  tr_Qhat <- as.data.frame(t(Qhat))
  tr_Qhat$sample <- rownames(tr_Qhat)
  long_tr_Qhata <- pivot_longer(data.frame(tr_Qhat), names_to='Breed',values_to='Proportion',cols=!sample)
  long_tr_Qhata$iter <- i
  print(head(long_tr_Qhata))
  all_iters <- rbind(all_iters, long_tr_Qhata)
  
}

colnames(twocolclust) <- c('sample','Breed')

icc_res <- data.frame(sample=as.character(), icc=as.numeric())
sd_res <- data.frame(matrix(nrow=65,ncol = 0))
for(samp in unique(all_iters$sample)){
  sdf <- all_iters[all_iters$sample==samp,]
  sdf_w <- pivot_wider(sdf ,values_from = Proportion,names_from = iter)
  sdf_w <- sdf_w[,-c(1,2)]
  print(class(sdf_w))
  stddev <- rowSds(as.matrix(sdf_w))
  print(stddev)
  sd_res[,samp] <- stddev
  intracc <- icc(sdf_w)
  print(intracc$value)
  re <- data.frame(sample=samp, icc=intracc$value)
  icc_res <- rbind(icc_res, re)
  
}

icc_res <- merge(twocolclust, icc_res, by='sample')
write.csv(icc_res,'icc_res.csv')

rownames(sd_res) <- freq$CLST
sd_res$Breed <- freq$CLST
sd_l <- pivot_longer(sd_res,cols=!Breed,names_to = 'sample',values_to = 'sd')
sd_l_b <- merge(sd_l, twocolclust, by=c('sample','Breed'))
colnames(sd_l_b) <- c('sample','TrueBreed','sd')
l_sd_correct <- merge(l_correct, sd_l_b, by=c('sample','TrueBreed'))

correct_plot_sd <- ggplot(l_sd_correct, aes(x=sample, y = Proportion, ymin=pmax(0,Proportion-sd), ymax=pmin(Proportion+sd,1)))+
  geom_bar(stat='identity') + facet_wrap(vars(TrueBreed), scales='free_x', ncol = 13)+geom_errorbar(width=0.2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylab('Proportion Correct')

correct_plot_sd

ggsave('QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350Boot.png',correct_plot_sd,width = 16,height=8,units = 'in')


### BSPEC SNP FST Distributions

setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\pretwo')
distfiles <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\pretwo", pattern="*-pre2.txt*", full.names = FALSE, no.. = TRUE)

distances <- data.frame(matrix(NA, nrow = 10000, ncol = 0))
for(file in distfiles){
  temp <- read.table(file, header=F)[,5]
  length(temp) <- 10000
  distances[,file] <- as.numeric(temp)
}
distances$temp <- 1
dlm <- cbind(str_split_fixed(colnames(distances),pattern='-',2)[,1],colMeans(distances))
colnames(dlm) <- c('Breed','m')
dl <- pivot_longer(distances, names_to = 'Breed', values_to = 'FST',cols=!temp)
dl$Breed <- str_split_fixed(dl$Breed,pattern='-',2)[,1]
dl <- merge(dl, dlm)
dl$m <- as.numeric(dl$m)
dvio <- ggplot(dl,aes(x=reorder(Breed, m), y=FST,fill=Breed))+geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
        #legend.title = element_text(size = 10), 
        #legend.text = element_text(size = 10), 
        legend.position = "none",
        axis.text.y = element_text(size=15,colour = 'black')) #+
  #guides(color = guide_legend(override.aes = list(size = 0.25))) +
  #guides(shape = guide_legend(override.aes = list(size = 0.25)))
  
ggsave('modFSTdistributions.png',dvio,width=16,height=12,units='in')

