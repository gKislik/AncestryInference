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
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')
Qhat <- read.table("filt-bcftools-ref-redone.imputed.filteredQhat.txt")
freq <- read.table("hf.txt",header=T)
rownames(Qhat) <- freq$CLST
fam <- read.table("filt-bcftools-ref-redone.imputed.filtered.fam")
clust <- read.table("used-clust.txt")
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
long_tr_Qhat$group <- "Estimate"
colnames(truth) <- c('sample','Breed','Proportion')
truth$group <- "Truth"
fulltotal <- rbind(truth, long_tr_Qhat)

#plotting

ftotal_plot <- ggplot(fulltotal, aes(x=sample,y=Proportion,fill=Breed))+
  geom_bar(position='stack',stat='identity')+theme_bw()+facet_wrap(~group,nrow=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))+
  scale_x_discrete(limits=clust$V1)
ftotal_plot_breed <- ggplot(fulltotal, aes(x=sample,y=Proportion,fill=Breed,group=Breed))+
  geom_bar(position='stack',stat='identity')+theme_bw()+facet_wrap(~group,nrow=2)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=12))+
  scale_x_discrete(limits=clust$V1,label=newlabs)
ftotal_plot_breed
ggsave("filt-bcftools-ref-redone.imputed.filteredQhatvTruth-Breed.png",width=24,height=6,units='in')
ftotal_plot
ggsave("filt-bcftools-ref-redone.imputed.filteredQhatvTruth.png",width=24,height=6,units='in')


### bootstrapping + differences
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
clust <- read.table("used-clust.txt")
twocolclust <- clust[,c(1,3)]
colnames(twocolclust) <- c('sample','Breed')

icc_res <- data.frame(sample=as.character(), icc=as.numeric())
sd_res <- data.frame(matrix(nrow=48,ncol = 0))
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

#get RMSE and Pearson R and Differences
fulltotal_w <- pivot_wider(fulltotal, values_from = 'Proportion', names_from = 'group')
stats_df <- data.frame(RMSE=as.numeric(),PearsonR=as.numeric(), Correct=as.numeric(),sample=as.character())

for (samp in unique(fulltotal$sample)){
  RMSE <- rmse(fulltotal_w$Truth[fulltotal_w$sample==samp],fulltotal_w$Estimate[fulltotal_w$sample==samp])
  Corr <- cor(fulltotal_w$Truth[fulltotal_w$sample==samp],fulltotal_w$Estimate[fulltotal_w$sample==samp])
  bdiff <- (fulltotal_w[fulltotal_w$sample==samp&fulltotal_w$Truth==1,'Estimate'])
  stats_vec <- data.frame(RMSE,Corr,bdiff,samp)
  colnames(stats_vec) <- c('RMSE','PearsonR','Correct','sample')
  stats_df <- rbind(stats_df, stats_vec)
}


stats_df <- merge(stats_df, sd_l_b,by='sample')
write.table(stats_df,"filt-bcftools-ref-redone.imputed.filteredQhat-ReferenceRecovery.txt",row.names = F)

differences <- ggplot(stats_df, aes(x=reorder(sample,desc(Correct)),y=Correct,fill=Breed))+geom_bar(stat='identity')+
  theme_bw()+geom_errorbar(aes(ymin=Correct-sd, ymax=pmin(1,Correct+sd)), width=.1,
                           position=position_dodge(.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=8),
        panel.background = element_rect(fill = "white")) +labs(x='sample')+
  scale_x_discrete(limits=clust$V1,label=newlabs)
ggsave('Differences.png',differences,width = 16,height=8,units = 'in')


# 1-IBS
dist <- read.table("filt-bcftools-ref-redone.imputed.filtered.mdist") 
IDs_dist <- read.table("filt-bcftools-ref-redone.imputed.filtered.mdist.id")
# IBS
#dist <- read.table("better-filt-ref-SRR8614076-0p03fst-dist.mibs")
#IDs_dist <- read.table("better-filt-ref-SRR8614076-0p03fst-dist.mibs.id")
# Allele counts for dist matrix
#dist <- read.table("better-filt-ref-SRR8614076-0p03fst-dist.dist")
#IDs_dist <- read.table("better-filt-ref-SRR8614076-0p03fst-dist.dist.id")
rownames(dist) <- IDs_dist[,1]
colnames(dist) <- IDs_dist[,1]

tree <- nj(as.matrix(dist))
tree <- root(tree, "SRR30817380")
twocolclust <- clust[,c(1,3)] #if using the 1-IBS, IBS, or AC methods
tree$edge.length <- tree$edge.length + 10 #to correct for weird negative branch/edge lengths, doesn't change overall topology
p <- ggtree(tree, layout = 'circular', branch.length = 'none') + theme_tree2() #+geom_tiplab(label = names, size = 1)
p
colnames(twocolclust) <- c('sample','breed')
updp_treeplot <- p %<+% twocolclust + geom_tiplab(aes(color=factor(breed)),geom='text',size=2)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
updp_breed_treeplot <- p %<+% twocolclust + geom_tiplab(aes(label=breed,color=factor(breed)),geom='text',size=2)+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.position = "none")
updp_treeplot
ggsave("filt-bcftools-ref-redone.imputed.filtered-phylocircle-1mIBS.png",height = 10,width = 10,units = 'in',plot = updp_treeplot)
updp_breed_treeplot
ggsave("filt-bcftools-ref-redone.imputed.filtered-phylocircle-breed-1mIBS.png", height = 10, width = 10, units= 'in',plot=updp_breed_treeplot)


# Using Plink2 Hudson FST
#dist_l <- read.table("better-filt-ref-SRR8614076-0p03fst-plink2-fst.fst.summary",header=F)
#colnames(dist_l) <- c('Pop1','Pop2','FST')
#dist_w <- data.frame(pivot_wider(dist_l, values_from = 'FST',names_from = 'Pop1'))
#rownames(dist_w) <- dist_w[,1]
#dist <- as.matrix(dist_w[,2:ncol(dist_w)])
#tdist <- t(dist)[,-nrow(dist)]
#dist[upper.tri(dist)] <- tdist[upper.tri(tdist, diag = F)]
#tree <- nj(as.matrix(1-dist))
#p <- ggtree(tree,layout='circular',branch.length='none') + theme_tree2() #+geom_tiplab(label = names, size = 1)
#p
#updp_treeplot <- p + geom_tiplab(geom='text',size=2, aes(color=label))+
#  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())
#updp_treeplot
#ggsave("better-filt-ref-SRR8614076-0p03fst-phylocircle-1mFSTsummary.png",height = 9,width = 12,units = 'in',plot = updp_treeplot)

#UMAP of plink standardized GRM (make-rel default output)
#trying to mimic how plink does their PCA, just with UMAP as the final step instead (using rel matrix)
#switched to 1-IBS similar to tree

#rel <- read.table("better-filt-ref-SRR8614076-0p03fst-rel.rel")
#rel_id <- read.table("better-filt-ref-SRR8614076-0p03fst-rel.rel.id")

set.seed(100)
rel <- read.table("filt-bcftools-ref-redone.imputed.filtered.mdist") 
rel_id <- read.table("filt-bcftools-ref-redone.imputed.filtered.mdist.id")

rel_umap <- as.data.frame(umap2(rel))
colnames(rel_umap) <- c("Component1","Component2")

rel_umap$sample <- rel_id$V1
rel_umap <- merge(rel_umap, twocolclust,by='sample')

UMAP_plot <- ggplot(rel_umap, aes(x=Component1, y=Component2, color=breed,label=breed))+
  geom_point() +theme_bw()+ theme(legend.position = 'none') +#xlim(-50,50)+ylim(-45,45)+
  geom_text_repel(show.legend = F, data = distinct(rel_umap,breed,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 100,max.iter = 1000)
UMAP_plot
ggsave("filt-bcftools-ref-redone.imputed.filtered-UMAP-text.png",height = 9,width = 12,units = 'in',plot = UMAP_plot)


#distribution violin plots
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')
distfiles <- list.files(path = "setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')", pattern="*-pre2.txt*", full.names = FALSE, no.. = TRUE)

distances <- data.frame(matrix(NA, nrow = 2500, ncol = 0))
for(file in distfiles){
  temp <- read.table(file, header=F)[,5]
  length(temp) <- 2500
  distances[,file] <- as.numeric(temp)
}
distances$temp <- 1
dl <- pivot_longer(distances, names_to = 'Breed', values_to = 'FST',cols=!temp)
dl$Breed <- str_split_fixed(dl$Breed,pattern='-',2)[,1]
dvio <- ggplot(dl,aes(x=Breed, y=FST,fill=Breed))+geom_boxplot()+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10),
        legend.title = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        axis.text.y = element_text(size=15,colour = 'black')) +
  guides(color = guide_legend(override.aes = list(size = 0.25))) +
  guides(shape = guide_legend(override.aes = list(size = 0.25)))
ggsave('modFSTdistributions.png',dvio,width=16,height=12,units='in')


### make diff & dist scatter plot 
avgFST <- aggregate(FST ~ Breed, data=dl, FUN=mean)
avgDiff <- aggregate(. ~ breed, data=stats_df[,-which(colnames(stats_df)=='sample')], FUN=mean)

colnames(avgDiff) <- c('Breed','RMSE','PearsonR','Difference')

FSTandDIFF <- merge(avgDiff,avgFST,by='Breed')
#FSTandDIFF <- merge(FSTandDIFF,icc_res[,-which(colnames(icc_res)=='sample')],by='Breed')
FSTandDIFF_long <- pivot_longer(FSTandDIFF,cols=!c(Breed,FST),names_to = 'method',values_to = 'measure')
FDscatter <- ggplot(FSTandDIFF_long,aes(x=FST,y=measure, color=Breed))+geom_point()+
  geom_smooth(method = "lm", se = FALSE, col = "blue")+theme_bw()+facet_wrap(vars(method))
ggsave('FSTandDiff.png',FDscatter,width=12,height=8,units='in')
