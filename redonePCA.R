### plot pca
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
library(uwot)
setwd("C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF")
finclusters <- read.table("EXP-G-corrected-clust.txt")
eigenvec <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.eigenvec")
twocolclust <- finclusters[,c(1,3)]
colnames(twocolclust) <- c('sample','breed')

dist_from_center <- function(umap_rel_mat){
  distances <- c()
  for(sample in unique(umap_rel_mat$sample)){
    b <- umap_rel_mat$breed[umap_rel_mat$sample == sample]
    center_PC1 <- mean(umap_rel_mat$Component1[umap_rel_mat$breed == b])
    center_PC2 <- mean(umap_rel_mat$Component2[umap_rel_mat$breed == b])
    s_PC1 <- umap_rel_mat$Component1[umap_rel_mat$sample == sample]
    s_PC2 <- umap_rel_mat$Component2[umap_rel_mat$sample == sample]
    distance <- sqrt((s_PC2-center_PC2)**2 + (s_PC1-center_PC1)**2)
    distances <- c(distances,distance)
  }
  umap_output <- umap_rel_mat
  umap_output$distance_from_center <- distances
  return(umap_output)
}

PCA_df <- merge(finclusters,eigenvec,by='V1')

PCA_plot <- ggplot(PCA_df , aes(x=V3.y,y=V4,color=V3.x,label=V3.x))+geom_point()+
  theme_bw()+labs(x='PC1',y='PC2')+ guides(fill=guide_legend(title="Breed"))+
  geom_text_repel(data = distinct(PCA_df,V3.x,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 50,max.iter = 2500)
#PCA_plot
PCA_plot_nt <- ggplot(PCA_df , aes(x=V3.y,y=V4,color=V3.x,label=V3.x))+geom_point()+
  theme_bw() + labs(x='PC1',y='PC2') +
  guides(fill=guide_legend(title="Breed"))#+geom_text_repel(data = distinct(PCA_df,V3.x,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 50,max.iter = 2500)
#PCA_plot_nt
ggsave("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.rel_PCA_CORRECTED.png",height = 9,width = 12,units = 'in',plot = PCA_plot)
ggsave("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.rel_PCA_CORRECTEDnotext.png",height = 9,width = 12,units = 'in',plot = PCA_plot_nt)


# UMAP for validation, mdist generated with below command:
# plink -bfile bcftools-no0-ref-redone-allvar-QUAL20-GENO0.1-MAF0.01-LD5_50_0.5 -chr-set 38 -allow-extra-chr -make-rel -out bcftools-no0-ref-redone-allvar-QUAL20-GENO0.1-MAF0.01-LD5_50_0.5

set.seed(100)
rel <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.rel") 
rel_id <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.rel.id")

rel_umap <- as.data.frame(umap2(rel, seed = 100, ret_model = FALSE, nn_method = 'annoy' ))
colnames(rel_umap) <- c("Component1","Component2")

rel_umap$sample <- rel_id$V1
rel_umap <- merge(rel_umap, twocolclust,by='sample')

UMAP_plot <- ggplot(rel_umap, aes(x=Component1, y=Component2, color=breed,label=breed))+
  geom_point() +theme_bw()+ theme(legend.position = 'none') +#xlim(-50,50)+ylim(-45,45)+
  geom_text_repel(show.legend = F, data = distinct(rel_umap,breed,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 100,max.iter = 2500)
UMAP_plot
ggsave("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350-UMAP-CORRECTEDtext.png",height = 9,width = 12,units = 'in',plot = UMAP_plot)

umap_distances <- dist_from_center(rel_umap)
write.csv(umap_distances,"QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350-UMAP-CORRECTEDdistances.csv",row.names = F,col.names = T,quote = F )
plot(density(umap_distances$distance_from_center))

###phylogenetic tree

dist <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.dist") 
IDs_dist <- read.table("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350.dist.id")

rownames(dist) <- IDs_dist[,1]
colnames(dist) <- IDs_dist[,1]

tree <- nj(as.matrix(dist))
tree <- root(tree, "SRR30817380")
twocolclust <- clust[,c(1,3)]
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
ggsave("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350-FST0.350-phylocircle.png",height = 10,width = 10,units = 'in',plot = updp_treeplot)
updp_breed_treeplot
ggsave("QUAL20-0-G-MAF0.01-GENO0.2-LD250_50_0.8-10kSNPs_FST0.350-FST0.350-phylocircle-breed.png", height = 10, width = 10, units= 'in',plot=updp_breed_treeplot)

