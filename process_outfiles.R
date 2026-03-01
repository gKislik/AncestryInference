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
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF')
#eigenvec <- read.table("mod-filt-ref-SRR8614076-fst0p03-pca.eigenvec")
files <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\synth", pattern="*out.txt*", full.names = FALSE, no.. = TRUE)
highCovfiles <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\out_highCoverage", pattern="*out.txt*", full.names = FALSE, no.. = TRUE)
#lowCovfiles <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\out_lowCoverage", pattern="*out.txt*", full.names = FALSE, no.. = TRUE)
downSamplingfiles <- list.files(path="C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\downsampled",pattern="*out.txt",full.names=FALSE,no.. = TRUE)
freq <- read.table('hf.txt',header = T)
clust <- freq$CLST
finclusters <- read.table("EXP-G-corrected-clust.txt")
#truth <- read.csv("C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\new_truth_matrix.csv",header=F)

makedf <- function(file_list,stats_name, delim){
  est_df <- data.frame(Breed=as.character(),Estimate=as.numeric(),Truth=as.numeric(),samp=as.character())
  colnames(est_df) <- 'Breed'
  stats_df <- data.frame(RMSE=as.numeric(),PearsonR=as.numeric(), sample=as.character())
  for(file in file_list){
    x <- read.table(file,header = T)
    x$Breed <- rownames(x)
    #print(head(x))
    colnames(x) <- c('Estimate','Breed')
    name <- strsplit(file,delim)[[1]]
    print(name[[1]])
    #print(truth[truth$V1==name,])
    props <- unlist(str_split(name[[1]],"_"))
    truthvec <- as.data.frame(cbind(clust,rep(0,65)))
    truthvec[,2] <- as.numeric(truthvec[,2])
    colnames(truthvec) <- c('Breed','Truth')
    print(props)
    #props <- props[nzchar(props)]
    for (i in seq(1,length(props)/2)){
      b <- props[i]
      print(b)
      per <- as.numeric(props[i+length(props)/2])
      print(per)
      #have to add in case the same breed was selected multiple times
      truthvec[substr(truthvec[,1],1,10) %like% substr(b,1,10),2] <- truthvec[substr(truthvec[,1],1,10) %like% substr(b,1,10),2] + per
    }
    x <- merge(x,truthvec,by='Breed') #as.vector(truth[truth$V1==name,-c(1,2)])
    x$samp <- name[[1]]  
    assign(name[[1]],x,.GlobalEnv)
    RMSE <- rmse(x$Truth, x$Estimate)
    Corr <- cor(x$Truth, x$Estimate)
    stats_vec <- data.frame(RMSE,Corr,name[[1]])
    colnames(stats_vec) <- c('RMSE','PearsonR','sample')
    stats_df <- rbind(stats_df, stats_vec)
    
    est_df <- rbind(est_df, x)
  }
  assign(stats_name,stats_df,.GlobalEnv)
  write.table(stats_df,paste0(stats_name,".txt"),row.names = F)
  long_est <- gather(est_df,key = 'group',value='proportion',-c(Breed,samp))
  return(long_est)
}


#est <- ggplot(est_df,aes(x=samp,y=Estimate,fill=Breed))+geom_bar(position='stack',stat='identity')+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))                                                                          
#truth <- ggplot(est_df,aes(x=samp,y=Truth,fill=Breed))+geom_bar(position='stack',stat='identity')+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))         
#both <- ggarrange(est,truth,nrow = 2,common.legend = T,legend='right')
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\synth')
synth_df <- makedf(files,'synth_stats', '_simgenotype')
synthplot <- ggplot(synth_df, aes(x=group,y=proportion,fill=Breed))+geom_bar(position='stack',stat='identity')+theme_bw()+
  facet_wrap(~factor(samp, levels = synth_stats$sample[order(synth_stats$RMSE)]), nrow = 10)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none" )+xlab('sample')
synthplot

setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\')
ggsave("synthplot.png",height = 9,width = 12,units = 'in',plot = synthplot)

setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\out_highCoverage')
high_df <- makedf(highCovfiles,'highCov_stats')
high_df$Coverage <- 'high'

#setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\out_lowCoverage')
#low_df <- makedf(lowCovfiles,'lowCov_stats')
#low_df$Coverage <- 'low'

#from the 96 samples
#coverage_df <- rbind(high_df,low_df)
coverage_df <- high_df
##remove rows where breed percentages are <3% to make sure it's easier to read
coverage_df$Breed[coverage_df$proportion < 0.03] <- '<3% breeds'
covplot <- ggplot(coverage_df, aes(x=factor(samp,level=highCov_stats$sample[order(highCov_stats$RMSE)]),y=proportion,fill=Breed))+geom_bar(position='stack',stat='identity')+
  theme_bw()+facet_wrap(~group+Coverage,nrow=2)+xlab('sample')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11))
covplot

setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')
#
ggsave("coverageplot.png",height = 9,width = 12,units = 'in',plot = covplot)

#for downsampling experiment
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF\\downsampled')
downsampled_df <- makedf(downSamplingfiles,'downSampling_stats', "_out.txt")
downsampled_df$props <- mapvalues(downsampled_df$samp,from = 
                                   c("MiniaturePoodle_Poodle_LabradorRetriever_CockerS_.422_.339_.091_.148",
                                     "MiniaturePoodle_Poodle_LabradorRetriever_CockerSp_.422_.339_.091_.148",
                                     "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpa_.422_.339_.091_.148",
                                     "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpan_.422_.339_.091_.148",
                                     "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpani_.422_.339_.091_.148",
                                     "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpanie_.422_.339_.091_.148",
                                     "Pomeranian_Chihuahua_IcelandicSh_.430_.400_.170",
                                     "Pomeranian_Chihuahua_IcelandicShe_.430_.400_.170",
                                     "Pomeranian_Chihuahua_IcelandicShee_.430_.400_.170",
                                     "Pomeranian_Chihuahua_IcelandicSheep_.430_.400_.170",
                                     "Pomeranian_Chihuahua_IcelandicSheepd_.430_.400_.170",
                                     "Pomeranian_Chihuahua_IcelandicSheepdo_.430_.400_.170",
                                     "Poodle_MiniaturePoodle_LabradorRetriever_CockerS_.402_.349_.104_.145",
                                     "Poodle_MiniaturePoodle_LabradorRetriever_CockerSp_.402_.349_.104_.145",
                                     "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpa_.402_.349_.104_.145",
                                     "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpan_.402_.349_.104_.145",
                                     "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpani_.402_.349_.104_.145",
                                     "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpanie_.402_.349_.104_.145"
                                     ),
                                 to = c("100% (2.5x)","75% (1.875x)","50% (1.25x)",
                                        "25% (0.625x)","10% (.25x)",
                                        "5% (.125x)","100% (2.5x)","75% (1.875x)","50% (1.25x)",
                                        "25% (0.625x)","10% (.25x)",
                                        "5% (.125x)","100% (2.5x)","75% (1.875x)","50% (1.25x)",
                                        "25% (0.625x)","10% (.25x)",
                                        "5% (.125x)"))

downsampled_df$names <- mapvalues(downsampled_df$samp,from = 
                                    c("MiniaturePoodle_Poodle_LabradorRetriever_CockerS_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSp_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpa_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpan_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpani_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpanie_.422_.339_.091_.148",
                                      "Pomeranian_Chihuahua_IcelandicSh_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicShe_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicShee_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicSheep_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicSheepd_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicSheepdo_.430_.400_.170",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerS_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSp_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpa_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpan_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpani_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpanie_.402_.349_.104_.145"
                                    ),
                                  to = rep(c('Labradoodle1','PomChi','Labradoodle2'), each=6))
downsampled_df$Breed2 <- downsampled_df$Breed
for(s in unique(downsampled_df$samp)){
  downsampled_df$Breed2[downsampled_df$samp == s] <- ifelse(downsampled_df$Breed %in% downsampled_df$Breed[downsampled_df$proportion > 0 & downsampled_df$group == 'Truth' & downsampled_df$samp == s],
                                                            downsampled_df$Breed, 'AnyIncorrectBreed')
}

downsampled_df$props[downsampled_df$group == 'Truth'] <- 'Truth'
downsamp_plot <- ggplot(downsampled_df, aes(x=props,y=proportion,fill=Breed2))+geom_bar(position='fill',stat='identity')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ggtitle("Downsampling Experiment with Mixed-Breed Samples") +
  scale_x_discrete(limits=c("Truth","100% (2.5x)","75% (1.875x)","50% (1.25x)","25% (0.625x)","10% (.25x)","5% (.125x)"))+
  #facet_grid( rows = vars(group), cols = vars(names))
  facet_wrap(vars(names))
downsamp_plot
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\G_0_REF')
ggsave("downsampledplot.png",height = 9,width = 12,units = 'in',plot = downsamp_plot)

downSampling_stats$names <- mapvalues(downSampling_stats$sample,from = 
                                    c("MiniaturePoodle_Poodle_LabradorRetriever_CockerS_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSp_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpa_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpan_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpani_.422_.339_.091_.148",
                                      "MiniaturePoodle_Poodle_LabradorRetriever_CockerSpanie_.422_.339_.091_.148",
                                      "Pomeranian_Chihuahua_IcelandicSh_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicShe_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicShee_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicSheep_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicSheepd_.430_.400_.170",
                                      "Pomeranian_Chihuahua_IcelandicSheepdo_.430_.400_.170",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerS_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSp_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpa_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpan_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpani_.402_.349_.104_.145",
                                      "Poodle_MiniaturePoodle_LabradorRetriever_CockerSpanie_.402_.349_.104_.145"
                                    ),
                                  to = rep(c('Labradoodle1','PomChi','Labradoodle2'), each=6))
avgr <- aggregate(PearsonR ~ names, data=downSampling_stats, FUN = mean)

#PCA plot
#finclusters$V1 <-finclusters[,1]
#PCA_df <- merge(finclusters,eigenvec,by='V1')

#PCA_plot <- ggplot(PCA_df , aes(x=V3.y,y=V4,color=V3.x,label=V3.x))+geom_point()+
#  xlim(-0.3,0.3) +theme_bw()+labs(x='PC1',y='PC2')+ guides(fill=guide_legend(title="Breed"))+
#  geom_text_repel(data = distinct(PCA_df,V3.x,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 50,max.iter = 1000)
#PCA_plot
#PCA_plot_nt <- ggplot(PCA_df , aes(x=V3.y,y=V4,color=V3.x,label=V3.x))+geom_point()+
#  xlim(-0.3,0.3) +theme_bw() + labs(x='PC1',y='PC2') +
#  guides(fill=guide_legend(title="Breed"))#+geom_text_repel(data = distinct(PCA_df,V3.x,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 50,max.iter = 1000)
#PCA_plot_nt
#ggsave("PCA_plot.png",height = 9,width = 12,units = 'in',plot = PCA_plot)
#ggsave("PCA_plot_nt.png",height = 9,width = 12,units = 'in',plot = PCA_plot_nt)
