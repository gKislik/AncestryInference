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
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')
eigenvec <- read.table("filt-bcftools-ref-redone.imputed.filtered.eigenvec")
files <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\synth", pattern="*out.txt*", full.names = FALSE, no.. = TRUE)
highCovfiles <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\out_highCoverage", pattern="*out.txt*", full.names = FALSE, no.. = TRUE)
lowCovfiles <- list.files(path = "C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\out_lowCoverage", pattern="*out.txt*", full.names = FALSE, no.. = TRUE)
downSamplingfiles <- list.files(path="C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\downsampled",pattern="*out.txt",full.names=FALSE,no.. = TRUE)
freq <- read.table('hf.txt',header = T)
clust <- freq$CLST
finclusters <- read.table("used-clust.txt")
#truth <- read.csv("C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\new_truth_matrix.csv",header=F)

makedf <- function(file_list,stats_name){
  est_df <- data.frame(Breed=as.character(),Estimate=as.numeric(),Truth=as.numeric(),samp=as.character())
  colnames(est_df) <- 'Breed'
  stats_df <- data.frame(RMSE=as.numeric(),PearsonR=as.numeric(), sample=as.character())
  for(file in file_list){
    x <- read.table(file,header = T)
    colnames(x) <- c('Estimate')
    x$Breed <- rownames(x)
    name <- strsplit(file,"_out.txt")[[1]]
    print(name[[1]])
    #print(truth[truth$V1==name,])
    props <- unlist(str_split(name,"X"))
    truthvec <- as.data.frame(cbind(clust,rep(0,48)))
    truthvec[,2] <- as.numeric(truthvec[,2])
    colnames(truthvec) <- c('Breed','Truth')
    print(props)
    props <- props[nzchar(props)]
    print(props)
    for (p in props){
      b <- (str_extract(p, "[aA-zZ]+"))
      print(b)
      per <- as.numeric(str_extract(p, "[0-9]+"))/1000
      print(per)
      #have to add in case the same breed was selected multiple times
      truthvec[substr(truthvec[,1],1,10) %like% substr(b,1,10),2] <- truthvec[substr(truthvec[,1],1,10) %like% substr(b,1,10),2] + per
    }
    x <- merge(x,truthvec,by='Breed') #as.vector(truth[truth$V1==name,-c(1,2)])
    x$samp <- name  
    assign(name,x,.GlobalEnv)
    RMSE <- rmse(x$Truth, x$Estimate)
    Corr <- cor(x$Truth, x$Estimate)
    stats_vec <- data.frame(RMSE,Corr,name)
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
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\synth')
synth_df <- makedf(files,'synth_stats')
synthplot <- ggplot(synth_df, aes(x=factor(samp,level=synth_stats$sample[order(synth_stats$RMSE)]),y=proportion,fill=Breed))+geom_bar(position='stack',stat='identity')+theme_bw()+
  facet_wrap(~group,nrow=2)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab('sample')
synthplot
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\out_highCoverage')
high_df <- makedf(highCovfiles,'highCov_stats')
high_df$Coverage <- 'high'

setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\out_lowCoverage')
low_df <- makedf(lowCovfiles,'lowCov_stats')
low_df$Coverage <- 'low'

#from the 96 samples
coverage_df <- rbind(high_df,low_df)
##remove rows where breed percentages are <3% to make sure it's easier to read
coverage_df$Breed[coverage_df$proportion < 0.03] <- '<3% breeds'
covplot <- ggplot(coverage_df, aes(x=factor(samp,level=highCov_stats$sample[order(highCov_stats$RMSE)]),y=proportion,fill=Breed))+geom_bar(position='stack',stat='identity')+
  theme_bw()+facet_wrap(~group+Coverage,nrow=2)+xlab('sample')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11))
covplot

setwd('/qmi_home/kislikgr/')
ggsave("synthplot.png",height = 9,width = 12,units = 'in',plot = synthplot)
ggsave("coverageplot.png",height = 9,width = 12,units = 'in',plot = covplot)

#for downsampling experiment
setwd('/qmi_home/kislikgr/downsampled')
downsampled_df <- makedf(downSamplingfiles,'downSampling_stats')
downsampled_df$samp <- mapvalues(downsampled_df$samp,from = 
                                   c("Poodle402XMiniaturePoodle349XLabradorRetriever104XCockerS145",
                                     "Poodle402XMiniaturePoodle349XLabradorRetriever104XCockerSp145",
                                     "Poodle402XMiniaturePoodle349XLabradorRetriever104XCockerSpa145",
                                     "Poodle402XMiniaturePoodle349XLabradorRetriever104XCockerSpan145",
                                     "Poodle402XMiniaturePoodle349XLabradorRetriever104XCockerSpani145",
                                     "Poodle402XMiniaturePoodle349XLabradorRetriever104XCockerSpanie145"
                                   ),
                                 to = c("100% (2.5x)","75% (1.875x)","50% (1.25x)",
                                        "25% (0.625x)","10% (.25x)",
                                        "5% (.125x)"))
downsamp_plot <- ggplot(downsampled_df, aes(x=samp,y=proportion,fill=Breed))+geom_bar(position='stack',stat='identity')+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ggtitle("Downampling Experiment with Mixed-Breed Sample") +
  scale_x_discrete(limits=c("100% (2.5x)","75% (1.875x)","50% (1.25x)",
                            "25% (0.625x)","10% (.25x)",
                            "5% (.125x)"))
downsamp_plot
setwd('C:\\Users\\gkisl\\Downloads\\Dog WGS - Pellegrini\\')
ggsave("downsampledplot.png",height = 9,width = 12,units = 'in',plot = downsamp_plot)

#PCA plot
finclusters$V1 <-finclusters[,1]
PCA_df <- merge(finclusters,eigenvec,by='V1')

PCA_plot <- ggplot(PCA_df , aes(x=V3.y,y=V4,color=V3.x,label=V3.x))+geom_point()+
  xlim(-0.3,0.3) +theme_bw()+labs(x='PC1',y='PC2')+ guides(fill=guide_legend(title="Breed"))+
  geom_text_repel(data = distinct(PCA_df,V3.x,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 50,max.iter = 1000)
PCA_plot
PCA_plot_nt <- ggplot(PCA_df , aes(x=V3.y,y=V4,color=V3.x,label=V3.x))+geom_point()+
  xlim(-0.3,0.3) +theme_bw() + labs(x='PC1',y='PC2') +
  guides(fill=guide_legend(title="Breed"))#+geom_text_repel(data = distinct(PCA_df,V3.x,.keep_all = TRUE),force_pull = 0.01,force = 20,max.overlaps = 50,max.iter = 1000)
PCA_plot_nt
ggsave("PCA_plot.png",height = 9,width = 12,units = 'in',plot = PCA_plot)
ggsave("PCA_plot_nt.png",height = 9,width = 12,units = 'in',plot = PCA_plot_nt)
