require(ggplot2)
#library(readr)
#library(ggrepel)
require(dplyr)
require(RColorBrewer)
library(GenomicRanges)
library(copynumber)
library(DNAcopy)
library(tidyr)


gg.manhattan <- function(df, thresholdGain, thresholdLoss, ylims, title){
#https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function
#Format GWAS results as you would for qqman: SNP CHR BP P (tab-delim)

  #~~~~~~~~~ format df
  df.tmp <- df %>%
   
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len=max(binMidPoint)) %>%
   
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
   
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
   
    # Add a cumulative position of each SNP
    arrange(CHR, binMidPoint) %>%
    mutate( binMidPointcum=binMidPoint+tot) %>%
   
    # Add highlight and annotation information
     mutate( is_annotate=ifelse(Zcov < thresholdGain | Zcov > thresholdLoss, "yes", "no"))
 
  #~~~~~~~ get chromosome center positions for x-axis
  axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(binMidPointcum) + min(binMidPointcum) ) / 2 )
 
  #~~~~~~~ ggplot
  ggplot(df.tmp, aes(x=binMidPointcum, y=Zcov)) +
 
     # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
   
     
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis
   
    # add plot and axis titles
    ggtitle(paste0(title)) +
    labs(x = "Chromosome") +
   
    # add threshold for gain and loss
    geom_hline(yintercept = c(thresholdGain, thresholdLoss), color='coral') +
    #geom_hline(yintercept = -log10(sugg), linetype="dashed") +
   
    # Add highlighted points
    #geom_point(data=subset(df.tmp, is_annotate=="yes"), color="orange", size=2) +
   
    # Add label using ggrepel to avoid overlapping
    ######## geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3) +
   
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}


gwZcov.graph <- function(df, na.rm = TRUE,  imma.thr.gain, imma.thr.loss,  ylim){
 
  # create list of id to loop over
  id_list <- unique(df$sampleID)
 
  # create for loop to produce ggplot2 graphs
  for (i in id_list) {
   
    # create plot for each id in df
   
    myp<- gg.manhattan (subset(df, sampleID==i), imma.thr.gain, imma.thr.loss,  ylim, title=i)
   
    ggsave(paste("gw.", i, ".png"), plot = myp, scale = 1, width = 40, height = 15 ,units = "cm", dpi = 300)
  }
}



######  INPUT imma.segments=segmentazione fatta con copynumber
# sampleID chrom arm start.pos end.pos n.probes    mean
# AS006     1   p    120908  867593       11 32.2964


#1) open file with coverage for all sample 
imma<-read.table("/home/silvia/misc/coverageCNV/input/allsamples.10.mean.tsv", header=T , sep="\t")
imma<-na.omit(imma)
#imma$CHR<-as.character(imma$CHR)
#imma$CHR<-as.numeric(strsplit(imma$CHR, split='chr')[[1]][2])

#2) WINDSORIZATION
imma.win <- winsorize(imma, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5,k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4,return.outliers = FALSE, save.res = FALSE  ,verbose = TRUE)

#3)Define gamma
imma.gamma=10
imma.kmin=5

#4) Segmentation
imma.segments <- pcf(data=imma.win,Y=imma,  gamma=imma.gamma , assembly="hg19", return.est=FALSE, save.res=FALSE,  normalize = FALSE)

#5)
summary(imma.segments$mean)

#6) Define treshold
nbsd=3
imma.thr.gain= mean(imma.segments$mean)+nbsd*sd(imma.segments$mean)
imma.thr.loss= mean(imma.segments$mean)-nbsd*sd(imma.segments$mean)


######  INPUT imma.segments=segmentazione fatta con copynumber
# sampleID chrom arm start.pos end.pos n.probes    mean
# AS006     1   p    120908  867593       11 32.2964

gwALL <- imma.segments %>% group_by(sampleID, chrom, start.pos, end.pos)  %>%  mutate(segSize=end.pos-start.pos , binName=paste(chrom, start.pos+50, sep="_") , CHR=chrom, binMidPoint=(end.pos-start.pos)/2+start.pos, Zcov=(mean - mean(imma.segments$mean))/sd(imma.segments$mean), pval=2*pnorm(-abs(Zcov), nprobes=n.probes ) %>% select(segSize, binName, CHR, binMidPoint, Zcov,pval, nprobes)

gwALL$CHR<- factor(gwALL$CHR, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12","13","14","15","16","17","18","19","20","21","22","X", "Y"))

nbsd=2
imma.thr.gain= mean(gwALL$Zcov)+nbsd*sd(gwALL$Zcov)
imma.thr.loss= mean(gwALL$Zcov)-nbsd*sd(gwALL$Zcov)
ylim=c(imma.thr.loss-1.5, imma.thr.gain+1.5)

write.table(gwALL, "/home/silvia/misc/coverageCNV/table/depth_mean_03X.tsv", quote=F, sep="\t", col.names=T, row.names=F)

gwZcov.graph(gwALL, na.rm = TRUE,  imma.thr.gain, imma.thr.loss,  ylim)
