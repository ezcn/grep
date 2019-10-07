library(GenomicRanges)
library(copynumber)
library(DNAcopy)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

#1) open file with coverage for all sample 
imma<-read.table("allsamples.depth.mean.tsv", header=T , sep="\t")

#2) WINDSORIZATION
imma.win <- winsorize(imma, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5,k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4,return.outliers = FALSE, save.res = FALSE  ,verbose = TRUE)

#3)Define gamma
imma.gamma=10

#4) Segmentation
imma.segments <- pcf(data=imma.win,Y=imma,  gamma=imma.gamma , assembly="hg19", return.est=FALSE, save.res=FALSE,  normalize = FALSE)

#5)
summary(imma.segments$mean)

#6) Define treshold
nbsd=3
imma.thr.gain= mean(imma.segments$mean)+nbsd*sd(imma.segments$mean)
imma.thr.loss= mean(imma.segments$mean)-nbsd*sd(imma.segments$mean)

#7)Plot data
png ("imma.cnv.pls.png", res=300, width=25 ,height=10, units="cm") 
plotAberration(segments=imma.segments, thres.gain=imma.thr.gain , thres.loss =imma.thr.loss)
dev.off()



























