imma=read.table("/home/enza/oogaProtocol/IMMA/2_arraycgh/array2/all.arraychr.head.tsv.forCopynumber", header=T , sep="\t" )

#### remove duplicates  (artifact from this particular experiment)
imma.noduplicat <- imma %>% distinct(chr, start, as_sample , .keep_all = TRUE) 

toexclude <- c("AS015_bad", "AS030_bad", "AS036_bad",  "AS065_bad", "AS078_bad" ,  "AS080_bad" , "AS093_bad") # , "AS032_3xchr22" ) 
# toexclude <- c("cicci" )
#### spread
imma.spread<- imma.noduplicat %>% filter( !(as_sample  %in% toexclude)  )   %>%  spread(as_sample , LogRatio )

## 1. WINSORIZATION 
###to save in files 
#imma.win <- winsorize(imma.spread, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5,k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4,return.outliers = FALSE, save.res = TRUE, file.names =c("all.arraychr.head.tsv.forCopynumber.winsorized", "all.arraychr.head.tsv.forCopynumber.outliers") ,verbose = TRUE) 
### to continue protocol 
imma.win <- winsorize(imma.spread, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5,k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4,return.outliers = FALSE, save.res = FALSE  ,verbose = TRUE) 

## 2. SEGMENTATION 
######## DO NOT USE NORMALIZATION!!!! 
### to save files 
#imma.segments <- pcf(data=imma.win, gamma=10, assembly="hg19", return.est=TRUE, save.res=TRUE , file.names=c("imma.win.pcf", "imma.win.segments"))
### to continue protocol ######## DO NOT USE NORMALIZATION!!!! 
imma.gamma=10 
imma.segments <- pcf(data=imma.win, gamma=imma.gamma , assembly="hg19", return.est=TRUE, save.res=FALSE,  normalize = FALSE )

summary(imma.segments$segments$mean) 
sd(imma.segments$segments$mean) 

## 3. CALLING 
imma.thr.gain= mean(imma.segments$segments$mean)+3*sd(imma.segments$segments$mean)
imma.thr.loss= mean(imma.segments$segments$mean)-3*sd(imma.segments$segments$mean)
plotAberration(segments=imma.segments, thres.gain=imma.thr.gain , thres.loss =imma.thr.loss)

imma.cnvCalls <- callAberrations(imma.segments$segments, imma.thr.gain, imma.thr.loss) 
summary((imma.cnvCalls$end.pos-imma.cnvCalls$start.pos)/1000000)
min((imma.cnvCalls$end.pos-imma.cnvCalls$start.pos))
min((imma.cnvCalls$n.probes))

imma.cnvCalls.gainloss <- imma.segments$segments %>% filter (mean<imma.thr.gain |  mean > imma.thr.gain )
imma.cnvCalls.gainloss$type="PLS"

## 4. COMPARISON WITH AGILENT 
## format agilent reference calls and add to seg compare 
mc=read.table("../array2/all.cyto.tsv.forcomparison", header=T , sep="\t" )
seg.agilent= cbind.data.frame(sampleID=mc$sampleid, chrom=mc$Chr,  arm=as.character(mc$Chr) ,  start.pos=mc$Start,  end.pos=mc$Stop_bp,  n.probes=as.numeric(mc$Probes), mean=mc$Amp.Gain.Loss.Del) 
seg.agilent$type="Agilent"

##  rbind 
all.cnvCalls=rbind(imma.cnvCalls.gainloss, seg.agilent) 
ggplot(all.cnvCalls, aes((end.pos-start.pos)/1000000, n.probes, color=type ) )+geom_point(alpha=0.6 ) +theme_bw() +facet_grid(type ~ . )



