imma=read.table("/home/enza/oogaProtocol/IMMA/2_arraycgh/array2/all.arraychr.head.tsv.forCopynumber", header=T , sep="\t" )

#### remove duplicates  (artifact from this particular experiment)
imma.noduplicat <- imma %>% distinct(chr, start, as_sample , .keep_all = TRUE) 

toexclude <- c("AS015_bad", "AS030_bad", "AS036_bad",  "AS065_bad", "AS078_bad" ,  "AS080_bad" , "AS093_bad", "AS032_3xchr22" ) 

#### spread
imma.spread<- imma.noduplicat %>% filter( !(as_sample  %in% toexclude)  )   %>%  spread(as_sample , LogRatio )



###to save in files 
#imma.win <- winsorize(imma.spread, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5,k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4,return.outliers = FALSE, save.res = TRUE, file.names =c("all.arraychr.head.tsv.forCopynumber.winsorized", "all.arraychr.head.tsv.forCopynumber.outliers") ,verbose = TRUE) 
### to continue protocol 
imma.win <- winsorize(imma.spread, pos.unit = "bp", arms = NULL, method = "mad", tau = 2.5,k = 25, gamma = 40, iter = 1, assembly = "hg19", digits = 4,return.outliers = FALSE, save.res = FALSE  ,verbose = TRUE) 

### to save files 
#imma.segments <- pcf(data=imma.win, gamma=10, assembly="hg19", return.est=TRUE, save.res=TRUE , file.names=c("imma.win.pcf", "imma.win.segments"))
### to continue protocol 
imma.segments <- pcf(data=imma.win, gamma=10, assembly="hg19", return.est=TRUE, save.res=FALSE )

summary(imma.segments$segments$mean) 
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-0.629600 -0.055500 -0.004500 -0.007237  0.048500  0.570400 
sd(imma.segments$segments$mean) 
[1] 0.1153323

imma.thr.gain= mean(imma.segments$segments$mean)+4*sd(imma.segments$segments$mean)
imma.thr.loss= mean(imma.segments$segments$mean)-4*sd(imma.segments$segments$mean)
plotAberration(segments=imma.segments, thres.gain=imma.thr.gain , thres.loss =imma.thr.loss)
