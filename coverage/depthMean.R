brary(tidyr)
library(dplyr)
allsamples <- data.frame(CHR=character(), BP=numeric(), DP=numeric(), BIN=numeric(), stringsAsFactors=FALSE)

setwd("/mpba0/vcolonna/silvia/coverage")
for (ss in c("AS006", "AS054", "AS064", "AS074", "AS090", "AS094")){
   filename=paste(ss, ".full.depthBin.bed", sep="")
   myd=read.table(filename, header=F )
    myd$ID=ss
    mydf<- data.frame(cbind(CHR=myd$V1, BP=myd$V2, DP=myd$V4, ID=myd$ID , BIN=myd$V8) )
    allsamples=rbind(allsamples, mydf)

}

allsamples$BP<-as.numeric(as.character(allsamples$BP))
allsamples$DP<-as.numeric(as.character(allsamples$DP))
allsamples$BIN<-as.numeric(as.character(allsamples$BIN))

myd<-allsamples %>% group_by(CHR,BIN,ID) %>% summarize(Mean = mean(DP))

mydM<-spread(myd,ID,Mean)

write.table(mydM ,"/mpba0/vcolonna/silvia/coverage/allsamples.depth.mean.tsv", quote=F, sep="\t", row.names=F)

