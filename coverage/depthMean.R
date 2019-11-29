library(tidyr)
library(dplyr)


allsamples <- data.frame(CHR=character(), BP=numeric(), DP=numeric(), BIN=numeric(), stringsAsFactors=FALSE)

setwd("/mpbastudies3/IMMA/samples/coverage/files")
for (ss in c("AS006", "AS054", "AS064", "AS074", "AS090", "AS094")){
	   filename=paste(ss, ".random360k.3X.depthBin.bed", sep="")
   myd=read.table(filename, header=F, col.names=c("CHR", "START", "END", "DP", "CHR", "START", "END", "BIN") )
       myd$ID=ss
       #mydf<- data.frame(cbind(CHR=myd$V1, BP=myd$V2, DP=myd$V4, ID=myd$ID , BIN=myd$V8))
       mydf<- myd %>% select(CHR,START,DP,BIN,ID)
           allsamples=rbind(allsamples, mydf)

}

#allsamples$BP<-as.numeric(as.character(allsamples$BP))
#allsamples$DP<-as.numeric(as.character(allsamples$DP))
#allsamples$BIN<-as.numeric(as.character(allsamples$BIN))

myd<-allsamples %>% group_by(CHR,BIN,ID) %>% summarize(Mean = mean(DP))

mydM<-spread(myd,ID,Mean)

write.table(mydM ,"/mpbastudies3/IMMA/samples/coverage/Rout/allsamples.random360K.3X.mean.tsv", quote=F, sep="\t", row.names=F)
