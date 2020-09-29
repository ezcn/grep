
################# TSTV

library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)
library(tidyverse)
library(reshape2)

qmystat74<- read.table("/home/gianluca/sequencinganalysis/stats/TSTV_AS074.chr22.bcfQUALstats.tsv",sep = "\t") 
qmystat90<- read.table("/home/gianluca/sequencinganalysis/stats/TSTV_AS090.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat94<- read.table("/home/gianluca/sequencinganalysis/stats/TSTV_AS094.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat06<- read.table("/home/gianluca/sequencinganalysis/stats/TSTV_AS006.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat54<- read.table("/home/gianluca/sequencinganalysis/stats/TSTV_AS054.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat64<- read.table("/home/gianluca/sequencinganalysis/stats/TSTV_AS064.chr22.bcfQUALstats.tsv",sep = "\t")

qmystat74$id<-"AS074"
qmystat90$id<-"AS090"
qmystat94$id<-"AS094"
qmystat06$id<-"AS006"
qmystat54$id<-"AS054"
qmystat64$id<-"AS064"

stat<- rbind(qmystat74,qmystat90,qmystat94,qmystat06,qmystat54,qmystat64)
stat$qual<-"qual>20"


mystat74<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS074.chr22.TSTV.tsv",sep = "\t") 
mystat90<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS090.chr22.TSTV.tsv",sep = "\t")
mystat94<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS094.chr22.TSTV.tsv",sep = "\t")
mystat06<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS006.chr22.TSTV.bcf-stats.tsv",sep = "\t")
mystat54<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS054.chr22.TSTV.bcf-stats.tsv",sep = "\t")
mystat64<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS064.chr22.TSTV.bcf-stats.tsv",sep = "\t")

mystat74$id<-"AS074"
mystat90$id<-"AS090"
mystat94$id<-"AS094"
mystat06$id<-"AS006"
mystat54$id<-"AS054"
mystat64$id<-"AS064"

stat2<- rbind(mystat74,mystat90,mystat94,mystat06,mystat54,mystat64)
stat2$qual<-"qual>0"

stat3<- rbind(stat,stat2)

colnames(stat3) <- c("ts", "tv", "ts/tv","id","qual")

pts<-ggplot(stat3, aes(x = id, y =ts , fill=qual)) + geom_bar(stat = "identity", position= "dodge") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ ggtitle("Transition")
ptv<-ggplot(stat3, aes(x = id, y =tv , fill=qual)) + geom_bar(stat = "identity", position= "dodge") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ggtitle("Transvertion")
pratio<-ggplot(stat3, aes(x = id, y =ts/tv , fill=qual)) + geom_bar(stat = "identity", position= "dodge") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+ggtitle("TSTV ratio")


ggsave("transition.png", plot= pts, device="png", width = 8, height = 8, units = "cm", dpi = 300)
ggsave("transversion.png", plot= ptv, device="png", width = 8, height = 8, units = "cm", dpi = 300)
ggsave("tstvratio.png", plot= pratio, device="png", width = 8, height = 8, units = "cm", dpi = 300)
