library(ggplot2)

qmystat74<- read.table("/home/gianluca/sequencinganalysis/stats/SN_AS074.chr22.bcfQUALstats.tsv",sep = "\t") 
qmystat90<- read.table("/home/gianluca/sequencinganalysis/stats/SN_AS090.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat94<- read.table("/home/gianluca/sequencinganalysis/stats/SN_AS094.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat06<- read.table("/home/gianluca/sequencinganalysis/stats/SN_AS006.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat54<- read.table("/home/gianluca/sequencinganalysis/stats/SN_AS054.chr22.bcfQUALstats.tsv",sep = "\t")
qmystat64<- read.table("/home/gianluca/sequencinganalysis/stats/SN_AS064.chr22.bcfQUALstats.tsv",sep = "\t")


qmystat74$id<-"AS074"
qmystat90$id<-"AS090"
qmystat94$id<-"AS094"
qmystat06$id<-"AS006"
qmystat54$id<-"AS054"
qmystat64$id<-"AS064"

stat<- rbind(qmystat74,qmystat90,qmystat94,qmystat06,qmystat54,qmystat64)

stat$qual<-"qual>20"

mystat74<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS074.chr22stat.tsv",sep = "\t") 
mystat90<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS090.chr22stat.tsv",sep = "\t")
mystat94<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS094.chr22stat.tsv",sep = "\t")
mystat06<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS006.chr22.SN.bcf-stats.tsv",sep = "\t")
mystat54<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS054.chr22.SN.bcf-stats.tsv",sep = "\t")
mystat64<- read.table("/home/gianluca/sequencinganalysis/stats/tsvStats/AS064.chr22.SN.bcf-stats.tsv",sep = "\t")

mystat74$id<-"AS074"
mystat90$id<-"AS090"
mystat94$id<-"AS094"
mystat06$id<-"AS006"
mystat54$id<-"AS054"
mystat64$id<-"AS064"

stat2<- rbind(mystat74,mystat90,mystat94,mystat06,mystat54,mystat64)

stat2$qual<-"qual>0"

stat3<- rbind(stat,stat2)


pindels<-ggplot(subset(stat3, V1=="number of indels:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("Number of indels")

ggsave("SN-numberofindels.png", plot= pstats3, device="png", width = 20, height = 15, units = "cm", dpi = 300)

precord<-ggplot(subset(stat3, V1=="number of records:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of records")

ggsave("SN-numberofrecords.png", plot= precord, device="png", width = 20, height = 15, units = "cm", dpi = 300)

pnoALT<-ggplot(subset(stat3, V1=="number of no-ALTs:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of no-ALTs")

ggsave("SN-numberofno-ALTs.png", plot= pnoALT, device="png", width = 20, height = 15, units = "cm", dpi = 300)

pnSNP<-ggplot(subset(stat3, V1=="number of SNPs:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of SNPs")

ggsave("SN-numberofSNPs.png", plot= pnSNP, device="png", width = 20, height = 15, units = "cm", dpi = 300)

pnMNP<-ggplot(subset(stat3, V1=="number of MNPs:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of MNPs")

ggsave("SN-numberofMNPs.png", plot= pnMNP, device="png", width = 20, height = 15, units = "cm", dpi = 300)

pnMAS<-ggplot(subset(stat3, V1=="number of multiallelic sites:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of multiallelic sites")

ggsave("SN-numberofmultiallelicsites.png", plot= pnMAS, device="png", width = 20, height = 15, units = "cm", dpi = 300)


pnMASS<-ggplot(subset(stat3, V1=="number of multiallelic SNP sites:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of multiallelic SNP sites")

ggsave("SN-numberofmultiallelicSNPsites.png", plot= pnMASS, device="png", width = 20, height = 15, units = "cm", dpi = 300)

pnOO<-ggplot(subset(stat3, V1=="number of others:"), aes(x = V1, y = V2, colour=qual))+ geom_jitter(stat = "identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle("number of others")

ggsave("SN-numberofothers.png", plot= pnOO, device="png", width = 20, height = 15, units = "cm", dpi = 300)







