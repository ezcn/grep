############### combine ARRAY, PROBES AND RANDOM #############################3

#1) array outcome file

myd<-read.table("/home/silvia/misc/coverageCNV/intervals/allcyto.hg38.bed", sep="\t", header=F)
colnames(myd)<-c("chrom", "start.pos", "end.pos", "sampleID", "Size_kb", "cytoband", "nprobes", "Zcov", "pval")
arr<-myd %>% select(sampleID, chrom, start.pos, end.pos, Zcov, pval, nprobes)
arr$type<-"array"
arr$coverage<-"array"
arr$variables="array"

arr$chrom<-as.character(arr$chrom)

#2) depth files (whith different coverage) gamma=10 kmin=5

setwd("/home/silvia/misc/coverageCNV/table")

probe<- data.frame(sampleID=character(), chrom=character(), start.pos=numeric(), end.pos=numeric(), Zcov=numeric(), pval=numeric(), nprobes=numeric(), type=character(), coverage=character(), variables=character(), stringsAsFactors=FALSE)


for (ss in c("30X", "15X", "1X")){
   filename=paste(ss, ".probes.depth.tsv", sep="") 
   myd=read.table(filename, header=T)
    myd$type="probes"
    myd$coverage=ss
    myd$variables="g10k5"
    mydf<- myd %>% select(sampleID,chrom,start.pos,end.pos,Zcov,pval,nprobes,type,coverage, variables)
    probe=rbind(probe, mydf)
   
}



random360<- data.frame(sampleID=character(), chrom=character(), start.pos=numeric(), end.pos=numeric(), Zcov=numeric(), pval=numeric(), nprobes=numeric(), type=character(), coverage=character(), variables=character(), stringsAsFactors=FALSE)


for (ss in c("30X", "15X", "1X")){
   filename=paste(ss, ".random360k.depth.tsv", sep="") 
   myd=read.table(filename, header=T)
    myd$type="random360k"
    myd$coverage=ss
    myd$variables="g10k5"
    mydf<- myd %>% select(sampleID,chrom,start.pos,end.pos,Zcov,pval,nprobes,type,coverage,variables)
    random360=rbind(random360, mydf)
   
}



random60<- data.frame(sampleID=character(), chrom=character(), start.pos=numeric(), end.pos=numeric(), Zcov=numeric(), pval=numeric(), nprobes=numeric(), type=character(), coverage=character(), variables=character(), stringsAsFactors=FALSE)


for (ss in c("30X", "15X", "1X")){
   filename=paste(ss, ".random60k.depth.tsv", sep="") 
   myd=read.table(filename, header=T)
    myd$type="random60k"
    myd$coverage=ss
    myd$variables="g10k5"
    mydf<- myd %>% select(sampleID,chrom,start.pos,end.pos,Zcov,pval,nprobes,type,coverage,variables)
    random60=rbind(random60, mydf)
   
}

#total<-rbind(probe,random,arr)



################################# Different gamma and kmin ######################################


#3) depth files with different segmentation (gamma=2 ;kmin=5)

probe2<- data.frame(sampleID=character(), chrom=character(), start.pos=numeric(), end.pos=numeric(), Zcov=numeric(), pval=numeric(), nprobes=numeric(), type=character(), coverage=character(), variables=character(), stringsAsFactors=FALSE)


for (ss in c("30X", "15X", "1X")){
   filename=paste(ss, ".probes.g2k5.depth.tsv", sep="") 
   myd=read.table(filename, header=T)
    myd$type="probes"
    myd$coverage=ss
    myd$variables="g2k5"
    mydf<- myd %>% select(sampleID,chrom,start.pos,end.pos,Zcov,pval,nprobes,type,coverage,variables)
    probe2=rbind(probe2, mydf)
   
}



random360_2<- data.frame(sampleID=character(), chrom=character(), start.pos=numeric(), end.pos=numeric(), Zcov=numeric(), pval=numeric(), nprobes=numeric(), type=character(), coverage=character(), variables=character(), stringsAsFactors=FALSE)


for (ss in c("30X", "15X", "1X")){
   filename=paste(ss, ".random360k.g2k5.depth.tsv", sep="") 
   myd=read.table(filename, header=T)
    myd$type="random360k"
    myd$coverage=ss
    myd$variables="g2k5"
    mydf<- myd %>% select(sampleID,chrom,start.pos,end.pos,Zcov,pval,nprobes,type,coverage,variables)
    random360_2=rbind(random360_2, mydf)
   
}

random60_2<- data.frame(sampleID=character(), chrom=character(), start.pos=numeric(), end.pos=numeric(), Zcov=numeric(), pval=numeric(), nprobes=numeric(), type=character(), coverage=character(), variables=character(), stringsAsFactors=FALSE)


for (ss in c("30X", "15X", "1X")){
   filename=paste(ss, ".random60k.g2k5.depth.tsv", sep="") 
   myd=read.table(filename, header=T)
    myd$type="random60k"
    myd$coverage=ss
    myd$variables="g2k5"
    mydf<- myd %>% select(sampleID,chrom,start.pos,end.pos,Zcov,pval,nprobes,type,coverage,variables)
    random60_2=rbind(random60_2, mydf)
   
}

#~~~~~~~~~~~~~~~~~~~~merge all files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
total<-rbind.data.frame(arr,probe,random60,random360,probe2,random60_2,random360_2)

total$coverage <- factor(total$coverage, levels = c("30X", "15X", "1X", "array"))
total$type <- factor(total$type, levels=c("array", "probes", "random60k", "random360k"))



mycol<-c("#090088", "#b0deff", "#f54291","#ffc6c7" )    # define colors for plot (ordine=blu,azzurro,rasa scuro,rosa chiaro)

ggplot(total, aes(end.pos-start.pos, fill=type, linetype=variables ) ) +geom_density(aes(alpha=0.4 ) )+facet_grid (coverage ~ . )+ scale_x_log10( ) + scale_fill_manual(values=mycol)


write.table(total, "/home/silvia/misc/coverageCNV/table/complete_CoverageVsArray.tsv", quote=F, sep="\t", row.names=F, col.names=T)

#ggplot(total2, aes(end.pos-start.pos, fill=type ) ) +geom_density(aes(alpha=0.4 ) )+facet_grid (coverage ~ . )+ scale_x_log10( )

#~~~~~~~~~~~~~~~~~~~~~~ compare

ggplot(subset(total, chrom==12 ) , aes(x=start.pos, y=Zcov, xend=end.pos, yend=Zcov) )+ geom_segment(aes(color=type  )) + ylim (-1, 2.5 )+facet_wrap(sampleID ~ . )



new<-total %>% group_by(type,variables) %>% mutate(plot=(ifelse(type=="array", "yes", ifelse(Zcov>quantile(Zcov, ,95) |  (Zcov<quantile(Zcov, ,5), "yes", "no"))))


#new<-total %>% group_by(type,variables) %>% mutate(plot=(ifelse(type=="array", "yes", ifelse(Zcov>quantile(Zcov, 95 & 5), "yes", "no")))) %>% select(sampleID, chrom, start.pos, end.pos, Zcov, pval, nprobes,type, coverage, variables, plot)



