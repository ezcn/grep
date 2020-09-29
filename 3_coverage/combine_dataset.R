############### combine ARRAY, PROBES AND RANDOM files#############################3

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

#~ Probes (same position and distribution of array CGH probes)

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


#~ 360k random probes

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

#~ 60k random probes

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


################################# Different gamma and kmin ######################################


#3) depth files with different segmentation (gamma=2 ;kmin=5)

#~ Probes

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

#~ 360k random probes

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

#~ 60k random probes

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

total$coverage <- factor(total$coverage, levels = c("30X", "15X", "1X", "array")) #new total$coverage order
total$type <- factor(total$type, levels=c("array", "probes", "random60k", "random360k")) #new total$type order



mycol<-c("#090088", "#b0deff", "#f54291","#ffc6c7" )    # define colors for plot (ordine=blu,azzurro,rasa scuro,rosa chiaro)

ggplot(total, aes(end.pos-start.pos, fill=type, linetype=variables ) ) +geom_density(aes(alpha=0.4 ) )+facet_grid (coverage ~ . )+ scale_x_log10( ) + scale_fill_manual(values=mycol)


write.table(total, "/home/silvia/misc/coverageCNV/table/complete_CoverageVsArray.tsv", quote=F, sep="\t", row.names=F, col.names=T)



#~~~~~~~~~~~~~~~~~~~~~~ compare

ggplot(subset(total, chrom==12 ) , aes(x=start.pos, y=Zcov, xend=end.pos, yend=Zcov) )+ geom_segment(aes(color=type  )) + ylim (-1, 2.5 )+facet_wrap(sampleID ~ . )

# new variable to plot only values >0.95 quantile e <0.05 

new<-total %>% filter(sampleID=="AS054" | sampleID=="AS064" |sampleID=="AS074" |sampleID=="AS094" |sampleID=="AS090" ) %>% group_by(type,newvar) %>% mutate(plot=(ifelse(type=="array", "yes", ifelse(Zcov>quantile(Zcov, .95) | Zcov<quantile(Zcov, .05), "yes", "no"))))

new$plotf=as.factor(new$plot)

ggplot(subset(new, plotf=="yes" & chrom==8 & sampleID=="AS064" ) , aes(x=start.pos, y=Zcov, xend=end.pos, yend=Zcov) )+ geom_segment(aes(color=type  )) +  facet_wrap(sampleID ~ . ) + geom_hline(yintercept=0, color="grey" )+theme_minimal()
