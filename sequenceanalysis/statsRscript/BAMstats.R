library(knitr)
library(tidyverse)
library(reshape2)
library(gridExtra)
library(ggplot2)

## Open bam file 
data <- readLines("/home/silvia/misc/align/bcf-stats/AS074.raw.sorted.stats.bam")

## Subset bam file and plot

#Summary numbers
sn <- grep("^SN",data, value=TRUE)
sn <- separate(data.frame(sn),col=1, into=c("ID", "Name","Value"), sep="\t")[,-1]
kable(sn, caption="Summary numbers")

#Insert sizes
is <- grep("^IS",data, value=TRUE)
is <- separate(data.frame(is),col=1, into=c("ID", "insert size","all pairs", "inward", "outward", "other"), sep="\t")[,-1]
ism <- melt(is,is="size")

g<-ggplot(data = is, aes(as.numeric(get("insert size")))) + geom_line(aes(y=as.numeric(get("inward"))),color="blue") +  geom_line(aes(y=as.numeric(get("outward"))),color="red") + geom_line(aes(y=as.numeric(get("other"))),color="green") + labs( x = "insert size", y = "all pairs", title ="Mapped insert sizes", subtitle = "All Pairs", caption = "all pairs insert size")
ggsave("insert_size.png", plot= g, device="png", width = 20, height = 15, units = "cm", dpi = 300)


#ACGT content per cycle
actg <- grep("^GCC",data, value=TRUE)
actg <- separate(data.frame(actg),col=1, into=c("ID", "cycle", "A", "C", "G", "T", "N", "O"), sep="\t")[,-1]
actgm <- melt(actg,id="cycle")

ic <- ggplot(actgm, aes(as.numeric(cycle), as.numeric(value), by=variable, colour=variable)) + geom_boxplot() + labs(x= "cycle", y= "count", title= "Base composition(AS074)")
ggsave("base_comp.png", plot= ic, device="png", width = 20, height = 15, units = "cm", dpi = 300)


#First Fragment Qualitites
q <- grep("^FFQ|^LFQ",data, value=TRUE)
fq <- separate(data.frame(fq),col=1, into=c("Pair", "Cycle", seq(43)), sep="\t")

#GC Content 
gc <- grep("^GCF|^GCL",data, value=TRUE)
gc <- separate(data.frame(gc),col=1, into=c("Pair", "GC", "Count"), sep="\t")

h<-ggplot(gc, aes(as.numeric(GC), as.numeric(Count)/sum(as.numeric(Count)),color=Pair)) + geom_line()  + labs( x = "GC", y = "percent", title ="GC content(AS074)")
ggsave("GCcontent.png", plot= h, device="png", width = 20, height = 15, units = "cm", dpi = 300)


#Indel distribution
id <- grep("^ID",data, value=TRUE)
id <- separate(data.frame(id),col=1, into=c("ID", "length", "insertion_count", "deletion_count"), sep="\t")[,-1]
id$length <- as.numeric(as.character(id$length))
idm <- melt(id,id="length")

k<-ggplot(idm, aes(x=length, y=value, color=variable)) + geom_line(size=1.5) + labs (x= "length", y= "indels", title="Indel distribution(AS074)") 
ggsave("indels.png", plot= k, device="png", width = 20, height = 15, units = "cm", dpi = 300)


#Indels per cycle
ic <- grep("^IC",data, value=TRUE)
ic <- separate(data.frame(ic),col=1, into=c("ID", "cycle", "ins_fwd", "ins_rev", "del_fwd", "del_rev"), sep="\t")[,-1]

#Coverage distribution
cov <- grep("^COV",data, value=TRUE)
cov <- separate(data.frame(cov),col=1, into=c("ID", "coverage_range", "coverage", "bases"), sep="\t")[,-1]

cov$coverage<-as.numeric(as.character(cov$coverage))
cov$bases<-as.numeric(as.character(cov$bases))

pcov<-ggplot(cov, aes(x = coverage, y = bases)) + geom_bar(stat = "identity")  + coord_cartesian(xlim=c(0,150)) + labs(title= "Coverage distribution(AS074)")
ggsave("coverageAS074.png", plot= pcov, device="png", width = 20, height = 15, units = "cm", dpi = 300)

#Make panel
myplot<- grid.arrange(pcov,pcov90,pcov94,pcov06,pcov54,pcov64, nrow = 2)
ggsave("coverage.png", plot = myplot, dpi=300, units="cm", width=35, height =30)

