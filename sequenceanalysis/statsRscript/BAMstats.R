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

#ACGT content per cycle
actg <- grep("^GCC",data, value=TRUE)
actg <- separate(data.frame(actg),col=1, into=c("ID", "cycle", "A", "C", "G", "T", "N", "O"), sep="\t")[,-1]

#First Fragment Qualitites
q <- grep("^FFQ|^LFQ",data, value=TRUE)
fq <- separate(data.frame(fq),col=1, into=c("Pair", "Cycle", seq(43)), sep="\t")

#GC Content 
gc <- grep("^GCF|^GCL",data, value=TRUE)
gc <- separate(data.frame(gc),col=1, into=c("Pair", "GC", "Count"), sep="\t")

#Indel distribution
id <- grep("^ID",data, value=TRUE)
id <- separate(data.frame(id),col=1, into=c("ID", "length", "insertion_count", "deletion_count"), sep="\t")[,-1]

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



