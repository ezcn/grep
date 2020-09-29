# DEPTH (quality not filtered)

cat AS090.chr22.bcf-stats | grep "^DP\|^# DP" | grep -v "DP," > AS090.chr22.depth.tsv

## Open vcf files
myd74q<- read.table("/home/silvia/misc/align/tsvStats/AS074.chr22.depth.tsv", header=T, sep="\t")
myd90q<- read.table("/home/silvia/misc/align/tsvStats/AS090.chr22.depth.tsv", header=T, sep="\t")
myd94q<- read.table("/home/silvia/misc/align/tsvStats/AS094.chr22.depth.tsv", header=T, sep="\t")
myd06q<- read.table("/home/silvia/misc/align/tsvStats/AS006.chr22.depth.tsv", header=T, sep="\t")
myd54q<- read.table("/home/silvia/misc/align/tsvStats/AS054.chr22.depth.tsv", header=T, sep="\t")
myd64q<- read.table("/home/silvia/misc/align/tsvStats/AS064.chr22.depth.tsv", header=T, sep="\t")


## Create ID column
myd74q$ID<-"AS074"
myd90q$ID<-"AS090"
myd94q$ID<-"AS094"
myd06q$ID<-"AS006"
myd54q$ID<-"AS054"
myd64q$ID<-"AS064"

## Merge datasets
mydq<- rbind(myd74q,myd90q,myd94q,myd06q,myd54q,myd64q)

## DEPTH (filtered quality >20)

## Open vcf files
myd74<- read.table("/home/silvia/misc/align/tsvStats/AS074.chr22QUAL.tsv", header=T, sep="\t")
myd90<- read.table("/home/silvia/misc/align/tsvStats/AS090.chr22QUAL.tsv", header=T, sep="\t")
myd94<- read.table("/home/silvia/misc/align/tsvStats/AS094.chr22QUAL.tsv", header=T, sep="\t")
myd06<- read.table("/home/silvia/misc/align/tsvStats/AS006.chr22QUAL.tsv", header=T, sep="\t")
myd54<- read.table("/home/silvia/misc/align/tsvStats/AS054.chr22QUAL.tsv", header=T, sep="\t")
myd64<- read.table("/home/silvia/misc/align/tsvStats/AS064.chr22QUAL.tsv", header=T, sep="\t")

## Create ID column
myd74$ID<-"AS074"
myd90$ID<-"AS090"
myd94$ID<-"AS094"
myd06$ID<-"AS006"
myd54$ID<-"AS054"
myd64$ID<-"AS064"

## Merge datasets
myd<- rbind(myd74,myd90,myd94,myd06,myd54,myd64)

## Create quality column
mydq$qual<-">0"
myd$qual<-">20"

## Merge the two datasets
myst<- rbind(mydq, myd)

## Plot and save
depqual<-ggplot(myst, aes(bin, number.of.sites) )+ geom_jitter(aes(color=qual)) + ggtitle("Depth")
ggsave("depth-quality.png", plot= depqual, device="png", width = 20, height = 15, units = "cm", dpi = 300)
