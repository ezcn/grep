library(ggplot2) 
library(dplyr) 
library(ggsci)
my1=read.table("freq.all.SNV1.out", header=T, sep ='\t' ) 
my05=read.table("freq.all.SNV05.out", header=T, sep ='\t' ) 
my1$thres=1
my05$thres=0.05
my=bind_rows(my1, my05)
ggplot(my, aes(type, Zmu, color=as.factor(thres)  ))+ geom_point() +  geom_errorbar(aes(ymin=Zmu-1.96*Zsd , ymax=Zmu+1.96*Zsd) , width = 0.1) +coord_flip()+ scale_color_ucscgb()+ theme_minimal() + labs(x="", y="standardized mean allele frequency in GREP", col="Maximum allele frequency\nin 1000G and gmomeAD", title="Single Nucleotide Variants")
