
library(tidyverse)

metadata = read.table("regions.metadata.txt", sep="\t" , header=F)

pcaTable = read.table("pca_chr22.txt", sep="\t" , header=F) #### pca.txt è file output di akt

myd=merge(pcaTable, metadata ,"V1")

write.table(myd, "myd.txt", row.names=F, col.names=F, quote=F)

mydPCA=read.table("myd.txt", header=F)

ggplot(mydPCA, aes(V2, V3, color=V22))  + geom_point() + scale_color_manual(values=c("GREP" = "red", "EUROPE" = "#ef962d", "AFRICA" = "#f3c623", "AMERICA" = "#687466", "CENTRAL_SOUTH_ASIA" = "#be79df", "EAST_ASIA" = "#e6739f", "MIDDLE_EAST" = "#8cba51", "OCEANIA" = "#550a46"))

ggsave("pca_chr22.pdf")