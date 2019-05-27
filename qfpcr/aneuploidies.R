myd=read.table("qfarray_summary.tsv_version1", header=T , sep="\t" ,  na ="na" )
ggplot(subset(myd, !is.na(qf_outcome )  |  qf_type!="contaminata") , aes(qf_type) ) +geom_bar(stat="count")+ coord_flip() +theme_bw()
