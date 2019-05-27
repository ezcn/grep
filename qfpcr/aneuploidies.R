myd=read.table("qfarray_summary.tsv_version1", header=T , sep="\t" ,  na ="na" )

ggplot(subset(myd, !is.na(qf_outcome )  |  qf_type!="contaminata") , aes(qf_type) ) +geom_bar(stat="count")+ coord_flip() +theme_bw() +ggtitle("Aneuploidies at chr XXXXXXXX ")



myqftype= subset(myd, !is.na(qf_type)) %>% filter(qf_type!="contaminata") %>%  group_by(qf_type) %>% tally() %>% mutate(ratio=scales::percent(n/sum(n)) )

myqftype2<-myqftype %>% arrange(desc(n))

pQfpcr <- ggplot(myqftype2, aes(x=reorder(qf_type, -n), y=n)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + ggtitle("Aneuploidies at chr 13,15,16,18,21,22,X, and Y") + ylab("count") + xlab("")
ggsave("qfpcrout.png", plot= pQfpcr, device="png", width = 15, height = 10, units = "cm", dpi = 300)
