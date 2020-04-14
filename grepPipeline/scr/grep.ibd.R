library(ggplot2)
library(ggsci) 
library(tidyr)
library(dplyr)


idlist=read.table('ID_category.tsv', header=F, col.names=c('id', 'category'), sep='\t')


# homozygosity by descent, runs of homozygosity 
temphom =read.table("memorial_exome.hbd", header=F , sep="\t",col.names=c('id1', 'hapid1' ,'id2', 'hapid2' , 'chr', 'start', 'end', 'lod', 'cM') )

hbd=merge(temphom, idlist, by.x='id1', by.y='id')

hbd %>% mutate(lenbp= end-start )%>%  group_by(id1, category) %>% summarize (totROH=sum(lenbp)) %>% ggplot(aes(id1,totROH/1000000 , fill=category))+geom_bar(stat='identity')+coord_flip()+theme_minimal()  +facet_wrap(category ~ ., scales='free_y' ) + ylab('genome in ROH (Mbp)') +xlab('')
ggsave('overview.png')


hbd$lenMb =(hbd$end-hbd$start)/10000000
hbd$rohtype =ifelse(hbd$lenMb<=0.5, 'ROH <= 0.5Mb', ifelse(hbd$lenMb>5, '>5 Mb', '0.5Mb > ROH =< 5Mb') ) 
 

hbd %>% filter (category=='PREIMPLANTATION DEVELOPMENT ARREST')   ggtitle('PREIMPLANTATION DEVELOPMENT ARREST')
ggsave('ROHlen.png')


temphom1=merge(temphom, idlist, by.x='id1', by.y='id')

hbd=merge(temphom1, idlist, by.x='id2', by.y='id')

#identity by descent 

tempibd=read.table('memorial_exome.ibd',header=F , sep="\t",col.names=c('id1', 'hapid1' ,'id2', 'hapid2' , 'chr', 'start', 'end', 'lod', 'cM') )

tempibd1=merge(tempibd, idlist, by.x='id1', by.y='id')
ibd=merge(temibd1, idlist, by.x='id2', by.y='id')
ibd$lenMb =(hbd$end-hbd$start)/10000000

ibd %>% filter (category.x=='PREIMPLANTATION DEVELOPMENT ARREST' & category.y=='PREIMPLANTATION DEVELOPMENT ARREST') %>% mutate (pair=paste(id1, hapid1, id2,  hapid2, sep='_')) %>% ggplot(aes(pair, lenMb, fill=lod)) + geom_bar(stat='identity') + facet_wrap(chr ~  . , scales ='free_y') +coord_flip() + scale_fill_gsea() + theme_bw()
ggsave('ibd.png')


 m4 %>% group_by(gene_symbol) %>% summarize(tot=sum(HIGH,LOW,MODERATE, na.rm=T  ))%>% filter (tot>20 ) %>% ggplot(aes(reorder(gene_symbol, tot ), tot)) +geom_bar(stat='identity')+coord_flip()


 m4 %>% group_by(gene_symbol) %>% summarize(tot=sum(HIGH,LOW,MODERATE, na.rm=T  ))%>% filter (tot>20 ) %>% ggplot(aes(reorder(gene_symbol, tot ), tot)) +geom_bar(stat='identity')+coord_flip()
