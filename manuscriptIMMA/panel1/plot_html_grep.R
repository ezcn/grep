library(tidyverse)
library(ggsci)
library (grid) 
library(lattice) 
library(gridExtra)

######## All consequence 

palette<-c("#381460", "#b21f66","#b0a160","#ffbd69","#5b8c85","#434e52","#fe346e","#ecce6d","#f8b195","#f67280","#c06c84","#6c567b")

pfolder="/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/plot/img/"

file="/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/scripts/output/varCons.out"

myd=read.table(file, header=T)

impact=read.table("/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/plot/code/lsoterm_impact.tsv", header=T)

mydd=merge(myd, impact, by="type")

######## ggplot(mydd, aes(id, count, fill=type)) + geom_bar(stat='identity')+ facet_wrap(. ~ impact , scales="free_y")


myddFilt= mydd %>% spread(id, count ) %>% mutate  (sum= rowSums(.[3:8 ])) %>% filter (sum >0 ) %>% gather (key= id , value= count ,-type, -impact, -sum) # %>%  ggplot( aes(id, count, fill=type)) + geom_bar(stat='identity')+ facet_wrap(. ~ impact , scales="free_y")

######### obtain frequences from vep impact

#mydd %>% group_by(id, impact) %>% summarise (tot=sum(count)) %>% spread(impact, tot ) %>% mutate(all=sum(HIGH, LOW, MODERATE, MODIFIER) , h=HIGH/all*100, mode=MODERATE/all*100, l=LOW/all*100, , modi=MODIFIER/all*100#)
#
########## ratio of LOW consequence
#
#lowCons=mydd %>% filter(impact=="LOW") %>% group_by(id, type) %>% summarise(tot=sum(count)) %>% spread(type,tot) %>% mutate(all=sum(incomplete_terminal_codon_variant, splice_region_variant, start_retained_variant, #stop_retained_variant, synonymous_variant), incomp=incomplete_terminal_codon_variant/all*100, start=start_retained_variant/all*100, splice=splice_region_variant/all*100, stop=stop_retained_variant/all*100, #synony=synonymous_variant/all*100)
#mean(lowCons$incomp)
#mean(lowCons$start)
#mean(lowCons$splice)
#mean(lowCons$stop)
#
#
######### ration of MODIFIER consequence
#
#modifierTemp=mydd %>% filter(impact=="MODIFIER") %>% group_by(id, type) %>% summarise(tot=sum(count)) %>% spread(type,tot) 
#names(modifierTemp)[names(modifierTemp) == '3_prime_UTR_variant'] <- 'three_prime_UTR_variant'
#names(modifierTemp)[names(modifierTemp) == '5_prime_UTR_variant'] <- 'five_prime_UTR_variant'
#modifierCons=modifierTemp %>% group_by(id) %>% mutate(all=sum(three_prime_UTR_variant, five_prime_UTR_variant, NMD_transcript_variant, TFBS_ablation, TFBS_amplification, TF_binding_site_variant, #coding_sequence_variant, downstream_gene_variant, feature_elongation, feature_truncation, intergenic_variant, intron_variant, mature_miRNA_variant, non_coding_transcript_exon_variant, #non_coding_transcript_variant, regulatory_region_amplification, regulatory_region_variant, upstream_gene_variant) , intron=intron_variant/all*100, intergenic=intergenic_variant/all*100)
#mean(modifierCons$intron)
#mean(modifierCons$intergenic)

#levels(myddFiltDef$impact) <- c("HIGH","MODERATE","LOW","MODIFIER")

#levels(myddFiltDef$type)<- c("splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","splice_region_variant","incomplete_terminal_codon_variant","stop_retained_variant","synonymous_variant","coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","non_coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TF_binding_site_variant","regulatory_region_variant","intergenic_variant")
#381460", "#b21f66","#fe346e","#ffbd69","#5b8c85","#434e52","#381460", "#b21f66","#fe346e","#ffbd69","#381460", "#b21f66","#fe346e","#ffbd69",



#ggplot(subset(myddFiltDef,id='AS006'), aes(id, count, fill=type)) + geom_bar(stat='identity')#+ facet_wrap(. ~ impact , scales="free_y")# +scale_fill_manual(values=palette)

#ggplot(mydd, aes(id, count, fill=type)) + geom_bar(stat='identity')+ facet_wrap(. ~ impact , scales="free_y")

#colnames(myddFilt)[colnames(myddFilt) == "sum"] <- "total"

#myHigh=myddFilt %>% filter(impact=="HIGH") %>% group_by(id) %>%  mutate(ratio=count/sum(count))

pHigh=ggplot(subset(myddFilt,impact=='HIGH'), aes(id, count, fill=type)) + geom_bar(stat='identity') + scale_fill_npg() + theme_bw() +scale_fill_manual(values=palette) + ggtitle("A.     High") + xlab("ID") + theme(legend.title=element_blank(),legend.text=element_text(size=20),plot.title = element_text(size=20))
#ggsave(paste(pfolder,'high_cons.png',sep=''),plot=last_plot())
pModerate=ggplot(subset(myddFilt,impact=='MODERATE'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw()+scale_fill_manual(values=palette)+ ggtitle("B.     Moderate") + xlab("ID") + theme(legend.title=element_blank(),legend.text=element_text(size=20),plot.title = element_text(size=20))
#ggsave(paste(pfolder,'moderate_cons.png',sep=''),plot=last_plot())
pLow=ggplot(subset(myddFilt,impact=='LOW'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw()+scale_fill_manual(values=palette)+ ggtitle("C.     Low") + xlab("ID") + theme(legend.title=element_blank(),legend.text=element_text(size=20),plot.title = element_text(size=20))
#ggsave(paste(pfolder,'low_cons.png',sep=''),plot=last_plot())
pModifier=ggplot(subset(myddFilt,impact=='MODIFIER'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw()+scale_fill_manual(values=palette)+ ggtitle("D.     Modifier") + xlab("ID") + theme(legend.title=element_blank(),legend.text=element_text(size=20),plot.title = element_text(size=20))
#ggsave(paste(pfolder,'modifier_cons.png',sep=''),plot=last_plot())

lay <- rbind(c(1,2),
             c(3,4))
pgrid=grid.arrange(pHigh, pModerate, pLow, pModifier, nrow = 2, layout_matrix = lay)
ggsave(paste(pfolder,'grid_cons.pdf',sep=''),plot=pgrid,width = 45, height = 20, units = "cm")

################################################ Consequence in coding region

file="/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/scripts/output/varCod.out"

myd=read.table(file, header=T)

impact=read.table("/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/plot/code/lsoterm_impact.tsv", header=T)

mydd=merge(myd, impact, by="type")

myddFilt= mydd %>% spread(id, count ) %>% mutate  (sum= rowSums(.[3:8 ])) %>% filter (sum >0 ) %>% gather (key= id , value= count ,-type, -impact, -sum)

#str(myddFilt)

pHigh=ggplot(subset(myddFilt,impact=='HIGH'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw() +scale_fill_manual(values=palette)
#ggsave(paste(pfolder,'high_cons-coding.png',sep=''),plot=last_plot())

pModerate=ggplot(subset(myddFilt,impact=='MODERATE'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw()+scale_fill_manual(values=palette)
#ggsave(paste(pfolder,'moderate_cons-coding.png',sep=''),plot=last_plot())

pLow=ggplot(subset(myddFilt,impact=='LOW'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw()+scale_fill_manual(values=palette)
#ggsave(paste(pfolder,'low_cons-coding.png',sep=''),plot=last_plot())

pModifier=ggplot(subset(myddFilt,impact=='MODIFIER'), aes(id, count, fill=type)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw()+scale_fill_manual(values=palette)
#ggsave(paste(pfolder,'modifier_cons-coding.png',sep=''),plot=last_plot())


lay <- rbind(c(1,2),
             c(3,4))
pgrid=grid.arrange(pHigh, pModerate, pLow, pModifier, nrow = 2, layout_matrix = lay)
ggsave(paste(pfolder,'grid_cons-coding.pdf',sep=''),plot=pgrid,width = 35, height = 20, units = "cm")

file="/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/scripts/output/varClass.out"

myd=read.table(file, header=T)

### obtain ration of variant classes and means among samples

ratioAndMean=myd %>% group_by(id) %>% mutate(ratio=count/sum(count)) %>% group_by(class) %>% mutate(mean=sum(ratio)/6)

write.table(x=ratioAndMean, "/Users/gianlucadamaggio/projects/miscarriage/html/html_grep/plot/code/ratioAndMean_variantClass.tsv", sep="\t",row.names=F)
###

vc=ggplot(myd, aes(id, count, fill=class)) + geom_bar(stat='identity')+ scale_fill_npg() + theme_bw() + xlab("ID") + ggtitle("Variant classes") + theme(legend.title=element_blank(),legend.text=element_text(size=15))


ggsave(paste(pfolder,'variantClass.png',sep=''), vc , width = 20, height = 15, units = "cm", dpi = 300)


