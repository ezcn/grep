library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggsci)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify two arguments: <inputFile> <uniqueCodeForResults>', call.=FALSE)
} 

code=args[2]
myd=read.table(args[1], header=T, sep='\t')
myd$impact= factor(myd$impact,levels=c('LOW', 'MODERATE','HIGH'))

#### summary stats 
sink(paste(code, '.output.txt', sep=''))
totnumberofuniqevariants=myd %>% select(index_x) %>% distinct() %>% tally() 
print(paste(code, 'Tot Number of Uniqe Variants', totnumberofuniqevariants, sep='\t'))
totalnumberuniqueNOB=myd %>% filter(rare=='NOB') %>% select (index_x) %>% distinct() %>% tally()
print(paste(code, 'Tot Number of NOB',totalnumberuniqueNOB, sep='\t'))

impdf<-myd %>% select (index_x , impact ) %>% distinct() %>% group_by(impact) %>% tally()
print(impdf)


avsdSites<- myd %>% select(index_x, sample) %>% distinct() %>% group_by(sample) %>% tally() %>% summarize(mean(n), sd(n))
print(paste(code, 'Average number of sites per sample', avsdSites[1], sep='\t'))
print(paste(code, 'Sd number of sites per sample', avsdSites[2], sep='\t'))
totnumberofuniqegenes=myd %>% select(gene_id) %>% distinct() %>% tally()
print(paste(code, 'Total number of Unique Genes', totnumberofuniqegenes, sep='\t'))
totnumberofuniqetranscripts=myd %>% select(element_id) %>% distinct() %>% tally()
print(paste(code, 'Total number of Unique Transcripts', totnumberofuniqetranscripts, sep ='\t'))
avesdGenes<- myd %>% select(gene_id, sample) %>% distinct() %>% group_by(sample) %>% tally() %>% summarize(mean(n), sd(n))
print(paste(code, 'Average number of genes per sample', avesdGenes[1], sep='\t'))
print(paste(code, 'Sd number of genes per sample', avesdGenes[2], sep='\t'))
avesdTranscripts<- myd %>% select(element_id, sample) %>% distinct() %>% group_by(sample) %>% tally() %>% summarize(mean(n), sd(n))
print(paste(code, 'Average number of transcripts per sample', avesdTranscripts[1], sep='\t'))
print(paste(code, 'Sd number of transcripts per sample', avesdTranscripts[2], sep='\t'))
sink()


####### 1. remove genes with too many variants 
### Multiple variants in the same gene

print('multiple variants same gene')
myd %>% select(index_x,sample, gene_symbol, impact) %>% distinct() %>% group_by( sample, gene_symbol, impact) %>% tally() %>% filter(n>1) %>% spread(impact, n) %>% write.table(paste(code, '.multiplevariants_pergene.tsv'), sep="\t", quote=F, row.names=F, col.names=T)

mycol=c("#FACF5A", "#49BEB7", "#FF5959")
myd %>% select(index_x,sample, gene_symbol, impact) %>% distinct() %>% group_by( sample, gene_symbol, impact) %>% tally() %>% filter(n>1) %>% ggplot(aes(reorder(gene_symbol, n) , n , fill=impact ))  +  geom_bar(stat='identity')+coord_flip() +theme_minimal() +scale_fill_manual(values=mycol)
ggsave('variantspergenemore2.png')

#myd %>% select(index_x,sample, gene_symbol, impact) %>% distinct() %>% group_by( gene_symbol, impact) %>% tally() %>% spread(impact, n) %>% summarize(tot=sum(LOW, MODERATE, HIGH, na.rm=T )) %>% ggplot(aes(tot)) +geom_density()+theme_minimal()+xlab('Number of variants per gene') + geom_vline( xintercept=quantile(tot, .95)) 


myg =myd %>% select(index_x,sample, gene_symbol, impact) %>% distinct() %>% group_by( gene_symbol, impact) %>% tally() %>% spread(impact, n) %>% summarize(tot=sum(LOW, MODERATE, HIGH, na.rm=T )) 

dens <- density(myg$tot)
df <- data.frame(x=dens$x, y=dens$y)
probs <- c(0,  0.95, 0.99, 1)
quantiles <- quantile(myg$tot, prob=probs)
df$quant <- factor(findInterval(df$x,quantiles))
ggplot(df, aes(x,y)) + geom_line() + geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) + scale_x_continuous(breaks=quantiles) + scale_fill_brewer(palette="Reds", guide="none") +theme_minimal()+labs(x='Number of variants per gene' , y='Density')
ggsave('variantspergenedensity.png')



####### 3. Results 
#######
print('~~ number')
#myd %>% select(sample, index_x, impact ) %>% distinct() %>% group_by (sample, impact) %>% tally() %>% spread(impact, n)
mycol=c("#FACF5A", "#49BEB7", "#FF5959")
number=myd %>% select(sample, index_x, impact ) %>% distinct() %>% group_by (sample, impact) %>% tally() 
#number$impact= factor(number$impact,levels=c('HIGH','MODERATE','LOW'))
ggplot(number, aes(as.factor(sample), n, fill=impact)) + geom_bar(stat='identity') + ggtitle(paste (code,'- Unique variants per sample', sep=' ')) + theme_minimal() + scale_fill_manual(values=mycol)+coord_flip() +xlab('') + ylab('Number of unique variants') 
ggsave(paste(code, '_number.png', sep=''))

######################
print('sumgene')
myd %>% select(sample, gene_id, sumGene)%>% distinct () %>% group_by(sample, sumGene)  %>% tally() %>% ggplot (aes(as.factor(sample), n, fill=sumGene)) + geom_bar(stat='identity')+ theme_minimal()  +coord_flip() +scale_fill_gradient(low = "#132B43",  high = "#56B1F7") +labs(fill = "Number of lists", y='Number of unique genes', x='' )  #+ ggtitle(paste(code,'- Genes in lists', sep=' '))             
ggsave(paste (code, '_sumGene.png', sep=''))  


#####################
print('most_severe_consequence')
csqtype=myd %>% select(sample, index_x, impact, most_severe_consequence) %>% distinct() %>% group_by (sample, impact, most_severe_consequence) %>% tally() 
#number$impact= factor(number$impact,levels=c('HIGH','MODERATE','LOW'))
nbcolors=length(levels(myd$most_severe_consequence)) 
mycolors <- colorRampPalette(brewer.pal(8,'Set2'))(nbcolors)
ggplot(csqtype, aes(as.factor(sample), n, fill=most_severe_consequence)) + geom_bar(stat='identity')  + theme_bw()+ facet_wrap (impact ~ ., scales='free')+ scale_fill_manual(values = mycolors) +coord_flip() +labs(fill = "Most severe Consequence", y='Number of unique variants', x='' ) #+ ggtitle(paste(code, '- most_severe_consequence', sep=' '))
ggsave(paste(code, '_most_severe_consequence.png', sep=''))



######
print('geneTranscripts')
genes =myd %>% select(sample, gene_id ) %>% group_by(sample) %>% distinct() %>% tally()
transcripts = myd %>% select(sample, element_id ) %>% group_by(sample) %>% distinct() %>% tally()
geneTranscripts=full_join(genes, transcripts, by='sample')
colnames(geneTranscripts)=c('sample','genes','transcripts')
geneTranscripts %>% gather ('type','n', c(genes, transcripts)) %>% ggplot(aes(as.factor(sample) , n , color=type) ) + geom_point() + ggtitle(paste (code, '- Genes and transcripts', sep=' ')) + theme_bw() + scale_color_nejm()
ggsave(paste (code, '_geneTranscripts.png', sep='')) 


########################
print('cadd pli ')
myd %>% ggplot(aes(CADD, color=impact)) + geom_density() + ggtitle(paste(code,'- CADD', sep=' ')) + theme_bw() + scale_color_nejm()
ggsave(paste(code, '_CADD.png', sep=''))

myd %>% ggplot(aes(pLI, color=impact)) + geom_density() + ggtitle(paste(code,'- pLI', sep=' ')) + theme_bw() + scale_color_nejm()
ggsave(paste(code, '_pLI.png', sep=''))


############ 
print('plicadsumgene')
myd %>% ggplot(aes(pLI, CADD, color=as.factor(sumGene))) + geom_point() + ggtitle(paste(code, '- pLI CADD sumGene', sep=' ')) + theme_bw() + scale_color_futurama()
ggsave(paste(code, '_pLiCADDsumGene.png', sep='')) 


##########################
print('plicadd')
#mycol=c('#ffbd69', '#fe346e', '#b21f66', '#381460') 
mycol=c('#FFBD69', '#2EF276','#FE346E', '#B21F66', '#381460')
datawithpliandcadd=myd %>% select (pLI, CADD, sumGene, impact) %>% distinct()%>% filter (!is.na(pLI) & !is.na(CADD)) %>% tally() 
myd %>% select (sample, pLI, CADD, sumGene, impact) %>% distinct() %>% ggplot(aes(pLI, CADD, color=as.factor(sumGene), alpha =0.7 )) + geom_point() + ggtitle(paste (code, 'pLI, CADD, gene lists \n# unique sites with data =',datawithpliandcadd, sep='')) + theme_bw() + scale_color_manual (values=mycol) + facet_wrap(sample ~ .)
ggsave(paste(code, '_plicadd.png', sep='')) 



### Aggregate analysis 
print('aggregate')
mycol=c("#FACF5A", "#49BEB7", "#FF5959")
myd %>% select (impact, gene_symbol, sample) %>%distinct() %>% group_by( impact, gene_symbol)  %>% tally() %>% filter(n>1) %>% ggplot( aes(reorder(gene_symbol, n), n, fill=impact)) + geom_bar(stat='identity', position='dodge') + xlab('genes') + ylab('number of samples') + scale_fill_manual(values=mycol) + coord_flip() + theme_minimal() + facet_wrap (impact ~ ., scales ='free') + ggtitle(paste(code, '- Genes shared among samples', sep=' ')) 
ggsave(paste(code, '_aggregate.png', sep='')) 


myd %>% select (index_x, impact, gene_symbol, sample) %>%distinct() %>% group_by( impact, gene_symbol)  %>% tally() %>% filter(n>1) %>% ggplot( aes(reorder(gene_symbol, n), n, fill=impact)) + geom_bar(stat='identity', position='dodge') + xlab('genes') + ylab('number of samples') + scale_fill_manual(values=mycol) + coord_flip() + theme_minimal() + facet_wrap (impact ~ ., scales ='free') + ggtitle(paste(code, '- Genes shared among samples', sep=' ')) 
ggsave(paste(code, '_aggregateindex_X.png', sep='')) 

