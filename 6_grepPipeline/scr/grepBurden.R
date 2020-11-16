library(dplyr) 
library(argparse)
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--countfile", action="store")
parser$add_argument("--numberCases", action="store")
parser$add_argument("--numberControls", action="store_true")
args <- parser$parse_args()

# lp95obs --> observed -log10(p-value) @95% of all genes 
# lp95exp --> expected -log10(p-value) @95% of all genes 
# lp0obs --> always zero, observed p-value of the gene with the max expected -log10(p-value) among genes with p-value==1
# lp0exp --> expected p-value of the gene with the max expected -log10(p-value) among genes with p-value==1
# INPUT: table og genes and counts 
# ca1H --> count of cases with at least one hit in the gene 
# co1H --> count of controls with at least one hit in the gene
# ca2H --> count of cases with two or more hits in the gene 
# co2H --> count of controls with two or more hits in the gene 


#totcase=numberCases
#totcon=numberControls
# read count file 
myd=read.table(countfile, header=T )
#pvalues recessive and dominant 
myd<- myd %>%  rowwise() %>% mutate(pdom=fisher.test(matrix(c(ca1H, numberCases-ca1H, co1H, numberControls-co1H), ncol=2, nrow=2 ) )$p.value , prec=fisher.test(matrix(c(ca2H, numberCases-ca2H, co2H ,  numberControls-co2H), ncol=2, nrow=2 ) )$p.value, log10pdom=(-log10(pdom) ), log10prec=(-log10(prec)))


#dominant
#calulate expected pval 
myd<-myd[order(myd$log10pdom),]
myd$log10pdom_exp<-sort(-log10(runif(nrow(myd)))) ## 'nrows' random values from uniform distribution, sorted   

lp95obs <- quantile( subset(myd, log10pdom!=0)$log10pdom, 0.95) #p95
lp95exp <- quantile( subset(myd, log10pdom!=0)$log10pdom_exp, 0.95) #u95
lp0exp <- maxlogp1 <- max(subset(myd, pdom==1)$log10pdom_exp) #u0
lp0obs <-0 #p0

lambda95 <-(lp95obs[[1]]-lp0obs) / (lp95exp[[1]]-lp0exp) #slope
yint<-lp95obs[[1]]-lambda95*lp95exp[[1]]

maxp<-ceiling(max(myd$log10pdom))

myd %>% ggplot (aes(log10pdom_exp, log10pdom) )+ geom_point() + geom_abline(intercept=yint, slope=lambda95, linetype=3, color='grey')+ geom_abline(intercept=0, slope=1, color='grey') + xlim(0, maxp)+ylim(0, maxp) + theme_minimal() +labs(x='expected -log10(p-value)' , y='observed -log10(p-value)' , title= paste('Dominant test - lambda95 =' ,round(lambda95, 3)  , sep= ' ' )  ) 

ggsave('dominant.png') 

















#myd %>% select (gene, pdom, log10pdom, log10pdom_exp)  %>% ggplot(aes( sample=log10pdom) ) +stat_qq() +geom_qq_line(color=1, linetype=3 ) +theme_minimal() +labs(y='-log10(p-valueDominant)')



myd %>% filter (pdom==1 ) %>% max(log10pdom_exp)




#recessive
myd<-myd[order(myd$log10prec),]
myd$log10prec_exp<-sort(-log10(runif(nrow(myd))))
