
# 1. create a bed file
paper_miscarriage.bed: Chr-start-end (+1)

# 2. create 2Mb intervals from paper coordinates
```
bedtools slop -i paper_miscarriage.bed -g /mpba0/vcolonna/IMMA/hg38p12/hg38.p12.fa.fai -b 2000000 > intervals_to_check.bed    # 
```
-b: add a fixed number of bases in each direction

# 3. Calculate Linkage Disequilibrium from HGDP database using PLINK using kore-plinkLD.sh

Only on chromosome 9 

```
singularity exec /mpba0/mpba-sw/biocontainers/plink.img plink --vcf hgdp_wgs.20190516.full.chr9.vcf.gz --r2 --out chr9_ld_snp --chr 9 --from-bp 42412110 --to-bp 44412110 
```

plink --bfile (chr21) --ld rs183453668  #potrei vedere se vi sono regioni di LD intorno a tale variante, causativa miscarrage (improbabile)
plink --bfile chr22_1000Gphase3_EUR_snp_maf_rmvsnp --r2 --out chr22_1000Gphase3_EUR_ldtable --ld-window-kb 1000 #LD in quel cromosoma

--r calculates and reports raw inter-variant allele count correlations (reports all results in table format) 

--r2 reports squared correlations (reports all results with a text matrix)

--ld-window-kb 10000000 (diecimila kb)

#Merge cromosomi: 
#~/bin/bcftools concat -o /mpba0/vcolonna/IMMA/samples/gatkConcat/${id}.chr21.fb.vep.vcf.gz 

#for id in AS006 AS054 AS064 AS074 AS090 AS094 ;do echo qsub -e /mpba0/vcolonna/flavia/$id.vep.err -o /mpba0/vcolonna/flavia/$id.vep.out -v #id="$id" -N vepmerge$id /mpba0/vcolonna/flavia/job/kore-vep.sh; done

#/mpba0/vcolonna/bin/bcftools merge /mpba0/vcolonna/IMMA/samples/fb/vep/${id}.chr21.fb.vep.vcf.gz -O z -o /mpba0/vcolonna/flavia/WGS/
#${id}.chr21.fb.vep.merge.vcf.gz


GRAFICO
library(tidyr)
library(dplyr)
path_chr21 <- ('/home/flavia/Desktop/GWAS VARIANTS /chr21__EUR_ldtable.ld')

mysnp=43412110
halfinterval=5000000
 	
data <- read.table(path_chr21, header = T) %>% select(BP_A, BP_B, R2)

#data <- read.table(path_chr21, header = T)
#data.filter = data %>% filter(BP_A>mysnp-halfinterval & BP_A<mysnp-halfinterval & BP_B>mysnp-halfinterval & BP_B<mysnp-halfinterval) 

subData <- data[1:1000,] #metto le prime 100 righe e tutte le colonne tutte
library(reshape2)
#convert long-to-wide

x <- dcast(subData, BP_A ~ BP_B, value.var = "R2")
# convert to matrix with column AND rownames

myM <- as.matrix(x[ , -1 ])
row.names(myM) <- x$BP_A

# I am converting all NAs to 0, reconsider if this is suitable in your case.

myM[ is.na(myM) ] <- 0
library(gplots)

heatmap.2(myM, Colv = NA, Rowv = NA, scale = "none", col = bluered(100), trace = "none", density.info = "none")

my_palette <- colorRampPalette(c("white", "yellow","orange", "red"))(n = 299)
heatmap.2(myM, Colv = NA, Rowv = NA, scale = "none", col = my_palette, trace = "none", density.info = "none")

ggsave("/mpba0/vcolonna/flavia/CDrate.png", plot= myplot, device="png", width = 20, height = 15, units = "cm", dpi = 300)


2)singularity exec /mpba0/mpba-sw/biocontainers/plink.img plink --vcf hgdp_wgs.20190516.full.chr21.vcf.gz --r2 --out chr21__threshold_ldtable --ld-window-kb 10000000 --ld-window-r2 0.5  #soglia R2
 
Variante chr21: rs183453668 #10 MB prima e dopo dalla variante



a)Importarlo in R 

library(tidyr)
library(dplyr)
path_chr21ld = ('/home/flavia/Desktop/GWAS_VARIANTS/chr21_ld_snp.ld')
dataLD <- read.table(path_chr21ld, header = T) %>% select(BP_A, BP_B, R2)
subDataLD <- dataLD[1:100,]
library(reshape2)
xLD <- dcast(subDataLD, BP_A ~ BP_B, value.var = "R2")
myMLD <- as.matrix(xLD[ , -1 ])
row.names(myMLD) <- xLD$BP_A
myMLD[ is.na(myMLD) ] <- 0
library(gplots)
my_palette <- colorRampPalette(c("white", "yellow","orange", "red"))(n = 299)
heatmap.2(myMLD, Colv = NA, Rowv = NA, scale = "none", col = my_palette, trace = "none", density.info = "none")








2)extract the variants which fall in the intervals linkage disequilibrium

tabix /mpba0/vcolonna/IMMA/samples/fb/vep/AS054.chr21.fb.vep.vcf.gz  -B intervals_to_check.bed (regioni in linkage disequilibrium) > chr21.fb.subset.vcf


Fenomeno per cui a livello di popolazione specifiche combinazioni di alleli a due o più loci concatenati tendono a trovarsi insieme sullo stesso cromosoma piu' frequentemente di quanto ci si attenda per caso.

