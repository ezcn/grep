
### 1. create a bed file
paper_miscarriage.bed: Chr-start-end (+1)

### 2. create 2Mb intervals from paper coordinates
```
bedtools slop -i paper_miscarriage.bed -g /mpba0/vcolonna/IMMA/hg38p12/hg38.p12.fa.fai -b 2000000 > intervals_to_check.bed    # 
```
-b: add a fixed number of bases in each direction

### 3. Calculate Linkage Disequilibrium from HGDP database using PLINK using [kore-plinkLD.sh]()

Only on chromosome 9 2Mb surrounding rs7859844  chr9:79063076-83063077 - 102856 variants - 155 individuals of European ancestry 


```
singularity exec /mpba0/mpba-sw/biocontainers/plink.img plink --vcf /mpba0/vcolonna/flavia/WGS/hgdp_wgs.20190516.full.chr9.vcf.gz --r2 --out /mpba0/vcolonna/flavia/ldchr9/hgdp_wgs.eur.2Mbrs7859844.ld.gz --chr 9 --from-bp 79063076 --to-bp 83063077 --keep EUR.list --ld-window-kb 200 --ld-window-r2 0.4
```

### 4. plot LD matrix  using [plotLDfromPlink.R]()

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gplots)

#####1.import data
path_chr9 <- ('hgdp_wgs.eur.2Mbrs7859844.ld.gz.ld')
data_9 <- read.table(path_chr9, header = T) %>% select(BP_A, BP_B, R2)

######2.convert long-to-wide
data_9_matrix <- dcast(data_9, BP_A ~ BP_B, value.var = "R2") # convert to matrix with column AND rownames
myM <- as.matrix(data_9_matrix[ , -1 ])
row.names(myM) <- data_9_matrix$BP_A
#I am converting all NAs to 0, reconsider if this is suitable in your case.
#myM[ is.na(myM) ] <- 0


######3.create a mask for loci of interest 
myMask = myM != 'X' # Hack per inizializzare un dataframe uguale a quello contenente i LD
myMask[myMask] = '' # Hack per mettere ovunque stringa vuota

#Inserimento coppie di interesse:
#- round ti esprime il valore con 2 cifre decimali
#- format fa si che 1 sia mostrato come 1.00 (sempre con 2 cifre decimali)

myMask['81062998', '81063077'] = format(round(myM ['81062998', '81063077'], 2), nsmall = 2)
myMask['81063001', '81063077'] = format(round(myM ['81063001', '81063077'], 2), nsmall = 2)
myMask['81063077', '81063204'] = format(round(myM ['81063077', '81063204'], 2), nsmall = 2)


#####4. plot 
my_palette <- colorRampPalette(c("white", "yellow","orange", "red"))(n = 299)

png("LD.png", width = 20, height = 15, units = "cm", res = 300)
heatmap.2(
 myM,
 key.xlab="LD",
 cellnote = myMask,
 notecol="black",
 Colv = NA, Rowv = NA,
 scale = "none",
 col = my_palette, 
 trace = "none",
 density.info = "none",
 margins = c(7,14),
 main = "Linkage Disequilibrium"
 )

dev.off()







