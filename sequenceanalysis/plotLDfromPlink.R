library(tidyr) 
library(dplyr) 
library(ggplot2) 
library(reshape2) 
library(gplots)

#####1.import data 
path_chr9 <- ('hgdp_wgs.eur.2Mbrs7859844.ld.gz.ld') 

data_9 <- read.table(path_chr9, header = T) %>% select(BP_A, BP_B, R2)

num_kb = 10
distance_b = num_kb * 1000

data_9s <- subset(data_9,(abs(BP_A - 81063077) <= distance_b) | (abs(BP_B - 81063077) <= distance_b))

######2.convert long-to-wide 
data_9_matrix <- dcast(data_9s, BP_A ~ BP_B, value.var = "R2") # convert to matrix with column AND rownames 
myM <- as.matrix(data_9_matrix[ , -1 ]) 
row.names(myM) <- data_9_matrix$BP_A #I am converting all NAs to 0, reconsider if this is suitable in your case. 
#myM[ is.na(myM) ] <- 0

######3.create a mask for loci of interest 
myMask = myM != 'X' # Hack per inizializzare un dataframe uguale a quello contenente i LD 
myMask[myMask] = '' # Hack per mettere ovunque stringa vuota

#Inserimento coppie di interesse: 
#- round ti esprime il valore con 2 cifre decimali 
#- format fa si che 1 sia mostrato come 1.00 (sempre con 2 cifre decimali)

myMask['81062998', '81063077'] = format(round(myM ['81062998', '81063077'], 2), nsmall = 2) 
myMask['81063001', '81063077'] = format(round(myM ['81063001', '81063077'], 2), nsmall = 2)

#####4. plot 
my_palette <- colorRampPalette(c("white", "yellow","orange", "red"))(n = 299)

png("LD.png", width = 20, height = 15, units = "cm", res = 300) 
heatmap.2(myM, key.xlab="LD", cellnote = myMask, notecol="black", Colv = NA, Rowv = NA, scale = "none", col = my_palette, trace = "none", density.info = "none", margins = c(7,14), main = "Linkage Disequilibrium")

dev.off()

#####5.job  [kore-plotLDfromPlink.sh]
singularity exec /mpba0/mpba-sw/biocontainers/r-bioconductor-base2.img Rscript /mpba0/vcolonna/flavia/scriptR/kore-plotLDfromPlink.R


