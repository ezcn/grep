
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
singularity exec /mpba0/mpba-sw/biocontainers/plink.img plink --vcf /mpba0/vcolonna/flavia/WGS/hgdp_wgs.20190516.full.chr9.vcf.gz --r2 --out /mpba0/vcolonna/flavia/ldchr9/hgdp_wgs.eur.2Mbrs7859844.ld.gz --chr 9 --from-bp 79063076 --to-bp 83063077 --keep EUR.list --ld-window-kb 200 
```

### 4. plot LD matrix using [plotLDfromPlink.R]()







