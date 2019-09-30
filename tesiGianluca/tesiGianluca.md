# 1. Merge VCF of multiple samples, by "chr", into one VCF 

```
for c in $(seq 1 22 ) ; do qsub -e /mpba0/vcolonna/gianluca/junkfile/mergChr$c.err -o /mpba0/vcolonna/gianluca/junkfile/mergChr$c.out -v chr="chr$c" -N MergeChr$c /mpba0/vcolonna/gianluca/job/kore-bcftoolsMerge.sh ; done
```
# 2. Prepare files for AFS analisys [AFS-script](../filtering/AFS_grepl.py)
```
for c in $(seq 1 22 ) ; do qsub -e /mpba0/vcolonna/gianluca/junkfile/AFS_Chr$c.err -o /mpba0/vcolonna/gianluca/junkfile/AFS_Chr$c.out -v chr="chr$c",output="-o /mpba0/vcolonna/gianluca/TESI/AFS/ScriptProcessed/merged.chr$c.AFS.fb.vep.vcf",error="-e /mpba0/vcolonna/gianluca/junkfile/AFSchr$c.err",csqimpact="-v /mpba0/vcolonna/gianluca/TESI/AFS/csqimpact.tsv" -N AFS_Chr$c /mpba0/vcolonna/gianluca/job/kore-scriptAFSpython.sh ; done
```
### 2.1 bgzip output
```
bgzip merged.chr$c.AFS.fb.vep.vcf
```
### 2.2 index of gz
```
tabix -h merged.chr$c.AFS.fb.vep.vcf.gz
```
# 3. Use VCF2TSV to convert a VCF into a TSV file
```
for c in $(seq 1 22 ) ; do qsub -e /mpba0/vcolonna/gianluca/junkfile/vcfChr$c.err -o /mpba0/vcolonna/gianluca/junkfile/vcfChr$c.out -v chr="chr$c" -N vcfChr$c /mpba0/vcolonna/gianluca/job/kore-vcf2tsv.sh ; done

```
### 3.1 Delete '#' at beginning 
```
for c in $(seq 1 22 ) ; do qsub -e /mpba0/vcolonna/gianluca/junkfile/sedChr$c.err -o /mpba0/vcolonna/gianluca/junkfile/sedChr$c.out -v chr="chr$c" -N sedChr$c /mpba0/vcolonna/gianluca/job/kore-sedRemove.sh ; done
```
### 4. Processing TSV file in R [plotAFS](plotAFS.R) 
```
qsub -e /mpba0/vcolonna/gianluca/junkfile/rplot.err -o /mpba0/vcolonna/gianluca/junkfile/rplot.out -N Rplots /mpba0/vcolonna/gianluca/TESI/AFS/kore-RscriptPlotAFS.sh
```






#                                       DA QUI IN POI : WORK IN PROGRESS  



# 4. Use splitVEP (BCFTOOLS plugin) for obtain VEP annotation as TSV file 

### 4.1 header-splitvep is obtained from :
```
bcftools +split-vep -f '%CHROM %POS %CSQ\n' -l -d -A tab merged.chr22.CSQfreq.vep.vcf.gz > header-splitvep
```
### 4.2 append to header-splitvep information from VEP annotation and rename
```
bcftools +split-vep -f '%CHROM %POS %CSQ\n'  -d -A tab merged.chr22.CSQfreq.vep.vcf.gz >> header-splitvep
```
```
mv header-splitvep merged.chr22.CSQfreq.splitvep.tsv
```
### 4.1 Import splitVEP.vcf in R
```
splitvep<-read.table("/home/gianluca/project/TESI/merged.chr22.CSQfreq.splitvep.tsv", sep = "\t", header = T)
```
# 5. Bind Merged.vcf and splitVEP.vcf
```
binding<-merge(merged_chr22,splitvep, by = "POS")
```



