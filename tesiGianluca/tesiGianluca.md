# 1. Merge VCF (Chr 22) of multiple samples into one VCF 

### 1.1 tabix of chr22 in our samples and filter with QUAL > 20 
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do tabix -h $id.fullvep.vcf.gz chr22 | vcffilter -f "QUAL > 20" > $id.chr22.vep.vcf ; done
```
### 1.2 merge VCF
```
for c in $(seq 1 22 ) ; do qsub -e /mpba0/vcolonna/gianluca/junkfile/mergChr$c.err -o /mpba0/vcolonna/gianluca/junkfile/mergChr$c.out -v chr="chr$c" -N MergeChr$c /mpba0/vcolonna/gianluca/job/kore-bcftoolsMerge.sh ; done
```
# 2. CSQ Allele freq using [CSQfreqAnnotation](../filtering/CSQfreqAnnotation.py)
```
python3 CSQfreqAnnotation.py -f /mpba0/vcolonna/gianluca/TESI/MergedFreqScript/merged.chr22.vep.vcf.gz -o merged.chr22.CSQfreq.vep.vcf -e freq.err
```
### 2.1 bgzip output
```
bgzip merged.chr22.CSQfreq.vep.vcf
```
### 2.2 index of gz
```
tabix -h merged.chr22.CSQfreq.vep.vcf.gz
```
# 3. Use VCF2TSV to convert a VCF into a TSV file
```
vcf2tsv -g merged.chr22.CSQfreq.vep.vcf.gz > merged.chr22.CSQfreq.vep.tsv
```
### 3.1 Import TSV file in R
```
merged_chr22<-read.table("/home/gianluca/project/TESI/merged.chr22.CSQfreq.vep.tsv", sep= "\t", header = T)
```
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



