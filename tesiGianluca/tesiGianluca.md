# 1. Merge VCF (Chr 22) of multiple samples into one VCF 

### 1.1 tabix of chr22 in our samples and filter with QUAL > 20 
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do tabix -h $id.fullvep.vcf.gz chr22 | vcffilter -f "QUAL > 20" > $id.chr22.vep.vcf ; done
```
### 1.2 merge VCF
```
qsub -e /mpba0/vcolonna/gianluca/TESI/merg.err -o /mpba0/vcolonna/gianluca/TESI/merg.out -N MergeVCF /mpba0/vcolonna/gianluca/job/kore-bcftoolsMerge.sh
```
### 1.3 CSQ Allele freq using [CSQfreq]( grepl/filtering/CSQfreqAnnotation.py )
```
python3 CSQfreqAnnotation.py -f /mpba0/vcolonna/gianluca/TESI/MergedFreqScript/merged.chr22.vep.vcf.gz -o mergedWithFreq.vcf -e freq.err
```

# 2. Evaluate the freq of allele with CSQ in VCF_merged in our population


# 3. Use a tool to convert a VCF into a TSV.file for import in R
