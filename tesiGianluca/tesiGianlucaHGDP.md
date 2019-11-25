# MAIN GOAL:  Want to obtain the positions that occurs in our Population (GREP) and in the control one (HGDP), after that we have to take only position that are shared between them.

# 1) Prepare 1st bed file from GREP and HGDP with [kore-positionGREP.sh](kore-positionGREP.sh) and [kore-positionHGDP](kore-positionHGDP.sh).
```
for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/gianluca/junkfiles/positionGREP.$c.err -o /mpba0/vcolonna/gianluca/junkfiles/positionGREP.$c.out -v chr="chr$c" -N GREPpos.$c /mpba0/vcolonna/gianluca/job/kore-positionGREP.sh ; done 
```
# 2) Obtain unique position from both bed.file with [kore-positionGREPuniq.sh](kore-positionGREPuniq.sh) and [kore-positionHGDPuniq.sh](kore-positionHGDPuniq.sh)
```
for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/gianluca/junkfiles/positionGREPuniq.$c.err -o /mpba0/vcolonna/gianluca/junkfiles/positionGREPuniq.$c.out -v chr="chr$c" -N GREPposU.$c /mpba0/vcolonna/gianluca/job/kore-positionGREPuniq.sh ; done 
```
### 2.1) Create one file with both Position.
```
cat postionGREP.chr22.bed > positionALLuniq.bed
```
```
cat positionHGDP.chr22.bed >> positionALLuniq.bed
```
### 2.3) Create a file with only shared position taken one time.
```
cat positionALLuniq.bed | sort | uniq -d > positionALLuniqB.chr22.bed
```
# 3) Recode VCF with (positionALLuniqB.chr22.bed) using [kore-recodeGREP.sh](kore-recodeGREP.sh) and [kore-recodeHGDP.sh](kore-recodeHGDP.sh).
```
for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/gianluca/junkfiles/recodeGREP.$c.err -o /mpba0/vcolonna/gianluca/junkfiles/recodeGREP.$c.out -v chr="chr$c" -N GREPrecode.$c /mpba0/vcolonna/gianluca/job/kore-recode.sh 
```
# 4) gzip file and use Python script for [GREP](../filtering/AFS-GREP_grepl.py) and [HGDP](../filtering/AFS-HGDP_random_grepl.py) for obtain frequency of alleles in VCF with [kore-grepPy.sh](kore-grepPy.sh) and [kore-hgdpPy.sh](kore-hgdpPy.sh).
```
for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/gianluca/junkfiles/pyGREP.$c.err -o /mpba0/vcolonna/gianluca/junkfiles/pyGREP.$c.out -N GREPpy.$c /mpba0/vcolonna/gianluca/job/kore-grepPy.sh 
```
# 5) Plot using [hgdp_grep_plot.R](hgdp_grep_plot.R) with [kore-hgdp_grep_plot.sh](kore-hgdp_grep_plot.sh).
```
qsub -e /mpba0/vcolonna/gianluca/junkfiles/plotR.err -o /mpba0/vcolonna/gianluca/junkfiles/plotR.out -N plotR /mpba0/vcolonna/gianluca/job/kore-hgdp_grep_plot.sh
```

### 5.1) CSQfreq wrap for Consequence and  fill by VariantClass [plot.png](plotWrapConsequence.png)

### 5.2) CSQfreq > 0 wrap for Consequence and  fill by VariantClass [plot.png](plotWrapConsequenceMTZ.png)

### 5.3) CSQfreq wrap for CSQrank and  fill by VariantClass [plot.png](plotWrapCSQrank.png)

### 5.4) CSQfreq > 0 wrap for CSQrank and  fill by VariantClass [plot.png](plotWrapCSQrankMTZ.png)


# 6) Instead of obtain a big file I have made another script in python [AFS-grepMean.py](AFS-grepMean.py) for our Population and [AFS-hgdpMean.py](AFS-hgdpMean.py) for the Control one.

### 6.1) [This](hgdp_grep_mean.tsv) is an Output file obtained for each Chromosome of Control-Population (HGDP) and in each row there is a "MEAN" for each Variant_class and Consequence of this Chromosome , repeted "n" time for "6" (our population is 6 samples) randomized samples in HGDP population. In the end we have also one line for each Chromosome for my GREP population.

# 7) [this](hgdp_grep_mean.tsv) file have to be processed with R for visualize an enrichment in specific Variant Class or Consequences.
