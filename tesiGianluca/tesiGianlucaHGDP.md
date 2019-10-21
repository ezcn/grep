Want to obtain the positions that occurs in our Population (GREP) and in the control one (HGDP), after that we have to take only 
position that are shared between them. For do that :

# 1) Prepare 1st bed file from GREP and HGDP with [kore-positionGREP.sh](kore-positionGREP.sh) and [kore-positionHGDP](kore-positionHGDP.sh)
```
qsub -e /mpba0/vcolonna/gianluca/junkfiles/positionGREP.err -o /mpba0/vcolonna/gianluca/junkfiles/positionGREP.out -N GREPpos /mpba0/vcolonna/gianluca/job/kore-positionGREP.sh 
```
# 2) Obtain unique position from both bed.file :
```
cat positionGREP.chr22.bed | sort | uniq -u > positionGREPuniq.chr22.bed
```
```
cat positionHGDP.chr22.bed | sort | uniq -u > positionHGDPuniq.chr22.bed
```
### 2.1) Create one file with both Position:
```
cat postionGREP.chr22.bed > positionALLuniq.bed
```
```
cat positionHGDP.chr22.bed >> positionALLuniq.bed
```
### 2.3) Create a file with only shared position taken one time :
```
cat positionALLuniq.bed | sort | uniq -d > positionALLuniqB.chr22.bed
```
# 3) Recode VCF with (positionALLuniqB.chr22.bed) using [kore-recodeGREP.sh](kore-recodeGREP.sh) and [kore-recodeHGDP.sh](kore-recodeHGDP.sh)
```
qsub -e /mpba0/vcolonna/gianluca/junkfiles/recodeGREP.err -o /mpba0/vcolonna/gianluca/junkfiles/recodeGREP.out -N GREPrecode /mpba0/vcolonna/gianluca/job/kore-recode.sh 
```
# 4) gzip file and use Python script for [GREP](filtering/AFS-GREP_grepl.py) for obtain frequency of allele in VCF [kore-grepPy.sh](kore-grepPy.sh) and [kore-hgdpPy.sh](kore-hgdpPy.sh)
```
qsub -e /mpba0/vcolonna/gianluca/junkfiles/pyGREP.err -o /mpba0/vcolonna/gianluca/junkfiles/pyGREP.out -N GREPpy /mpba0/vcolonna/gianluca/job/kore-grepPy.sh 
```
# 5) Use [hgdp_grep_plot.R](hgdp_grep_plot.R) with [kore-hgdp_grep_plot.sh](kore-hgdp_grep_plot.sh) for plot :

### 5.1) CSQfreq wrap for Consequence and  fill by VariantClass [plot.png](plotWrapConsequence.png)

### 5.2) CSQfreq > 0 wrap for Consequence and  fill by VariantClass [plot.png](plotWrapConsequenceMTZ.png)

### 5.3) CSQfreq wrap for CSQrank and  fill by VariantClass [plot.png](plotWrapCSQrank.png)

### 5.4) CSQfreq > 0 wrap for CSQrank and  fill by VariantClass [plot.png](plotWrapCSQrankMTZ.png)
