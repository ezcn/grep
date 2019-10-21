Want to obtain the positions that occurs in our Population (GREP) and in the control one (HGDP), after that we have to take only 
position that are shared between them. For do that :

### 1) Prepare 1st bed file from GREP and HGDP with [kore-positionGREP.sh](kore-positionGREP.sh) and [kore-positionHGDP](kore-positionHGDP.sh)
```
qsub -e /mpba0/vcolonna/gianluca/junkfiles/positionGREP.err -o /mpba0/vcolonna/gianluca/junkfiles/positionGREP.out -N GREPpos /mpba0/vcolonna/gianluca/job/kore-positionGREP.sh 
```
### 2) After that have to obtain unique position from both bed.file :
```
cat positionGREP.chr22.bed | sort | uniq -u > positionGREPuniq.chr22.bed
```
```
cat positionHGDP.chr22.bed | sort | uniq -u > positionHGDPuniq.chr22.bed
```
# 2.1) From this two files have to create one file with both Position:
```
cat postionGREP.chr22.bed > positionALLuniq.bed
```
```
cat positionHGDP.chr22.bed >> positionALLuniq.bed
```
# 2.3) Now have to create a file with only shared position take one time :
```
cat positionALLuniq.bed | sort | uniq -d > positionALLuniqB.chr22.bed
```
### 3) Recode VCF with (positionALLuniqB.chr22.bed) using 
