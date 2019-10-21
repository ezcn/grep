Want to obtain the positions that occurs in our Population (GREP) and in the control one (HGDP), after that we have to take only 
position that are shared between them. For do that :

### 1) Prepare 1st bed file from GREP and HGDP with [kore-positionGREP.sh](kore-positionGREP.sh)

zcat /mpba0/vcolonna/IMMA/samples/fb/vep/mergedByChr/merged.chr22.fb.vep.vcf.gz | grep -v '#' | cut -f1,2 | awk '{print $1,$2,$2}'| tr " " "\t" > /mpba0/vcolonna/gianluca/TESI/hgdp/positionGREP.${chr}.bed
