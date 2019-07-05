# Variatn Effect Prediction

### 1) Use [VEP](kore-vep.sh) for predict variant effect 
```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ;do echo qsub -e /mpba0/vcolonna/gianluca/$id.vep.err -o /mpba0/vcolonna/gianluca/$id.vep.out -v id="$id" -N vep$id /mpba0/vcolonna/gianluca/kore-vep.sh; done
```
### 2) Use Bcftools split-vep for expand information inside output file of VEP analysis
```
bcftools +split-vep -f '%CHROM %POS %CSQ\n' -d -A tab [filename].vcf | less -S
```
https://samtools.github.io/bcftools/howtos/plugin.split-vep.html

### 3) Genes correlate to embryo development 

3.1) Find genes: http://amigo.geneontology.org/amigo/term/GO:0009790#display-lineage-tab  

3.2) Transform UniProtKB AC/ID in Ensembl ID: https://www.uniprot.org/mapping/

3.3) Use BioMart to find Start/End position of Ensembl ID [(Script)](biomartScript/biomaRt.R) , the output is a BED file.


### 4) Intersect the output of VEP analysis with embryo development bed file using [Bedtools](kore-bedintersect.sh)

```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ;do qsub -e /mpba0/vcolonna/gianluca/$id.int.err -o /mpba0/vcolonna/gianluca/$id.int.out -v id="$id" -N $id.int kore-bedintersect.sh; done
```
