# Variatn Effect Prediction

### 1) Use a tool (like VEP) on file.VCF for predict variant effect
```
./vep --af --appris --biotype --buffer_size 500 --check_existing --distance 5000 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --input_file [input_data] --output_file [output_file]
```
### 2) Use Bcftools split-vep for expand information inside output file of VEP analysis
```
bcftools +split-vep -f '%CHROM %POS %CSQ\n' -d -A tab OgzFBNb3VDOntEp9.vcf | less -S
```
https://samtools.github.io/bcftools/howtos/plugin.split-vep.html

### 3) Find gene correlate to embryo development 

http://amigo.geneontology.org/amigo/term/GO:0009790#display-lineage-tab

transform UniProtKB AC/ID in Ensembl ID:

https://www.uniprot.org/mapping/

Use BioMart for find Start/End position of Ensembl ID [script](biomartScript/biomaRt.R)


