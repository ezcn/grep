### 1. Create a BED file with genomic position of Pseudogene and CDS regions.



### 2. Merge VCF (Chr 22) of multiple samples into one VCF 
'''for id in AS054 AS064 AS074 AS090 AS094; do tabix -h $id.fullvep.vcf.gz chr22 | vcffilter -f "QUAL > 20" > $id.chr22.vep.vcf ; done'''


### 3. Annotation of our samples with BED.file 


### 4. Evaluate the freq of allele with CSQ in VCF_merged in our population


### 5. Use a tool to convert a VCF into a TSV.file for import in R
