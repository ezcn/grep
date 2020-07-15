
#### 1. format data base files per chromosome 
- merge gene lists: merge_allGenes.py  
Input gene list: tab-delimited two column file with ENSGENEID, category. Genes in more than one category can be repeated 
Output: merged file 
- cadd 
- fathmm 


#### 2. per chromosome, per merged vcf file run vep produce table
all.sample.chrx.vcf -> all.sample.chrx.vep.tsv
all.control.chrx.vcf -> all.control.chrx.vep.tsv  

```
vep --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /data/biocontainers/vepcache --offline --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,AF,AFR_AF,AMR_AF,ASN_AF,EUR_AF,EAS_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,CADD_RAW,CADD_PHRED" --force_overwrite --variant_class -i all.sample.chrx.vcf --plugin CADD,/lustre/home/enza/CADD/whole_genome_SNVs.tsv.gz -o  all.sample.chrx.vep.tsv && touch tabOk/all.sample.chrx.table_ok
```

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/


#### 3. per chromosome annotate sites with grep.csq.sql.py 
Annotate samples 
```
python3 /grepPipeline/scr/grep.csq.sql.py -db /db/sample.chrx.db -i vep/all.sample.chrx.vep.tsv -g testdata/db/all_geneList.tsv -pli testdata/db/pLIscore.tsv -fathmmcoding /lustre/home/enza/fathmm_database/coding/fathmm_xf_coding.$(chr).tsv -fathmmnc /lustre/home/enza/fathmm_database/noncoding/fathmm_xf_noncoding.$(chr).tsv -o /lustre/home/enza/grep/csq/all.sample.chrx.sql.csq
```

Annotate controls
```
python3 /grepPipeline/scr/grep.csq.sql.py -db /db/control.chrx.db -i vep/all.control.chrx.vep.tsv -g testdata/db/all_geneList.tsv -pli testdata/db/pLIscore.tsv -fathmmcoding /lustre/home/enza/fathmm_database/coding/fathmm_xf_coding.$(chr).tsv -fathmmnc /lustre/home/enza/fathmm_database/noncoding/fathmm_xf_noncoding.$(chr).tsv -o /lustre/home/enza/grep/csq/all.control.chrx.sql.csq
```
#### 4. extract allele count from vcf for each individual and change file formatting
Extract counts with vcftools
```
vcftools  --gzvcf all.sample.chrx.vcf.gz --out $(chr)/$(id).$(chr)_counts  --counts --indv $(id)
```
Change file formatting with python script
```
python3 grepPipeline/scr/altCounts.py -i /$(chr)/$(id).$(chr)_counts.frq.count -id $(id)
```

4. per chromosome, per individual filter with grep.filter.py
python3 scr/grep.filter.hgdp.py -ccsq testdata/con/test.control.chr22.csq -cvcf  testdata/con/test.control.vcf.gz  -cl testdata/con/test.control.id  -n 4 -i 5  -gt .5  -scsq testdata/sam/test.samples.chr22.csq  -svcf testdata/sam/test.samples.chr22.vcf.gz  -sl testdata/sam/test.samples.id -r False  -pli 0 -g 0 -cadd 0  -o cicci.out -e cicci.err  -ac 1

  
5. grep.curation.py 



