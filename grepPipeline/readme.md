
#### 1. format data base files per chromosome 
- merge gene lists: merge_allGenes.py  
Input gene list: tab-delimited two column file with ENSGENEID, category. Genes in more than one category can be repeated 
Output: merged file 
- fathmm 


#### 2.  Annotate variants using `V.E.P`
merged vcf per chromosome. output is a table

all.sample.chrx.vcf -> all.sample.chrx.vep.tsv
all.control.chrx.vcf -> all.control.chrx.vep.tsv  
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/

```
vep --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /data/biocontainers/vepcache --offline --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,IMPACT,SYMBOL,STRAND,SIFT,PolyPhen,EXON,AF,AFR_AF,AMR_AF,ASN_AF,EUR_AF,EAS_AF,SAS_AF,AA_AF,EA_AF,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF,gnomAD_SAS_AF,MAX_AF,CADD_RAW,CADD_PHRED" --force_overwrite --variant_class -i all.sample.chrx.vcf --plugin CADD,/lustre/home/enza/CADD/whole_genome_SNVs.tsv.gz -o  all.sample.chrx.vep.tsv && touch tabOk/all.sample.chrx.table_ok
```

#### 2.1  Remove Header of file VEP
for i in {1..22} X ; do grep -v "##" ${i}.vep.tsv > ${i}.vep.noheader.tsv; done
for i in {1..22} X; do sed -i '/Uploaded_variatio/s/^#//g' ${i}.vep.noheader.tsv ; done

#### 3. Parsing VEP annotations with `grep.csq.sql.py` 
 per chromosome
```
python3 /grepPipeline/scr/grep.csq.sql.py -db /db/sample.chrx.db -i vep/all.sample.chrx.vep.tsv -g testdata/db/all_geneList.tsv -pli testdata/db/pLIscore.tsv -fathmmcoding /lustre/home/enza/fathmm_database/coding/fathmm_xf_coding.$(chr).tsv -fathmmnc /lustre/home/enza/fathmm_database/noncoding/fathmm_xf_noncoding.$(chr).tsv -o /lustre/home/enza/grep/csq/all.sample.chrx.sql.csq
```


#### 4. Extract individual's allele count from vcf 
1. Extract counts with vcftools 
```
vcftools  --gzvcf all.sample.chrx.vcf.gz --out $(chr)/$(id).$(chr)_counts  --counts --indv $(id)
```
2. Change file formatting with `altCounts.py`
```
python3 grepPipeline/scr/altCounts.py -i /$(chr)/$(id).$(chr)_counts.frq.count -id $(id)
```

#### 5. Filter variable sites with `grep.filter.py`
per chromosome 
```
python3 scr/grep.filter.hgdp.py -ccsq testdata/con/test.control.chr22.csq -cvcf  testdata/con/test.control.vcf.gz  -cl testdata/con/test.control.id  -n 4 -i 5  -gt .5  -scsq testdata/sam/test.samples.chr22.csq  -svcf testdata/sam/test.samples.chr22.vcf.gz  -sl testdata/sam/test.samples.id -r False  -pli 0 -g 0 -cadd 0  -o cicci.out -e cicci.err  -ac 1
```

#### 6. Format results with grep.curation.py 



