




1. format data base files per chromosome 
- merge gene lists: merge_allGenes.py  
Input gene list: tab-delimited two column file with ENSGENEID, category. Genes in more than one category can be repeated 
Output: merged file 
- cadd 
- fathmm 


2. per chromosome, per merged vcf file run vep produce json 
all.sample.chrx.vcf -> all.sample.chrx.vep.json
all.control.chrx.vcf -> all.sample.chrx.vep.json  

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/


3. per chromosome annotate sites with grep.csq.py 
annotate samples 
python3 scr/grep.csq.hgdp.py  -j testdata/all.sample.chrx.vep.json -g testdata/all_geneList.tsv  -pli testdata/pLIscore.tsv  -cadd testdata/CADD.chr22.tsv.gz  -fathmm testdata/fathmmxf.chr22.tsv.gz  -v testdata/csqimpact.tsv  -c 1 -r 0.05 -o test.all.sample.chrx.csq.out -e test.all.sample.chrx.csq.err -fathmmcoding testdata/fathmmxf.chr22.tsv.gz

annotate controls 
python3 scr/grep.csq.hgdp.py  -j testdata/all.control.chrx.vep.json -g testdata/all_geneList.tsv  -pli testdata/pLIscore.tsv  -cadd testdata/CADD.chr22.tsv.gz  -fathmm testdata/fathmmxf.chr22.tsv.gz  -v testdata/csqimpact.tsv  -c 1 -r 0.05 -o test.all.control.chrx.csq.out -e test.all.control.chrx.csq.err -fathmmcoding testdata/fathmmxf.chr22.tsv.gz

4. per chromosome, per individual filter with grep.filter.py


5. grep.curation.py 



