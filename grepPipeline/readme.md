
1. merge gene lists: merge_allGenes.py  

Input gene list: tab-delimited two column file with ENSGENEID, category. Genes in more than one category can be repeated 
Output: merged file 

2. run vep produce json 


3. annotate vcf:   grep.csq.py

python3 scr/grep.csq.py  -j testdata/AS006.chr22.vep.json -g testdata/all_geneList.tsv  -pli testdata/pLIscore.tsv  -cadd testdata/CADD.chr22.tsv.gz  -fathmm testdata/fathmmxf.chr22.tsv.gz  -v testdata/csqimpact.tsv  -c 1 -r 0.05 -o test.csq.out -e test.csq.err -fathmmcoding testdata/fathmmxf.chr22.tsv.gz

4. grep.filter.py
5. grep.curation.py 



