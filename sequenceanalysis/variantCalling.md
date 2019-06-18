
#### download hg38 

#### index the reference


#### align reads to reference genome [kore-align](kore-align.sh)
```
 qsub -o /mpba0/vcolonna/silvia/AS074.out -e /mpba0/vcolonna/silvia/AS074.err -v gz1="AS074-Av-L_R1.fastq.gz",gz2="AS074-Av-L_R2.fastq.gz",idout="AS074",dir=AS074-Av-L -N AS074  kore-align.sh
 
 ```
 
 #### sort reads bam file [kore-sort.sh](kore-sort.sh)
 ```
 qsub -o /mpba0/vcolonna/silvia/AS090.sort.out -e /mpba0/vcolonna/silvia/AS090.sort.err -v idrbam="AS090" -N sortAS090  kore-sort.sh
 
  ```
 
 #### remove PCR duplicates
  ```
 qsub -o /mpba0/vcolonna/silvia/AS090.mkdup.out -e /mpba0/vcolonna/silvia/AS090.mkdup.err -v idrbam="AS090" -N mkdupAS090  kore-markdup.sh
 
 ```
 
 #### index bam files
