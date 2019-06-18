
#### download hg38 

#### index the reference


#### align reads to reference genome [kore-align](kore-align.sh)
```
 qsub -o /mpba0/vcolonna/silvia/AS074.out -e /mpba0/vcolonna/silvia/AS074.err -v gz1="AS074-Av-L_R1.fastq.gz",gz2="AS074-Av-L_R2.fastq.gz",idout="AS074",dir=AS074-Av-L -N AS074  kore-align.sh
 
 ```
 
 #### sort reads bam file and make bam index  [kore-sort](kore-sort.sh)
 ```
 qsub -o /mpba0/vcolonna/silvia/AS090.sort.out -e /mpba0/vcolonna/silvia/AS090.sort.err -v idrbam="AS090" -N sortAS090  kore-sort.sh
 
  ```
 
 #### remove PCR duplicates [kore-markdup](kore-markdup.sh)
  ```
 qsub -o /mpba0/vcolonna/silvia/AS090.mkdup.out -e /mpba0/vcolonna/silvia/AS090.mkdup.err -v idrbam="AS090" -N mkdupAS090  kore-markdup.sh
 
 ```
 
 #### variant calling 
  ```
 qsub -o /mpba0/vcolonna/silvia/AS090.fby.out -e /mpba0/vcolonna/silvia/AS090.fby.err -v id="AS090",chr="chr1" -N fbyAS090  kore-freebayes.sh
 
 ```
 
 cycle for chr 1-22 
  ```
 for c in $(seq 1 22); do  echo qsub -o /mpba0/vcolonna/silvia/AS090.fby.out -e /mpba0/vcolonna/silvia/AS090.fby.err -v id="AS090",chr="chr$c" -N fbyAS090  kore-freebayes.sh; done 
 
 ```


 cycle for ID for  chr 1-22 
  ```
 for id in AS006 AS074 ; do for c in $(seq 1 22); do  echo qsub -o /mpba0/vcolonna/silvia/$id.fby.out -e /mpba0/vcolonna/silvia/$id.fby.err -v id="$id",chr="chr$c" -N fby$id  kore-freebayes.sh; done; done  
 
 ```

