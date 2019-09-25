
#### 1. download hg38 

#### 2. index the reference


#### 3. align reads to reference genome [kore-align](kore-align.sh)
```
 qsub -o /mpba0/vcolonna/silvia/AS074.out -e /mpba0/vcolonna/silvia/AS074.err -v gz1="AS074-Av-L_R1.fastq.gz",gz2="AS074-Av-L_R2.fastq.gz",idout="AS074",dir=AS074-Av-L -N AS074  kore-align.sh
 
 ```
 
 #### 4. sort reads bam file and make bam index  [kore-sort](kore-sort.sh)
 ```
 qsub -o /mpba0/vcolonna/silvia/AS090.sort.out -e /mpba0/vcolonna/silvia/AS090.sort.err -v idrbam="AS090" -N sortAS090  kore-sort.sh
 
 ```

#### 5. make stats from bam file  [kore-stats](kore-stats.sh)
 ```
 qsub -v id="AS064" -N statsAS064 kore-stats.sh
 
  ```
 
 #### 6. remove PCR duplicates [kore-markdup](kore-markdup.sh)
  ```
 qsub -o /mpba0/vcolonna/silvia/AS090.mkdup.out -e /mpba0/vcolonna/silvia/AS090.mkdup.err -v idrbam="AS090" -N mkdupAS090  kore-markdup.sh
 
 ```
 
 #### 7. variant calling [kore-freebayes](kore-freebayes.sh)
  ```
 qsub -o /mpba0/vcolonna/silvia/AS090.fby.out -e /mpba0/vcolonna/silvia/AS090.fby.err -v id="AS090",chr="chr1" -N fbyAS090  kore-freebayes.sh
 
 ```
 
 cycle for chr 1-22 
  ```
 for c in $(seq 1 22); do  echo qsub -o /mpba0/vcolonna/silvia/AS090.chr$c.fby.out -e /mpba0/vcolonna/silvia/AS090.chr$c.fby.err -v id="AS090",chr="chr$c" -N AS090.chr$c.fby  kore-freebayes.sh; done 
 
 ```


 cycle for ID for  chr 1-22 
  ```
  for id in AS006 AS054 AS064; do for c in $(seq 1 22); do qsub -o /mpba0/vcolonna/gianluca/$id.chr$c.fby.out -e /mpba0/vcolonna/gianluca/$id.chr$c.fby.err -v id="$id",chr="chr$c" -N $id.chr$c.fby  /mpba0/vcolonna/gianluca/kore-freebayes.sh; done; done
 
 ```

#### 8. variant calling chromosomes X and Y [kore-fb_chrXY](kore-fb_chrXY.sh)

 cycle for ID for  chrY
  ```
 for id in AS006 AS054 AS064 AS074 AS090 AS094 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrX.fby.out -e /mpba0/vcolonna/silvia/err/$id.chrX.fby.err -v id="$id",chr="chrX" -N $id.chrX.fb  kore-fb_chrXY.sh ; done
 
 ```
 cycle for ID chrY
 ```
 for id in AS054 AS090 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrY.fby.out -e /mpba0/vcolonna/silvia/err/$id.chrY.fby.err -v id="$id",chr="chrY" -N $id.chrY.fb  kore-fb_chrXY.sh ; done
 
 ```
 
 #### 9. make index from VCF file  [kore-vcfindex](kore-vcfindex.sh)
 ```
for id in AS006 AS054 AS064 AS090 AS094 AS074; do for c in $(seq 1 22); do  qsub -o /mpba0/vcolonna/gianluca/error_out/$id.chr$c.index.out -e /mpba0/vcolonna/gianluca/error_out/$id.chr$c.index.err -v id="$id",chr="chr$c" -N $id.chr$c.index  /mpba0/vcolonna/gianluca/kore-vcfindex.sh; done; done
 
  ```
 #### 10. VCF filter for quality >20 [kore-vcfFilterfb](kore-vcfFilterfb.sh)
 cycle for ID for chr 1-22
 ```
 for id in AS006 AS054 AS064 AS090 AS094 AS074; do for c in $(seq 1 22); do  qsub -o /mpba0/vcolonna/silvia/out/$id.chr$c.fb.filt.out -e /mpba0/vcolonna/silvia/err/$id.chr$c.filt.err -v id="$id",chr="chr$c" -N $id.chr$c.fb.filt  /mpba0/vcolonna/silvia/job/kore-vcfFilterfb ; done; done

```
cycle for ID for chr X 
```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrX.fb.filt.out -e /mpba0/vcolonna/silvia/err/$id.chrX.fb.filt.err -v id="$id",chr="chrX" -N $id.chrX.fb.filt  kore-vcfFilterfb.sh  ; done

```
cycle for ID for chr Y
```
for id in AS054 AS090 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrY.fb.filt.out -e /mpba0/vcolonna/silvia/err/$id.chrY.fb.filt.err -v id="$id",chr="chrY" -N $id.chrY.fb.filt  kore-vcfFilterfb.sh ; done

```

 #### 11. VCF normalize [kore-vtnormalizeFb](kore-vtnormalizeFb.sh)
 cycle for ID for chr 1-22
 ```
for id in AS006 AS054 AS064 AS090 AS094 AS074; do for c in $(seq 1 22); do  qsub -o /mpba0/vcolonna/silvia/out/$id.chr$c.fb.norm.out -e /mpba0/vcolonna/silvia/err/$id.chr$c.fb.norm.err -v id="$id",chr="chr$c" -N $id.chr$c.fb.norm  /mpba0/vcolonna/silvia/job/kore-vtnormalize.sh ; done; done

```
cycle for ID for chr X
```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrX.fb.norm.out -e /mpba0/vcolonna/silvia/err/$id.chrX.fb.norm.err -v id="$id",chr="chrX" -N $id.chrX.fb.norm  kore-vtnormalize.sh  ; done`

```
cycle for ID for chr Y
```
for id in AS054 AS090 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrY.fb.norm.out -e /mpba0/vcolonna/silvia/err/$id.chrY.fb.norm.err -v id="$id",chr="chrY" -N $id.chrY.fb.norm  kore-vtnormalize.sh  ; done

```


 #### 12. VCF Decompose (decomposes biallelic block substitutions into its constituent SNPs) [kore-decomposeblock](kore-decomposeblock.sh)
 cycle for ID for chr 1-22
 ```
 for id in AS006 AS054 AS064 AS090 AS094 AS074; do for c in $(seq 1 22); do  qsub -o /mpba0/vcolonna/silvia/out/$id.chr$c.fb.norm.dec.out -e /mpba0/vcolonna/silvia/err/$id.chr$c.fb.norm.dec.err -v id="$id",chr="chr$c" -N $id.chr$c.fb.dec  /mpba0/vcolonna/silvia/job/kore-decomposeblock.sh ; done; done
 
 ```
 cycle for ID for chr X
 ```
 for id in AS006 AS054 AS064 AS074 AS090 AS094 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrX.fb.dec.out -e /mpba0/vcolonna/silvia/err/$id.chrX.fb.dec.err -v id="$id",chr="chrX" -N $id.chrX.fb.dec  kore-decomposeblock.sh  ; done
 
```
cycle for ID for chr Y
```
for id in AS054 AS090 ;  do qsub -o /mpba0/vcolonna/silvia/out/$id.chrY.fb.dec.out -e /mpba0/vcolonna/silvia/err/$id.chrY.fb.dec.err -v id="$id",chr="chrY" -N $id.chrY.fb.dec  kore-decomposeblock.sh  ; done

```
 
 
#### 13. make stats from VCF file  [kore-bcfstats](kore-bcfstats.sh) - [kore-bcfQUAL](kore-bcfQUALstats.sh)
 ```
qsub -o /mpba0/vcolonna/gianluca/AS006.bcfstats.out -e /mpba0/vcolonna/gianluca/AS006.bcfstats.err -v id="AS006",chr="chr1" -N AS006bstats /mpba0/vcolonna/gianluca/kore-bcfstats.sh

 
  ```
