#### 0a Rationale for the nb of probes/intervals 
60,000 probes is the number of probes of the arrayCGH 

- same nb as array CGH: 60,000 probes of 100bp --> 6,000,000 bp --> 6Mbp/3,600Mbp = 0.0016 of the genome 
- 10% of the genome: 3,600,000 probes of 100bp --> 360,000,000 bp --> 360Mbp/3,600Mbp = 0.10 of the genome 

probes generation bedtools Random 

#### 0b Rationale for coverage 
- 30X --> 100% of reads 
- 15X --> 50% of reads 
- 1X --> 3.3% of reads 


#### 1. estimate read depth for all samples (only for positions of CGH probes)[kore-chrdepth.sh](jobs/kore-chrdepth.sh)
for chr 1-22
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/silvia/err/$id.chr$c.depth.err -o /mpba0/vcolonna/silvia/out/$id.chr$c.depth.out -v id="$id",chr="chr$c" -N $id.chr$c.depth kore-chrdepth.sh; done; done

```
for chrX
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrX.depth.err -o /mpba0/vcolonna/silvia/out/$id.chrX.depth.out -v id="$id",chr="chrX" -N $id.chrX.depth kore-chrdepth.sh; done

```
for chrY
```
for id in AS054 AS090; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrY.depth.err -o /mpba0/vcolonna/silvia/out/$id.chrY.depth.out -v id="$id",chr="chrY" -N $id.chrY.depth kore-chrdepth.sh; done

```
#### 2.  post process samtoolsdepth output to make a bed file[kore-bedDepth](jobs/kore-bedDepth.sh)
for chr 1-22
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/silvia/err/$id.chr$c.bedDepth.err -o /mpba0/vcolonna/silvia/out/$id.chr$c.bedDepth.out -v id="$id",chr="chr$c" -N $id.chr$c.beddepth kore-bedDepth.sh; done; done

```
for chrX
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrX.bedDepth.err -o /mpba0/vcolonna/silvia/out/$id.chrX.bedDepth.out -v id="$id",chr="chrX" -N $id.chrX.beddepth kore-bedDepth.sh; done

```
for chrY
```
for id in AS054 AS090; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrY.bedDepth.err -o /mpba0/vcolonna/silvia/out/$id.chrY.bedDepth.out -v id="$id",chr="chrY" -N $id.chrY.beddepth kore-bedDepth.sh; done

```

#### 3. postprocess agilent probes file to make a bed with bin identifier (midpoint)
```
cat agilentProbes.bed  | awk '{print $0, ($3-$2)/2+$2} ' | tr " " "\t"  > agilentProbesBin.bed
 
 ```
 
 #### 4. intersect bed file agilent probes and depth bed files[kore-bedtoolsIntersect](jobs/kore-bedtoolsIntersect.sh)
for chr 1-22
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/silvia/err/$id.chr$c.intersect.err -o /mpba0/vcolonna/silvia/out/$id.chr$c.intersect.out -v id="$id",chr="chr$c" -N $id.chr$c.int kore-bedtoolsIntersect.sh; done; done

```
for chrX
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrX.intersect.err -o /mpba0/vcolonna/silvia/out/$id.chrX.intersect.out -v id="$id",chr="chrX" -N $id.chrX.int kore-bedtoolsIntersect.sh; done; done

```
for chrY
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrY.intersect.err -o /mpba0/vcolonna/silvia/out/$id.chrY.intersect.out -v id="$id",chr="chrY" -N $id.chrY.int kore-bedtoolsIntersect.sh; done; done

```

#### 5. merge all sample by chr

#### 6. run R script to have files with mean depth for each sample[depthMean.R](depthMean.R)

#### 7. analize and plot depth data with R[coverageCNV.R](coverageCNV.R)
   ## args list: 
      [1]= input file
      [2]= gamma
      [3]=kmin
      [4]=experimentset(plot title)
      [5]=plot name in ggsave
      [6]=output file




