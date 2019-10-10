#### 1. estimate read depth for all samples (only for positions of CGH probes)[kore-chrdepth.sh](jobs/kore-chrdepth.sh)
for chr 1-22
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/silvia/err/$id.chr$c.depth.err -o /mpba0/vcolonna/silvia/out/$id.chr$c.depth.out -v id="$id",chr="chr$c" -N $id.chr$c.depth kore-chrdepth.sh; done; done

```
for chrX
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrX.depth.err -o /mpba0/vcolonna/silvia/out/$id.chrX.depth.out -v id="$id",chr="chrX" -N $id.chrX.depth kore-chrdepth.sh; done; done

```
for chrY
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrY.depth.err -o /mpba0/vcolonna/silvia/out/$id.chrY.depth.out -v id="$id",chr="chrY" -N $id.chrY.depth kore-chrdepth.sh; done; done

```
#### 2.  post process samtoolsdepth output to make a bed file[kore-bedDepth](jobs/kore-bedDepth.sh)
for chr 1-22
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do for c in $(seq 1 22); do qsub -e /mpba0/vcolonna/silvia/err/$id.chr$c.bedDepth.err -o /mpba0/vcolonna/silvia/out/$id.chr$c.bedDepth.out -v id="$id",chr="chr$c" -N $id.chr$c.beddepth kore-bedDepth.sh; done; done

```
for chrX
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrX.bedDepth.err -o /mpba0/vcolonna/silvia/out/$id.chrX.bedDepth.out -v id="$id",chr="chrX" -N $id.chrX.beddepth kore-bedDepth.sh; done; done

```
for chrY
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/silvia/err/$id.chrY.bedDepth.err -o /mpba0/vcolonna/silvia/out/$id.chrY.bedDepth.out -v id="$id",chr="chrY" -N $id.chrY.beddepth kore-bedDepth.sh; done; done

```

#### 3. postprocess agilent probes file to make a be dwith bin identifier (midpoint)
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

#### 5. make stats from bam file  [kore-stats](jobs/kore-stats.sh)
