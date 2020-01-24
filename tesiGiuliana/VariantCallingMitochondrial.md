#### 1. variant calling 
for all samples [kore-freebayes](job/kore-freebayes_all_mt.sh)

```
qsub -e /mpba0/vcolonna/giuliana/fb.err -o /mpba0/vcolonna/giuliana/fb.out -N fb-chrM /mpba0/vcolonna/giuliana/job/kore-freebayes_all.sh
```

for single sample [kore-freebayes](job/kore-freebayes_each_mt_.sh)
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/giuliana/errorFile/$id.err -o /mpba0/vcolonna/giuliana/outputFile/$id.out -v id="$id" -N $id.chrM /mpba0/vcolonna/giuliana/kore-freebayes_each.sh; done
```


#### 2. VCF filter for quality >20 
for all samples [kore-vcfFilterfb](job/kore-vcfFilterfb_all_mt_.sh)
```
qsub -e /mpba0/vcolonna/giuliana/filt.err -o /mpba0/vcolonna/giuliana/filt.out -N filt-chrM /mpba0/vcolonna/giuliana/job/kore-vcfFilterfb.sh 
```
for single sample  [kore-vcfFilterfb](job/kore-vcfFilterfb_each_mt_.sh)
```
for id in AS006 AS054 AS064 AS074 AS090 AS094; do qsub -e /mpba0/vcolonna/giuliana/errorFile/$id_filt.err -o /mpba0/vcolonna/giuliana/outputFile/$id_filt.out -v id="$id" -N $id.chrM /mpba0/vcolonna/giuliana/job/kore-vcfFilterfb_each.sh; done
```

#### 3. VCF normalize 
for all samples [kore-vtnormalizeFb](job/kore-vtnormalize_all_mt_.sh)
```
qsub -e /mpba0/vcolonna/giuliana/norm.err -o /mpba0/vcolonna/giuliana/norm.out -N ALLnorm_chrM /mpba0/vcolonna/giuliana/job/kore-vtnormalize_all.sh
```
for single sample [kore-vtnormalizeFb](job/kore-vtnormalize_each_mt_.sh)
```
for id in AS006 AS074 AS090 AS094 AS054 AS064; do qsub -e /mpba0/vcolonna/giuliana/errorFile/$id.norm.err -o /mpba0/vcolonna/giuliana/outputFile/$id.norm.out -v id="$id" -N $id.norm /mpba0/vcolonna/giuliana/job/kore-vtnormalize_each.sh ; done
```
#### 4. VCF Decompose (decomposes biallelic block substitutions into its constituent SNPs) [kore-decomposeblock](job/kore-decomposeblock_mt.sh)
```
for id in ALLgrep AS006 AS074 AS090 AS094 AS054 AS064; do qsub -e /mpba0/vcolonna/giuliana/errorFile/$id.deco.err -o /mpba0/vcolonna/giuliana/outputFile/$id.deco.out -v id="$id" -N $id.norm /mpba0/vcolonna/giuliana/job/kore-decomposeblock.sh ; done
```
 #### 5. make stats from VCF file  [kore-bcfstats](job/kore-bcfstats_mt_decompose.sh) 
 ```
 for id in ALLgrep AS006 AS074 AS090 AS094 AS054 AS064; do qsub -e /mpba0/vcolonna/giuliana/errorFile/$id.decompose.bcfstats.err -o /mpba0/vcolonna/giuliana/outputFile/$id.decompose.bcfstats.out -v id="$id" -N $id.bcfstats /mpba0/vcolonna/giuliana/job/kore-bcfstats.mt.decompose.sh  ; done
 ```
 
