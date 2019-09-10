#### 1. download hg38 

#### 2. index the reference

#### 3. make dictionary (.dict) of the reference sequence [kore-dict](kore-dict.sh)
```
qsub -o /mpba0/vcolonna/silvia/dict.out -e /mpba0/vcolonna/silvia/dict.err -N dictref  kore-dict.sh

```

#### 4. align reads to reference genome [kore-align] (https://github.com/ezcn/imma/edit/master/sequenceanalysis/variantCalling.md)

#### 5. sort reads bam file and make bam index  [kore-sort] (https://github.com/ezcn/imma/edit/master/sequenceanalysis/variantCalling.md)

#### 6. variant calling [kore-gatk1](kore-gatk1.sh)
```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ; do for c in $(seq 1 22); do qsub -o /mpba0/vcolonna/silvia/$id.chr$c.gatk.out -e /mpba0/vcolonna/silvia/$id.chr$c.gatk.err -v id="$id",chr="chr$c" -N $id.chr$c.gatk  kore-gatk1.sh; done 

```



