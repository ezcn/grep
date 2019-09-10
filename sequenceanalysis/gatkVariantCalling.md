#### 1. download hg38 

#### 2. index the reference

#### 3. make dictionary (.dict) of the reference sequence [kore-dict.sh]
```
qsub -o /mpba0/vcolonna/silvia/dict.out -e /mpba0/vcolonna/silvia/dict.err -N dictref  kore-dict.sh

```

#### 4. align reads to reference genome [kore-align] (https://github.com/ezcn/imma/edit/master/sequenceanalysis/variantCalling.md)

#### 5. sort reads bam file and make bam index  [kore-sort] (https://github.com/ezcn/imma/edit/master/sequenceanalysis/variantCalling.md)

#### 6. variant calling


