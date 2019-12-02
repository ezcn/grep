# Variant Effect Prediction

### Use [VEP](jobs/kore-vep.sh) for predict variant effect

```
for id in ALLgrep AS006 AS054 AS064 AS074 AS090 AS094 ;do qsub -e /mpba0/vcolonna/giuliana/errorFile/$id.vep.err -o /mpba0/vcolonna/giuliana/outputFile/$id.vep.out -v id="$id" -N vep$id /mpba0/vcolonna/giuliana/job/kore-vep.sh; done
```
