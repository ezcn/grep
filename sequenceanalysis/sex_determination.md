
### 1. Compute the depth at each position [kore-depth](jobs/kore-depth.sh)
cycle for ID
```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ; do qsub -o /mpba0/vcolonna/silvia/out/$id.depth.out -e /mpba0/vcolonna/silvia/err/$id.depth.err -v id=$id -N $id.depth kore-depth.sh ; done 

```
### 2. Awk script [sexDetermination.awk](sexDetermination.awk)

### 3. Run awk script [kore-awkscript](jobs/kore-awkscript.sh)
cycle for ID
```
for id in AS006 AS054 AS064 AS074 AS090 AS094 ; do qsub -o /mpba0/vcolonna/silvia/out/$id.awk.out -e /mpba0/vcolonna/silvia/err/$id.awk.err -v id=$id -N $id.awk kore-awkscript.sh ; done 

```
### 4. Sex determination from txt file
if xCoverage=autCoverage sample is female  
if xCoverage is half of autCoverage sample is male

our samples sex :  
AS006 FEMALE  
AS054 MALE  
AS064 FEMALE  
AS074 FEMALE  
AS090 MALE  
AS094 FEMALE  



