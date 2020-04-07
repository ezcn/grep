After [GFAtoVCFodgi.py](/tesiFlavia/GFAtoVCFodgi.py) with detect bubble.
Now I convert Simulation to VCF without going through the graph for test result of MS.

I simulated two populations that are divided into three different times.


#### 1. Simulation sequences (MS) 

```
./ms 80 100 -t 11.2 -I 2 40 40 -g 1 44.36 -n 2 0.05 -eg 0.03125 1 0.0 -ej 0.09375 2 1 > outms40popT3  #15.000 generazioni 
```
```
./ms 80 100 -t 11.2 -I 2 40 40 -g 1 44.36 -n 2 0.05 -eg 0.03125 1 0.0 -ej 0.03125 2 1   > outms40T1pop  #5.000 generazioni
```
```
./ms 80 100 -t 11.2 -I 2 40 40 -g 1 44.36 -n 2 0.05 -eg 0.03125 1 0.0 -ej 0.0625 2 1  > outms40T2pop  #10.000 generazioni
```

#### 2. MStoVCF

[ms2vcf.py](/tesiFlavia/ms2vcf.py)  
 
#### 3. Bgzip and Tabix

bgzip and tabix for (ouputs of ms2vcf.py) each replicated: 

```
for f in ./*.vcf; do bgzip "$f"; done   
tabix: for f in ./*.vcf.gz; do tabix -h "$f"; done
```
#### 4. Extract individuals from the populations
```
parallel bcftools view -S pop_1.txt  ms_rep{1}.vcf.gz -o ms_rep{1}.bcf.vcf.gz ::: {1..100} #pop1
parallel bcftools view -S pop_2.txt ms_rep{1}.vcf.gz -o ms_rep{1}.bcf.vcf ::: {1..100}  #pop2 
```
#### 5. Calculate Fst statistics in pop1 and pop2
```
parallel vcftools --vcf ms_rep{1}.bcf.vcf --freq --out ms_rep{1}.bcf.pop1.vcf ::: {1..100}   #pop1
parallel vcftools --vcf ms_rep{1}.bcf.vcf --freq --out ms_rep{1}.bcf.pop2.vcf ::: {1..100}  #pop2
```
#### 6. Fst script
[CalculateFst.py](/tesiFlavia/CalculateFst.py)  

CalculateFSt for three different time.

#### 7. Plot Fst
```
w <- as.data.frame(rbind(
  t1, t2, t3
))
w$time <- c("T1", "T2", "T3")

library(ggplot2)
library(reshape2)
w.plot <- melt(w) 

p <- ggplot(aes(x=value, colour=time), data=w.plot)
p + geom_density()
```
The longer they separated, the greater the distribution of Fst. Simulations work well.




 
