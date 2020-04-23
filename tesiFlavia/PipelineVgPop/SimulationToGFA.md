After [SimulationToVcf.md](SimulationToVcf.md) I convert MStoGfa for analysis of population genetics.

#### 1. Simulation sequences (MS) 

```
./ms 4 1 -t 11.2 -I 2 2 2 -g 1 44.36 -n 2 0.05 -eg 0.03125 1 0.0 -ej 0.03125 2 1 > out2pop.ms
```
 [MstoGfa.py](/tesiFlavia/MstoGfa.py)

I get the GFA but I need the whole sequence rebuilt and the links between the bubbles for use "odgi".


#### 2. Reconstruct sequence

```
./ms 4 1 -T -t 11.2 -I 2 2 2 -g 1 44.36 -n 2 0.05 -eg 0.03125 1 0.0 -ej 0.03125 2 1 > tree.ms 
```
I have Tree with the family history, but the reconstructed sequence is missing.
 
 ```
 seq-gen -mHKY -l 40 -s .2 -wa -z 783763255346462154 <tree.ms> seqwa.seqgen
```

```
 seq-gen -mHKY -l 40 -s .2 -z 783763255346462154 <tree.ms> seq.seqgen
```
 I use -z for set no random seed to try even without wa. 
 
 I get the same result, but -wa write Ancestral Sequences that was not in the tree of ms.
 
 #### 3. Seq-gene to GFA

 [SeqGenetoGfa.py](/tesiFlavia/SeqgenToGfa.py)
 
 
 #### 4. Calculate Allele Frequency on GFA :)
 
In this example there are only biallelic alleles (or reference or a different nucleotide). 
For calculate AF I count ATGC in a position and I check the reference base and calculate for each allele the frequency (count/num_haplotype).
 
 [calculateAlleleFrequency.py](/tesiFlavia/calculateAlleleFrequency.py)
 
 #### 5. Calculate Fst
 
