#!/bin/sh


#id=${id}
#chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

singularity exec -B /mpbastudies3/IMMA/samples/fb/:/data -B /mpbastudies3/IMMA/hg38p12/:/ref /mpba0/mpba-sw/biocontainers/vt.img vt normalize /data/mitochondrial/filtered/ALLgrep.chrM.fb.filt.vcf.gz -r /ref/hg38.p12.fa -o /data/mitochondrial/normalized/ALLgrep.chrM.fb.norm.vcf.gz
