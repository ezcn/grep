#!/bin/sh


id=${id}
chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

singularity exec -B /mpbastudies3/IMMA/samples/fb/:/data /mpba0/mpba-sw/biocontainers/vt.img vt  decompose_blocksub  /data/normalized/${id}.${chr}.fb.norm.vcf.gz  -o /data/normalized/decomposed/${id}.${chr}.fb.norm.decompose.vcf.gz
