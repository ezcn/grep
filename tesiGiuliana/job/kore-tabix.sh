#!/bin/sh


id=${id}
#chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=16G

# number of threads in multi-threaded jobs 
#$ -pe smp 1



~/bin/tabix -p vcf /mpba0/vcolonna/IMMA/samples/fb/mitochondrial/${id}.chrM.fb.vcf.gz
