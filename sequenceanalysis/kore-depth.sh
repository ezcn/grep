#!/bin/sh


id=${id}
#chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 




/mpba0/mpba-sw/samtools depth -q20 -Q20 -a /mpba0/vcolonna/IMMA/samples/bam/${id}.bam
