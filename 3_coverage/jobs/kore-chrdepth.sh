#!/bin/sh


id=${id}
chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 


/mpba0/mpba-sw/samtools depth -b /mpba0/vcolonna/silvia/coverage/agilentProbes/agilentProbes.bed -r ${chr} /mpba0/vcolonna/IMMA/samples/bam/${id}.bam >/mpba0/vcolonna/silvia/coverage/depthChr/${id}.${chr}.depth
