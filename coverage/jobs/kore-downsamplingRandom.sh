#!/bin/sh

id=${id}
#chr=${chr}

# set some parameters for qsub
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=16G

# number of threads in multi-threaded jobs 
#$ -pe smp 1

/mpba0/mpba-sw/samtools view -h -s 0.033 -b /mpbastudies3/IMMA/samples/bam/${id}.bam >/mpbastudies3/IMMA/samples/bam/${id}.random.1X.bam 
