#!/bin/sh

#id=${id}
#chr=${chr}

# set some parameters for qsub
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=16G

# number of threads in multi-threaded jobs 
#$ -pe smp 1

/mpba0/mpba-sw/bedtools random -n 360000 -l 100 -seed 239 -g /mpba0/vcolonna/silvia/coverage/hg38.p12.bed >/mpba0/vcolonna/silvia/random_1p

#/mpba0/mpba-sw/bedtools random -n 3600000 -l 100 -seed 415 -g /mpba0/vcolonna/silvia/coverage/hg38.p12.bed >/mpba0/vcolonna/silvia/random_3.6M

#/mpba0/mpba-sw/bedtools random -n 60000  -seed 899 -g /mpba0/vcolonna/silvia/hg38.p12.bed >/mpba0/vcolonna/silvia/hg38p12.bed
