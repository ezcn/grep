#!/bin/sh


id=${id}
#chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

/usr/bin/cat /mpba0/vcolonna/silvia/coverage/${id}.depth.out | /usr/bin/awk -f /mpba0/vcolonna/silvia/job/sexDetermination.awk > /mpba0/vcolonna/silvia/coverage/${id}.covcount.txt
