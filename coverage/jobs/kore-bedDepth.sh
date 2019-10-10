#!/bin/sh

id=${id}
chr=${chr}

# set some parameters for qsub 
#$ -q bld.q 

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G
    
# number of threads in multi-threaded jobs 
#$ -pe smp 4
    
    
/usr/bin/cat /mpbastudies3/IMMA/samples/coverage/${id}.${chr}.depth | /usr/bin/awk '{print $1, $2, $2, $3 }' | /usr/bin/tr " " "\t">/mpbastudies3/IMMA/samples/coverage/${id}.${chr}.depth.bed
