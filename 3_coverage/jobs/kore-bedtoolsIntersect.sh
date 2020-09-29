d=${id}
chr=${chr}

# set some parameters for qsub
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4



/mpba0/mpba-sw/bedtools intersect -a /mpba0/vcolonna/silvia/coverage/depthChr/${id}.${chr}depth.bed -b /mpba0/vcolonna/silvia/coverage/depthChr/agilentProbesBin.bed  -wa -wb >/mpba0/vcolonna/silvia/coverage/files/${id}.${chr}.depthBin.bed

