#!/bin/sh


id=${id}
chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts
#$ -M silvia.buonaiuto@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

singularity exec -B /mpba0/vcolonna/silvia/:/data -B /mpbastudies3/IMMA/:/reference /mpba0/mpba-sw/biocontainers/gatk4.img gatk HaplotypeCaller -I /reference/samples/${id}.bam -L ${chr} -O /reference/samples/gatkVCF/${id}.${chr}.g.vcf.gz -R /reference/hg38p12/hg38.p12.fa
