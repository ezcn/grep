#!/bin/sh

idsample=${id}
chr=${chr}

# set some parameters for qsub
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 16

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

~/bin/freebayes -f /mpbastudies3/IMMA/hg38/hg38.p12.fa -r ${chr} --gvcf -g 2000  /mpbastudies3/IMMA/samples/${idsample}.bam | bgzip > /mpbastudies3/IMMA/samples/${idsample}.vcf.gz
