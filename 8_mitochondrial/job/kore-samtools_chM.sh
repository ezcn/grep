#!/bin/sh
id=${id}
#chr=${chr}
# -q name of the queue to use
#$ -q bld.q
# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G
# number of threads in multi-threaded jobs
#$ -pe smp 2
# -M emailaddress@organization.xx, where to send email alerts
#$ -M giulianadang@gmail.com
# condition occurred under which email is to be sent, es=end,suspend
#$ -m es


/mpba0/mpba-sw/samtools view -b -h /mpba0/vcolonna/IMMA/samples/bam/${id}.bam chrM > /mpba0/vcolonna/giuliana/chM_bam/${id}_chrM.bam 
