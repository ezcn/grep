#!/bin/sh


id=${id}
chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es


~/bin/vcfstats /mpbastudies3/IMMA/samples/vcf/${id}.${chr}.vcf.gz > /mpba0/vcolonna/gianluca/stats/${id}.${chr}.stats
