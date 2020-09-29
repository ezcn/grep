#!/bin/sh


id=${id}


# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts


# condition occurred under which email is to be sent, es=end,suspend
#$ -m es


~/bin/bcftools stats /mpbastudies3/IMMA/samples/fb/mitochondrial/${id}.chrM.fb.vcf.gz > /mpba0/vcolonna/IMMA/samples/fb/mitochondrial/VCFstats_mt/${id}.chrM.bcf-stats
