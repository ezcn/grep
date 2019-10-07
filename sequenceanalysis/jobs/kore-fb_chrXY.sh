#!/bin/sh

id=${id}
chr=${chr}

# set some parameters for qsub
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4

# -M emailaddress@organization.xx, where to send email alerts
#$ -M silvia.buonaiuto@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

~/bin/freebayes -f /mpbastudies3/IMMA/hg38p12/hg38.p12.fa -r ${chr} --cnv-map /mpba0/vcolonna/silvia/fb_XY.bed -b /mpbastudies3/IMMA/samples/bam/${id}.bam | ~/bin/bgzip > /mpbastudies3/IMMA/samples/fb/raw/${id}.${chr}.fb.vcf.gz
