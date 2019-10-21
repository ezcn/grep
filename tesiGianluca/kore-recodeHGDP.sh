#!/bin/sh


#id=${id}
#chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

~/bin/vcftools --gzvcf /mpba0/vcolonna/IMMA/samples/hgdp_wgs/vep/hgdp_wgs.20190516.full.chr22.vep.vcf.gz --positions /mpba0/vcolonna/gianluca/TESI/hgdp/positionALLuniqB.chr22.bed --recode --recode-INFO-all --out /mpba0/vcolonna/gianluca/TESI/hgdp/merged.chr22.B.fb.vep.vcf
