#!/bin/sh


#id=${id}
chr=${chr}

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

~/bin/vcftools --gzvcf /mpba0/vcolonna/IMMA/samples/fb/vep/mergedByChr/merged.${chr}.fb.vep.vcf.gz --positions /mpba0/vcolonna/gianluca/TESI/hgdp/positionALL/positionALLuniqB.${chr}.bed --recode --recode-INFO-all --out /mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/recodeGREP/merged.${chr}.B.fb.vep.vcf
