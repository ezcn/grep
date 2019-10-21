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


zcat /mpba0/vcolonna/IMMA/samples/fb/vep/mergedByChr/merged.chr22.fb.vep.vcf.gz | grep -v '#' | cut -f1,2 | awk '{print $1,$2,$2}'| tr " " "\t" > /mpba0/vcolonna/gianluca/TESI/hgdp/positionGREP.${chr}.bed

