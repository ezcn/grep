#!/bin/sh


id=${id}


# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=1G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es



/mpba0/mpba-sw/bedtools intersect -header -wa -wb -a /mpba0/vcolonna/gianluca/vcfFiltered/${id}.filteredVep -b /mpbastudies3/IMMA/amigo2/embryodev.bed /mpbastudies3/IMMA/amigo2/DDD.bed /mpbastudies3/IMMA/amigo2/gnomAD.bed /mpbastudies3/IMMA/amigo2/repro.bed > /mpbastudies3/IMMA/samples/intersectVep/${id}.intersect.vcf

