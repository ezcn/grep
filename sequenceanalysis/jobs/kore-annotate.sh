#!/bin/sh


id=${id}


# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4

# -M emailaddress@organization.xx, where to send email alerts
#$ -M silvia.buonaiuto@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es


~/bin/vcfannotate -b /mpba0/vcolonna/IMMA/amigo2/embryo.ok.bed -k annotation /mpba0/vcolonna/gianluca/vcfFiltered/${id}.filteredVep.vcf >/mpba0/vcolonna/IMMA/samples/annotateVEP/${id}.vepannotation.vcf
