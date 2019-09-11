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



zcat /mpba0/vcolonna/IMMA/samples/gatkConcat/${id}.g.fullvcf.gz | ~/bin/vcffilter -f "QUAL > 20" > /mpba0/vcolonna/IMMA/samples/gatkFiltered/${id}.g.fullfiltered.vcf
