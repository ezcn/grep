#!/bin/sh


idfile=${id}

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


/mpba0/mpba-sw/samtools stats /mpbastudies3/IMMA/samples/${idfile}.raw.sorted.bam  > /mpba0/vcolonna/gianluca/${idfile}.stats.text

