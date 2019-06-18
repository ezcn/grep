#!/bin/sh

## flr=file_row
idflrbam=${idrbam}

# set some parameters for qsub

### -N name of the job
###$ -N al$idsample

# -q name of the queue to use
#$ -q bld.q

# number of threads in multi-threaded jobs 
#$ -pe smp 16

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

~/bin/sambamba sort ${idflrbam}.raw.bam
##~/bin/sambamba markdup ${idflrbam}.raw.sorted.bam /mpbastudies3/IMMA/samples/${idflrbam}.bam
