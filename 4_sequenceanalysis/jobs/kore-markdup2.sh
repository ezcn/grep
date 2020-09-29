#!/bin/sh

## flr=file_row
idflrbam=${idrbam}

# set some parameters for qsub

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 16

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

##~/bin/sambamba sort ${idflrbam}.raw.bam
~/bin/sambamba markdup -t 16 -p --tmpdir /scratch /mpbastudies3/IMMA/samples/${idflrbam}.raw.sorted.bam /mpbastudies3/IMMA/samples/${idflrbam}.bam
