#!/bin/sh

# set some parameters for qsub

# -N name of the job
#$ -N indexhg38

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

## dobbiamo prima indicizzare 
bwa index /mpbastudies3/IMMA/hg38/hg38.p12.fa
