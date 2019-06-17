#!/bin/sh

# set some parameters for qsub

# -N name of the job
#$ -N alhg38

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
##bwa index hg38.p12.fa

## per allineare dobbiamo eseguire questo comando che alla fine ci da un file bam

bwa mem -t 16 -R '@RG\tID:AS006\tSM:AS006' /mpbastudies3/IMMA/hg38/hg38.p12.fa /mpbastudies3/181113_Vincenza-Colonna/AS006-Av-L/AS006-Av-L_R1.fastq.gz  /mpbastudies3/181113_Vincenza-Colonna/AS006-Av-L/AS006-Av-L_R2.fastq.gz | samtools view -b - > /mpbastudies3/IMMA/samples/AS006.raw.bam

## questo file possiamo sortarlo e marcare quelli che sembrano duplicati ridondanti della PCR  

##sambamba sort "nomefileoutput".raw.bam
##sambamba markdup "fileoutputsenzaduplicati".raw.sorted.bam "".bam
