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


python3 /mpba0/vcolonna/gianluca/pythonScript/AFS-GREP_grepl.py -f /mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/merged.chr22.B.fb.vep.vcf.gz -o /mpba0/vcolonna/gianluca/TESI/hgdp/bothPos/grep.chr22.PyProcessed.tsv -v /mpba0/vcolonna/gianluca/pythonScript/csqimpact.tsv -e /mpba0/vcolonna/gianluca/junkfile/grep.py.err                     
