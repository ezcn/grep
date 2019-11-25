#!/bin/sh                                                                                                                                             
#id=${id}
#chr=${chr}

# -q name of the queue to use                                                                                                                  
#$ -q bld.q                                                                                                                                          
# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.                                                               
#$ -l mf=9G

# number of threads in multi-threaded jobs
#$ -pe smp 8

# -M emailaddress@organization.xx, where to send email alerts
#$ -M gianluca.damaggio@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

singularity exec /mpba0/mpba-sw/biocontainers/r-bioconductor-base2.img Rscript /mpba0/vcolonna/gianluca/TESI/hgdp/hgdp_grep_plot.R
