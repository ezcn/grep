#!/bin/sh


#id=${id}
#chr=${chr}

# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=2G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts
#$ -M flavia.villani@igb.cnr.it

# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

singularity exec /mpba0/mpba-sw/biocontainers/plink.img plink --vcf /mpba0/vcolonna/flavia/WGS/hgdp_wgs.20190516.full.chr9.vcf.gz --r2 --out /mpba0/vcolonna/flavia/WGS/WGS_Variants/HGDP_Soglia/chr9_Ld_2MB_Soglia --chr 9 --from-bp 79063076 --to-bp 83063077 --keep /mpba0/vcolonna/flavia/ldchr9/EUR.list --ld-window-kb 200
  
