#!/bin/sh

id=${id}
#chr=${chr}

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

~/bin/bcftools concat -o /mpba0/vcolonna/IMMA/samples/gatkConcat/${id}.g.fullvcf.gz -O z --threads 4 /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr1.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr2.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr3.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr4.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr5.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr6.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr7.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr8.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr9.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr10.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr11.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr12.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr13.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr14.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr15.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr16.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr17.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr18.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr19.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr20.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr21.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chr22.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chrX.g.vcf.gz /mpba0/vcolonna/IMMA/samples/gatkVCF/${id}.chrY.g.vcf.gz
