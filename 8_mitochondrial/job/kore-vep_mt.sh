#!/bin/sh


id=${id}


# -q name of the queue to use
#$ -q bld.q

# -l mf=amount of memory requested (this is a MANDATORY parameter), use carefully.
#$ -l mf=4G

# number of threads in multi-threaded jobs 
#$ -pe smp 4 

# -M emailaddress@organization.xx, where to send email alerts


# condition occurred under which email is to be sent, es=end,suspend
#$ -m es

singularity run -B /mpbastudies3/IMMA/samples/:/data /mpba0/mpba-sw/biocontainers/vep-96.3.img --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --minimal --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /mpba0/vcolonna/vepcache/ --offline --vcf --variant_class -i /data/fb/mitochondrial/normalized/decomposed/${id}.chrM.fb.norm.decompose.vcf.gz -o /data/fb/mitochondrial/vep/${id}.chrM.fb.norm.decompose.vep.vcf
