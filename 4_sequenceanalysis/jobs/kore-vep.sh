#!/bin/sh


id=${id}


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


singularity run -B /mpbastudies3/IMMA/samples/:/data /mpba0/mpba-sw/biocontainers/vep-96.3.img --af_1kg --af_gnomad --appris --biotype --buffer_size 5000 --check_existing --distance 5000 --fork 4 --polyphen b --pubmed --regulatory --sift b --species homo_sapiens --symbol --tsl --cache --dir_cache /mpba0/vcolonna/vepcache/ --offline --vcf --variant_class -i /data/vcfconcat/${id}.fullVcf -o /data/vepVCF/${id}.fullvep


###for id in AS006 AS054 AS064 AS074 AS090 AS094 ;do echo qsub -e /mpba0/vcolonna/gianluca/$id.vep.err -o /mpba0/vcolonna/gianluca/$id.vep.out -v id="$id" -N vep$id /mpba0/vcolonna/gianluca/kore-vep.sh; done
