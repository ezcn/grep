



cat /lustre/home/enza/variantcallingFB/ID_to_process | while read line ; do for c in $(seq 1 22) ; do echo condor_submit -name ettore -a "id=$line" -a "chr=chr$c"; done ; done



NB per ogni passaggio dopo il primo vcf, ogni file deve essere bgzip e tabix per poter funzionare il prossimo passaggio

~~~~~~~~~~~~~~~~~~~~~~~ step 0

lista degli ID dei campioni da usare di HGDP /lustre/home/enza/variantcallingFB/ID_to_process

~~~~~~~~ step 2

trasformare i CRAM in BAM: "condor-samtoolsCram2Bam.job"

~~~~~~~~ step 3

indicizzare i BAM per ottenere i BAI: "condor-samtoolsIndex.job"

~~~~~~~~ step 4

variant calling con freebayes di tutti i campioni HGDP e GREP usando le stesse informazioni : "condor-freebayes.job"

~~~~~~~~ step 5 

effettuare filtro per qualità: "condor-bcftoolsFilter.job"

~~~~~~~~ step 6 

merge di ogni ID per chr fino ad ottenere la cartella "/lustre/home/enza/variantcallingFB/merged": "condor-bcftoolsMerge.job" 
necessito di una folder contenente file di testo con i path di ogni file per ogni chr "/lustre/home/enza/variantcallingFB/ID_path"

~~~~~~~~ step 7 

è necessario fare VTnormalize e vcfuniq dei file merged : "condor-vtNormalize.job" & "condor-vcfuniq.job"

~~~~~~~~ step 8

forzare il variant calling dai file merged normalizzati e considerando solo i siti unici (vcfuniq) : "condor-freebayes_force.job"

~~~~~~~~ step 9 

fare lo step 6





VariantCalling HGDP - GREP

cat /lustre/home/enza/hgdpgrep/variantcallingFB/listID | while read line ;do for c in $(seq 1 22) ; do condor_submit -name ettore -a "id=$line" -a "chr=chr$c" -a "pop=hgdp" condor-freebayesV2_touch.job ; done; done 

cat /lustre/home/enza/hgdpgrep/variantcallingFB/listID | while read line ;do for c in $(seq 1 22) ; do condor_submit -name ettore -a "id=$line" -a "chr=chr$c" -a "pop=grep" condor-freebayesV2_touch.job ; done; done 

