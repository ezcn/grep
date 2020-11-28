#!/bin/bash

#Command to run the pipeline:
'''
~/repositories/grep/pipeline_grep_1.sh 
		SAMPLE_NAME 
		/data/resources/genome/FASTA_genome/genome_hg38.fasta 
		/path/to/output/directory
		/full/path/to/R1/AS0006_R1.fastq.gz 
		/full/path/to/R2/AS0006_R2.fastq.gz
'''

# Optional Arguments:
############################
stage=1                   ##
forcing=false             ##
target_bed="not_present"  ##
############################

while getopts "t:f:s:" flag; do
case "$flag" in
    f) forcing=$OPTARG;; ## this is to overwrite a file if exists
    s) stage=$OPTARG;; ## this is to skip some steps
	t) target_bed=$OPTARG;; ##this is if in the analysis is avaiable a bed target file for variant calling
esac
done

# Positional arguments
############################
sample_name=${@:$OPTIND:1} #
ref=${@:$OPTIND+1:1}       #
out_dir=${@:$OPTIND+2:1}   #
pairs1=${@:$OPTIND+3:1}    #
pairs2=${@:$OPTIND+4:1}    #
############################

## Auxiliar functions
log(){
     echo "[$(date --rfc-3339=seconds)]: $*"
}

# Starting pipeline
log "Sample name:         $sample_name"
log "Reference genome:    $ref"
log "Output directory:    $out_dir"
log "pair1 files:         $pairs1"
log "pair2 files:         $pairs2"
log "Overwriting files:   $forcing"

# outdir
log "Creating directory $out_dir"
mkdir -p $out_dir

# log dir
log_dir=$out_dir/log
mkdir -p $log_dir
log_sample=$log_dir/${sample_name}

# fastq files
pairs1_array=(`echo $pairs1 | sed 's/,/ /g'`)
pairs2_array=(`echo $pairs2 | sed 's/,/ /g'`)

# checking number of files
if [ ${#pairs1_array[@]} -ne ${#pairs2_array[@]} ]; then
    log "[ERROR] The number of files for pair1 and pair2 must be the same (${#pairs1_array[@]} vs ${#pairs2_array[@]})"
    exit
fi

log "Starting pipeline....................."

# lanes
lanes=${#pairs1_array[@]}
log "Number of lanes: $lanes"

# 1 Align reads to reference genome
#log "Starting alignment"
if [ "$stage" -eq 1 ]; then

    log "Doing stage "$stage
    create_folder_step1=$out_dir/bam
    mkdir -p $create_folder_step1
    out_step1=$out_dir/bam/${sample_name}

    if [ -f "$out_step1.raw.bam" ] && [ $forcing == "false" ]; then
        log "> File: $out_step1.raw.bam exist"
    else
        log "Alignment"
        bwa mem -M -t 30 -R "@RG\tID:${sample_name}\tSM:${sample_name}" $ref $pairs1 $pairs2 | samtools view -bSh - > $out_step1.raw.bam 
    fi

    if [ -f "$out_step1.raw.sorted.bam" ] && [ $forcing == "false" ]; then
        log "> File: $out_step1.raw.sorted.bam exist"
    else
        log "Sambamba sort bed file"
        sambamba sort -t 30 -m 25G --tmpdir /tmp -o $out_step1.raw.sorted.bam $out_step1.raw.bam
    fi

    if [ -f "$out_step1.bam" ] && [ $forcing == "false" ]; then
        log "> File: $out_step1.bam exist"
    else
        log "Sambamba remove PCR duplicates"
        sambamba markdup -t 30 -p --tmpdir /tmp --overflow-list-size 500000 $out_step1.raw.sorted.bam $out_step1.bam
    fi

fi # end stage


# 2 variant calling
#log "Starting variant calling"
if [ "$stage" -lt 4 ]; then

    log "Doing stage 2"
    create_folder_step2=$out_dir/vcf
    mkdir -p $create_folder_step2
    out_step2=$out_dir/vcf/${sample_name}

    for c in {1..22} X Y; do
        if [ -f "$out_step2.chr${c}.vcf.gz" ] && [ $forcing == "false" ]; then
            log "> File: $out_step2.chr${c}.vcf.gz exist"
        else
        	if [ $target_bed != "not_present" ]; then
            	log "freebayes chr${c}"
            	freebayes -f $ref -r chr${c} --min-base-quality 3 -F 0.15 -g 2000 -b $out_step1.bam | bgzip > $out_step2.chr${c}.vcf.gz &
            else
            	log "freebayes chr${c}"
            	freebayes -f $ref -r chr${c} --min-base-quality 3 -F 0.15 -g 2000 -t $target_bed -b $out_step1.bam | bgzip > $out_step2.chr${c}.vcf.gz &
            fi
        fi
    done
    wait
    for c in {1..22} X Y; do
        if [ -f "$out_step2.chr${c}.fb.filt.vcf" ] && [ $forcing == "false" ]; then
            log "> File: $out_step2.chr${c}.fb.filt.vcf exist"
        else
            log "vcffilter quality > 20 on chr${c}"
            zcat $out_step2.chr${c}.vcf.gz | vcffilter -f "QUAL > 20" | bgzip > $out_step2.chr${c}.fb.filt.vcf.gz &
        fi
    done
    wait
    for c in {1..22} X Y; do
        if [ -f "$out_step2.chr${c}.fb.filt.vcf" ] && [ $forcing == "false" ]; then
            log "> File: $out_step2.chr${c}.fb.filt.vcf exist"
        else
            log "normalization chr${c}"
            vt normalize -n $out_step2.chr${c}.fb.filt.vcf.gz -r $ref -o $out_step2.chr${c}.fb.norm.vcf.gz &
            #tabix -p vcf $out_step2.chr${c}.fb.norm.vcf.gz
        fi
    done
    wait
    for c in {1..22} X Y; do
        if [ -f "$out_step2.chr${c}.fb.norm.decompose.vcf.gz" ] && [ $forcing == "false" ]; then
            log "> File: $out_step2.chr${c}.fb.norm.decompose.vcf.gz exist"
        else
            log "decomposing chr${c}"
            vt decompose_blocksub  $out_step2.chr${c}.fb.norm.vcf.gz  -o $out_step2.chr${c}.fb.norm.decompose.vcf.gz &
            #tabix -p vcf $out_step2.chr${c}.fb.norm.decompose.vcf.gz
        fi
    done
    wait
    log "deleting tmp files..."
    rm $out_step2.chr*.fb.filt.vcf.gz
    rm $out_step2.chr*.fb.norm.vcf.gz
    log "Stage 2 Completed!"
    bcftools concat -Oz -o $out_step2.whole.vcf.gz $out_step2.chr*.fb.norm.decompose.vcf.gz
    rm $out_step2.chr*.fb.norm.decompose.vcf.gz
fi
#end stage


#Author info
#mail: adriano.demarino@gmail.com
