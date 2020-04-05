# PCA

# 1. Obtain control data

### 1.1 Download control VCF (HGDP) from sanger institute.

```
wget ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr22.vcf.gz . 
```

# 2. Merge GREP and HGDP

### 2.1 There is an incompatibility between our file (GREP) and HGDP one. In our file the VariantCalling is performed with FREEBAYES, while HGDP is analyzed with GATK. We use VCFKEEPINFO & VCFKEEPGENO (from vcflib repo : https://github.com/vcflib/vcflib) to obtain two VCF that is possible to merge.

```
vcfkeepinfo /lustrehome/gianluca/PCA/data/hgdp_wgs.20190516.full.chr22.vcf.gz AC > /lustrehome/gianluca/PCA/data/keepinfo/hgdp_AC.vcf
```
``` 
vcfkeepgeno /lustrehome/gianluca/PCA/data/keepinfo/hgdp_AC.vcf.gz GT > /lustrehome/gianluca/PCA/data/keepgeno/hgdp_chr22.vcf
```

### 2.2 Indexing with BCFTOOLS

```
/bin/singularity exec -B /lustrehome/gianluca /lustre/home/enza/biocontainers/bcftools-1.9.img bcftools index /lustrehome/gianluca/PCA/data/keepgeno/hgdp_chr22.vcf.gz
```

### 2.3 Merge with BCFTOOLS

```
/bin/singularity exec -B /lustrehome/gianluca /lustre/home/enza/biocontainers/bcftools-1.9.img bcftools merge -0 -o /lustrehome/gianluca/PCA/data/merge/ALL_chr22.vcf.gz -O z /lustrehome/gianluca/PCA/data/keepgeno/hgdp_chr22.vcf.gz /lustrehome/gianluca/PCA/data/keepgeno/grep_chr22.vcf.gz
```

