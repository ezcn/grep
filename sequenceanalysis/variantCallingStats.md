# Summary statistics from mapping - (from bam files) [script.R](statsRscript/BAMstats.R)

### [Insert size](img/insert_size.png)
### [ACGT content per cycle](img/base_comp.png)
### [GC Content](img/GCcontent.png)
### [Indel distribution](img/indels.png)
### [Coverage distribution](img/coverageAS074.png)

# Summary statistics of variant calling (form vcf files) [script.R](statsRscript/SNstatsVCF.R)

### [Number of MNPs - full](img/fullSN-numberofMNPs.png) : number of rows with a MNP, such as CC>TT
### [Number of MNPs - chr 22](img/SN-numberofMNPs.png)
### [Number of SNPs - full](img/fullSN-numberofSNPs.png) :  number of rows with a SNP
### [Number of SNPs - chr 22](img/SN-numberofSNPs.png)
### [Number of Indels - full](img/fullSN-numberofindels.png) : number of rows with an indel
### [Number of Indels - chr 22](img/SN-numberofindels.png) 
### [Number of multiallelic SNPs sites - full](img/fullSN-numberofmultiallelicSNPsites.png) number of rows with multiple alternate alleles, all SNPs
### [Number of multiallelic SNPs sites - chr 22](img/SN-numberofmultiallelicSNPsites.png) 
### [Number of multiallelic sites - full](img/fullSN-numberofmultiallelicsites.png) number of rows with multiple alternate alleles
### [Number of multiallelic sites - chr 22](img/SN-numberofmultiallelicsites.png) 
### [Number of no-ALTs - full](img/fullSN-numberofno-ALTs.png) reference-only sites, ALT is either "." or identical to REF
### [Number of no-ALTs - chr 22](img/SN-numberofno-ALTs.png)
### [Number of Others - full](img/fullSN-numberofothers.png) number of rows with other type, for example a symbolic allele or a complex substitution, such as ACT>TCGA
### [Number of Others - chr 22](img/SN-numberofothers.png)
### [Number of Records - full](img/fullSN-numberofrecords.png) number of data rows in the VCF
### [Number of Records - chr 22](img/SN-numberofrecords.png) 

# Substitution types (from vcf files) [script.R](statsRscript/STstatsVCF.R)

### [Substitution types -full] miss
### [Substitution types - chr 22](img/ST-Substitutiontypes.png) 

# Transitions/transversions (from vcf files) [script.R](statsRscript/TSTVstatsVCF.R)

### [Transitions - full](img/fulltransition.png)
### [Transitions - chr 22](img/transition.png)
### [Transversion - full](img/fulltransversion.png)
### [Transversion - chr 22](img/transversion.png)
### [TSTV ratio - full](img/fulltstvratio.png)
### [TSTV ratio - chr 22](img/tstvratio.png)

# Depth (from vcf files) [script.R](statsRscript/DEPTHstatsVCF.R)


### [Depth - chr 22](img/depth-quality1.png)  number of reads, on average, thet are likely to be aligned at a given reference base position.
 
