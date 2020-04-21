#### 0a Download hg38 1000genomes reference vcf and plink map files

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}

```

```
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip 

```


#### 0b  Manipulate plink files (add 'chr' tag to the beginning of the line and store the output with suitable filename)

```
for CHR in {1..23}; do 
    cat plink.chr${CHR}.GRCh38.map | \
    sed 's/^/chr/' \
    > beagle_chr${CHR}_b38.map
done 
```
 
#### 1 Manipulate 1000genomes reference panel files 
- Generate a chromosome renaming file

```
for CHR in {1..23} X ; do 
    echo ${CHR} chr${CHR}
done >> chr_names.txt
```

- Multiple processing commands piped together

```
for CHR in {1..22} X; do
    bcftools annotate --rename-chrs chr_names.txt \
        ALL.chr${CHR}_GRCh38.genotypes.20170504.vcf.gz -Ou | \
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' -Ou | \
    bcftools norm -m -any -Ou | \
    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Ou | \
    bcftools norm -f Homo_sapiens_assembly38.fasta -d none -Ou | \
    bcftools view -m 2 -M 2 -Ou | \
    bcftools view -g ^miss -Oz -o 1000GP_chr${CHR}.vcf.gz"
done
```

- Fix the chromosome X ploidy to phased diploid (Requires a ploidy.txt file containing space-separated CHROM,FROM,TO,SEX,PLOIDY )

```
echo "chrX 1 156040895 M 2" > ploidy.txt
bcftools +fixploidy \
    1000GP_chrX.vcf.gz -Ov -- -p ploidy.txt | \
    sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | \
bcftools view -Oz -o 1000GP_chr23.vcf.gz"
```
- Remove duplicates

```
for CHR in {1..23}; do
    bcftools query -f '%ID\n' 1000GP_chr${CHR}.vcf.gz | \
    sort | uniq -d > 1000GP_chr${CHR}.dup_id

    if [[ -s 1000GP_chr${CHR}.dup_id ]]; then
    	bcftools view -e ID=@1000GP_chr${CHR}.dup_id \
    	1000GP_chr${CHR}.vcf.gz \
        -Oz -o 1000GP_filtered_chr${CHR}.vcf.gz
    else 
    	mv 1000GP_chr${CHR}.vcf.gz \
        1000GP_filtered_chr${CHR}.vcf.gz
    fi
done
```
-  Generate a list of the reference panel sample IDs

```
bcftools query -l 1000GP_AF_chr22.vcf.gz \
    > 1000GP_sample_IDs.txt

```
- Add grep variants to 1000genomes reference panel file

```
bcftools merge -Oz -0 -o /lustre/home/enza/referenceBeagle/1000GP_grep_$(chr).vcf.gz /lustre/home/enza/referenceBeagle/1000GP_filtered_$(chr).vcf.gz /lustre/home/enza/grep/variantcalling/allGrep.decomp.$(chr).vcf.gz

``` 

- Remove grep samples from merged vcf and leave 1000genomes individuals

```
view -S /lustre/home/enza/referenceBeagle/1000GP_sample_IDs.txt -o /lustre/home/enza/referenceBeagle/newref_$(chr).vcf /lustre/home/enza/referenceBeagle/1000GP_grep_$(chr).vcf.gz

```
- Replace missing genotypes with 0|0

```
/usr/bin/sed -i 's/\//\|/g' /lustre/home/enza/referenceBeagle/newref_$1.vcf && touch /lustrehome/silvia/junk/$1.reference.ok

```
- bgzip and tabix all vcf files

#### 2. Run Beagle for phasing


```
java -jar /lustrehome/gianluca/beagle/beagle.25Nov19.28d.jar gt=/lustre/home/enza/grep/variantcalling/allGrep.decomp.vcf.gz ref=/lustre/home/enza/referenceBeagle/newref_$(chr).vcf.gz out=/lustre/home/enza/referenceBeagle/phasing2/grep.$(chr).phased map=/lustre/home/enza/referenceBeagle/plinkMaps/beagle_$(chr)_b38.map

```

