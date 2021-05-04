This rep will keep all info about the phylogenetic analyses I will do for Timema.

##########################################################################################################################3

Pipeline
0 - get reads (OMG, this was a bit messy)
1 - filter reads
2 - map reads against Tce assembly


###################
## 0 - get reads ##

Reads were from mixed sources. Some reads were RADseq done on Chloe samples, others were Nosil GBS or whole genome sequencing.
GBS samples from the Nosil group were more or less 3x less (in terms total number of basis) than the Schwander group. 
(from a quick and dirty estimate)
So I will pool the Nosil GBS samples in groups of 3.

Nosil reads were trimmed to 60 bp
Schwander reads trimmed to 80 bp

```bash
cat inds_trim80 | while read line; do co reads/"$line"* . ; done
```


## Call SNPs ##

I will use VarScan (because there are several samples)

```bash
screen -S Tdi_tree
module load Bioinformatics/Software/vital-it 
module load UHTS/Analysis/samtools/1.10
samtools mpileup -f 1_Tdi_b3v08.fasta *filtered.bam > Tdi_tree.mpileup
java -jar VarScan.v2.3.9.jar mpileup2snp Tdi_tree.mpileup --min-coverage 10 --min-reads2 5 --output-vcf 1 > Tdi_tree.vcf
```
- VarScan options
USAGE: java -jar $VARSCAN/VarScan.jar mpileup2snp [mpileup file] OPTIONS
        mpileup file - The SAMtools mpileup file

        OPTIONS:
        --min-coverage  Minimum read depth at a position to make a call [8]
        --min-reads2    Minimum supporting reads at a position to call variants [2]
        --min-avg-qual  Minimum base quality at a position to count a read [15]
        --min-var-freq  Minimum variant allele frequency threshold [0.01]
        --min-freq-for-hom      Minimum frequency to call homozygote [0.75]
        --p-value       Default p-value threshold for calling variants [99e-02]
        --strand-filter Ignore variants with >90% support on one strand [1]
        --output-vcf    If set to 1, outputs in VCF format
        --variants      Report only variant (SNP/indel) positions (mpileup2cns only) [0]



## Filter by DP ##
First I need to get the number of reads per position, so I will plot the read distribution per variant and decide the cutoff according to the 10 and 99 quantiles.


In screen :
```bash
module add Bioinformatics/Software/vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.0.0
module load UHTS/Analysis/picard-tools/2.9.0

picard-tools CreateSequenceDictionary REFERENCE=1_Tdi_b3v08.fasta OUTPUT=1_Tdi_b3v08.dict

GenomeAnalysisTK VariantsToTable \
 -R 1_Tdi_b3v08.fasta \
 -V Tdi_tree.vcf \
 -F CHROM -F POS -GF GT -GF DP \
 --output Tdi.tree.DP.table
```

There was a problem while running GATK:
"The provided VCF file is malformed at approximately line number 15735: Duplicate allele added to VariantContext"
So I deleted the problematic line and re-run GATK (I had to deleted a few lines):

 ```bash
sed -e '15735d' Tdi_tree.vcf > Tdi_treeclean.vcf
sed -e '55609d' Tdi_treeclean.vcf > Tdi_treeclean2.vcf
sed -e '64952d' Tdi_treeclean2.vcf > Tdi_treeclean.vcf
sed -e '65870d' Tdi_treeclean.vcf > Tdi_treeclean2.vcf
sed -e '153305d' Tdi_treeclean2.vcf > Tdi_treeclean.vcf
sed -e '201687d' Tdi_treeclean.vcf > Tdi_treeclean2.vcf
sed -e '202618d' Tdi_treeclean2.vcf > Tdi_treeclean.vcf
sed -e '280256d' Tdi_treeclean.vcf > Tdi_treeclean2.vcf
```


```bash
GenomeAnalysisTK VariantsToTable \
 -R 1_Tdi_b3v08.fasta \
 -V Tdi_treeclean.vcf \
 -F CHROM -F POS -GF DP \
 --output Tdi.tree.DP.table

GenomeAnalysisTK VariantsToTable \
 -R 1_Tdi_b3v08.fasta \
 -V Tdi_treeclean2.vcf \
 -F CHROM -F POS -GF DP \
 --output Tdi.tree.DP.table
```

Check coverage to set the DP filter interval

```bash
for ((i=3; i<=96; i +=2)); do cut -f $i,$((i+1)) Tdi.tree.DP.table | awk '$1 != "./." {print $2}' > $i.DP; done
```
in R
```R
# define the file names list (10 samples here) 
nameList <- c()
for (i in 3:96) { # 21 - odd number for 10 samples 
  if (i %% 2 ==1) nameList <- append(nameList, paste(i, ".DP", sep=""))
}

qlist <- matrix(nrow = 47, ncol = 3) # define number of samples (10 samples here)
qlist <- data.frame(qlist, row.names=nameList)
colnames(qlist)<-c('5%', '10%', '99%')

png("Tdi_tree.DP.png", height=2000, width=1600)
par(mar=c(1, 1, 1, 1), cex=1.5, mfrow=c(12,4)) # define number of plots for your sample
for (i in 1:48) {
  DP <- read.table(nameList[i], header = T)
  qlist[i,] <- quantile(DP[,1], c(.05, .1, .99), na.rm=T)
  d <- density(DP[,1], from=0, to=100, bw=1, na.rm =T)
  plot(d, xlim = c(0,100), main=nameList[i], col="blue", xlab = dim(DP)[1], lwd=2)
  abline(v=qlist[i,c(1,3)], col='red', lwd=3)
}
dev.off()

write.table(qlist, "Tdi.tree.DP.percentiles.txt")
```

And now, filter the DP

```bash
## but first we need to create a dictionary for the fasta reference with Picard
module load UHTS/Analysis/picard-tools/2.9.0
picard-tools CreateSequenceDictionary REFERENCE=3_Tce_b3v08.fasta OUTPUT=3_Tce_b3v08.dict
picard-tools CreateSequenceDictionary REFERENCE=3_Tcm_b3v08.fasta OUTPUT=3_Tcm_b3v08.dict

## and samtools faidx
#samtools faidx 3_Tce_b3v08.fasta
#samtools faidx 2_Tcm_b3v08.fasta


##### - keep reads DP > 8 and DP < 200
#### gatk ####
GenomeAnalysisTK VariantFiltration \
 -R 1_Tdi_b3v08.fasta \
 -V Tdi_tree_samplename.vcf \
 --genotype-filter-expression 'DP<8||DP>200' \
 --genotype-filter-name 'DP_8-200' \
 --set-filtered-genotype-to-no-call\
 --output Tdi_tree_DPfilter.vcf

## OPTIONS:
## --set-filtered-genotype-to-no-call : used this to set the filtered genotypes (the ones that less than 8 > DP > 200)
## to ./.
## https://gatkforums.broadinstitute.org/gatk/discussion/7577/removing-variants-based-on-the-format-field

## I had a problem with this script and I was getting constant 'invalid argument' errors.
## I solved it by removing all spaces in between the formal arguments...
# http://selectvariants1.rssing.com/chan-9504484/all_p2.html
```

And now, finally remove low coverage individuals:

Check if some individuals have very low coverage
```bash
vcftools --vcf Tdi_tree_DPfilter.vcf --depth --out Tdi
```

And exclude the ones that have
```bash
Tps_04_45	0	-nan
Tps_04_60	0	-nan
TpsF_LMA_1_O15	0	-nan
TpsF_LMC_86_RO15	0	-nan
TpsF_LMC_88_RRO15	0	-nan
TpsM_LMA_10_M4	217	12.8525
TpsM_LMA_5_M3	0	-nan
TpsM_LMA_6_M1	191	12.911
TpsM_LMA_9_M1	621	13.5813
TpsM_LMC_87_RM1	434	12.1452
TpsM_LMC_89_RRM1	678	11.7168
Tps_04_27	2013	12.4471
Tps_04_5	2653	13.4365
TpsF_LMA_3_M5	2398	17.1914
TpsM_LMA_7_M1	2813	16.0128

vcftools --recode --recode-INFO-all --vcf Tdi_tree_DPfilter.vcf --remove-indv Tps_04_45 --remove-indv Tdi14_P1 --remove-indv Tps_04_60 --remove-indv TpsF_LMA_1_O15  --remove-indv TpsF_LMC_86_RO15 --remove-indv TpsF_LMC_88_RRO15 --remove-indv TpsM_LMA_10_M4 --remove-indv TpsM_LMA_5_M3 --remove-indv TpsM_LMA_6_M1 --remove-indv TpsM_LMA_9_M1 --remove-indv TpsM_LMC_87_RM1 --remove-indv TpsM_LMC_89_RRM1 --remove-indv Tps_04_27 --remove-indv Tps_04_5 --remove-indv TpsF_LMA_3_M5 --remove-indv TpsM_LMA_7_M1 --out Tdi.missfilt
```
```bash
vcftools --recode --recode-INFO-all --vcf Tdi.missfilt.recode.vcf --max-missing 0.75 --out Tdi_final
```


# Tree with Chloe/Guillaume dataset


17-1811

17-2076

17-1927

17-1810

17-2165

17-1793

17-2037

17-2005

17-1770



Files were trimmed with different parameters, so they have different endnames. I will run bwa in three separate screens regarding how the files were named.




```bash
screen -S trim80mapping
module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.8

### Create a BWA index in the genomic reference
# gunzip 1_Tdi_b3v08.fasta.gz

#samtools faidx 1_Tps_b3v08.fasta
#bwa index 1_Tps_b3v08.fasta
srun -p wally --nodes 1 --ntasks=4 --mem=20GB --time 10-00:00:0 --account=tschwand_default --pty bash 
for i in $(cat trim80mapping); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 1_Tps_b3v08.fasta $i.trim80.fq.gz $i-map.sam;
        echo $i;
done
```

```bash
screen -S tpsmapping
module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.8
samtools faidx 1_Tps_b3v08.fasta

### Create a BWA index in the genomic reference
# gunzip 1_Tdi_b3v08.fasta.gz

samtools faidx 1_Tps_b3v08.fasta
bwa index 1_Tps_b3v08.fasta

for i in $(cat Tpsmapping); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 1_Tps_b3v08.fasta $i.fq.gz-trim.fq.gz $i-map.sam;
        echo $i;
done
```
```bash
screen -S trim60mapping
module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.8
samtools faidx 1_Tps_b3v08.fasta

### Create a BWA index in the genomic reference
# gunzip 1_Tdi_b3v08.fasta.gz

samtools faidx 1_Tps_b3v08.fasta
bwa index 1_Tps_b3v08.fasta

for i in $(cat trim60mapping); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 1_Tps_b3v08.fasta $i.trim60.fq.gz $i-map.sam;
        echo $i;
done

# convert individual sam to bam, sort and index
for i in $(cat trim60mapping); do
        samtools view -S -h $i-map.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b -h -F 0x100 -q 30 - | samtools sort -o $i-filtered.bam - ;
        samtools index $i-filtered.bam;
done

# grep -v -e 'XA:Z:' -e 'SA:Z:' this will exclude all possible multi mapped reads.
# this step has to be done on the sam file, and not on the (binary) bam file.
# Info taken from here: https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
# -S	   ignored (input format is auto-detected)
#      --input-fmt-option OPT[=VAL]
#               Specify a single input file format option in the form of OPTION or OPTION=VALUE

# -q INT Skip alignments with MAPQ smaller than INT [0].
# -u Output uncompressed BAM. This option saves time spent on compression/decompression and is thus preferred when the output is piped to another samtools command
# -h Include the header in the output.

# -b Output in the BAM format.
# -F  Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
# -q Skip alignments with MAPQ smaller than INT [0].


# test quality of alignments
for i in $(cat trim60mapping); do
        samtools flagstat $i-filtered.bam > $i-flagstat.txt ;
done

```


