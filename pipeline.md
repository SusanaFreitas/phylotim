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

Individuals:
17-1811
17-2076
17-1927
17-1810
17-2165
17-1793
17-2037
17-2005
17-1770



Files were trimmed with different parameters, so they have different endnames. I will run bwa in three separate screens regarding how the files were named:

### Map the reads

- TRIMTpsMAPPING

```bash
screen -S tpsmapping
module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.8
# samtools faidx 1_Tps_b3v08.fasta

# cd /scratch/wally/FAC/FBM/DEE/tschwand/nephus/phylotim/02-mapping/reads

### Create a BWA index in the genomic reference
# gunzip 1_Tdi_b3v08.fasta.gz

samtools faidx 1_Tps_b3v08.fasta
bwa index 1_Tps_b3v08.fasta

for i in $(cat Tpsmapping); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 1_Tps_b3v08.fasta $i.fq.gz-trim.fq.gz > $i-map.sam;
        echo $i;
done

# convert individual sam to bam, sort and index
for i in $(cat Tpsmapping); do
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
for i in $(cat Tpsmapping); do
        samtools flagstat $i-filtered.bam > $i-flagstat.txt ;
done
```



- TRIM60MAPPING


```bash
screen -S trim60mapping
module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.8
# samtools faidx 1_Tps_b3v08.fasta

# cd /scratch/wally/FAC/FBM/DEE/tschwand/nephus/phylotim/02-mapping/reads


### Create a BWA index in the genomic reference
# gunzip 1_Tdi_b3v08.fasta.gz

#samtools faidx 1_Tps_b3v08.fasta
#bwa index 1_Tps_b3v08.fasta
srun -p wally --nodes 1 --ntasks=4 --mem=20GB --time 10-00:00:0 --account=tschwand_default --pty bash 

for i in $(cat trim60mapping); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 1_Tps_b3v08.fasta $i.trim60.fq.gz > $i-map.sam;
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

- TRIM80MAPPING

```bash
screen -S tpsmapping
module add Bioinformatics/Software/vital-it
module load UHTS/Aligner/bwa/0.7.15
module load UHTS/Analysis/samtools/1.8
# samtools faidx 1_Tps_b3v08.fasta

# cd /scratch/wally/FAC/FBM/DEE/tschwand/nephus/phylotim/02-mapping/reads

### Create a BWA index in the genomic reference
# gunzip 1_Tdi_b3v08.fasta.gz

#samtools faidx 1_Tps_b3v08.fasta
#bwa index 1_Tps_b3v08.fasta

srun -p wally --nodes 1 --ntasks=4 --mem=20GB --time 10-00:00:0 --account=tschwand_default --pty bash 

for i in $(cat trim80mapping); do
        bwa mem -t 4 -M -R "@RG\tID:$i\tSM:$i" 1_Tps_b3v08.fasta $i.trim80.fq.gz > $i-map.sam;
        echo $i;
done

# convert individual sam to bam, sort and index
for i in $(cat trim80mapping); do
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
for i in $(cat trim80mapping); do
        samtools flagstat $i-filtered.bam > $i-flagstat.txt ;
done
```


### Call SNPs

I will use VarScan (because there are several samples)

```bash
screen -S Tdi_tree
# cd /scratch/wally/FAC/FBM/DEE/tschwand/nephus/phylotim/02-mapping/reads
module load Bioinformatics/Software/vital-it 
module load UHTS/Analysis/samtools/1.10
# samtools mpileup -f 1_Tdi_b3v08.fasta *filtered.bam --vcf-sample-list Tdi.list > Tdi_tree.mpileup
samtools mpileup -f 1_Tdi_b3v08.fasta *filtered.bam > Tdi_tree.mpileup
java -jar VarScan.v2.3.9.jar mpileup2snp Tdi_tree.mpileup --min-coverage 10 --min-reads2 5 --output-vcf 1 > Tdi_tree.vcf
#
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
        --vcf-sample-list Add a list of sample names to use in the VCF header. This list should be in plain text, one sample per line, in the order that samples appear in the raw mpileup input. This option is only available after 1.13 version




Working directory:
```bash
cd /scratch/wally/FAC/FBM/DEE/tschwand/nephus/phylotim/02-mapping/reads 

##### - keep reads DP > 8 and DP < 200
#### gatk ####
GenomeAnalysisTK VariantFiltration \
 -R 1_Tdi_b3v08.fasta \
 -V Tdi_tree.vcf \
 --genotype-filter-expression 'DP<8||DP>200' \
 --genotype-filter-name 'DP_8-200' \
 --set-filtered-genotype-to-no-call\
 --output Tdi_tree_DPfilter.vcf
``` 
 

Vcf editing (change the samplenames):
```bash
# get specific column number with awk:

awk -F " " '{print $9}' name.file > samples

I got a list of sample names. I want to substitute the newline by a tab (or whatever it is between the sample names in the vcf file)

tr "\n" "\t" < samples > tabsamples



head -30 Tdi_tree_DPfilter.vcf > head.Tdi.DPfilter
cat head.Tdi.DPfilter new_line_sample_names 
cat head.Tdi.DPfilter new_line_sample_names > head2.Tdi.DPfilter

If you want to delete lines 1 through 31:

sed -e '1,31d' Tdi_tree_DPfilter.vcf > 
cat head2.Tdi.DPfilter baseTdi_DPfilter > Tdi_DPfilter_ newnames.vcf
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
Tps_LMA_62	0	-nan
TpsM_LMA_5_M3	0	-nan
Tps_SRA_78	0	-nan
Tps_SRA_79	0	-nan
Tps_SRA_80	0	-nan


vcftools --recode --recode-INFO-all --vcf Tdi_DPfilter_newnames.vcf --remove-indv Tps_04_45 --remove-indv Tps_04_60 --remove-indv TpsF_LMA_1_O15 --remove-indv TpsF_LMC_86_RO15 --remove-indv TpsF_LMC_88_RRO15 --remove-indv Tps_LMA_62 --remove-indv TpsM_LMA_5_M3 --remove-indv Tps_SRA_78 --remove-indv Tps_SRA_79 --remove-indv Tps_SRA_80 --out Tdi.missfilt
```

Keep only SNPs with max allowed missing data:
```bash
vcftools --recode --recode-INFO-all --vcf Tdi.missfilt.recode.vcf --max-missing 0.75 --out Tdi_final75
# After filtering, kept 29179 out of a possible 433645 Sites

vcftools --recode --recode-INFO-all --vcf Tdi.missfilt.recode.vcf --max-missing 0.50 --out Tdi_final50
# After filtering, kept 108896 out of a possible 433645 Sites

```
Now we will make the pairwise distance tree, in R:

```R
################################################################
######################## ADEGENET ##############################
################################################################

library(vcfR)
library(adegenet)
#library(adegraphics)
#library(pegas)
#library(StAMPP)
#library(lattice)
#library(gplots)
#library(ape)
#library(ggmap)
library(poppr)

#### set working directory
setwd("/home/cravinhos/Documents/cryptic_gene_flow/tree_guillaume/3_distancetree_R/")
## 75% max missing data

tdivcf <- read.vcfR("Tdi_final75.recode.vcf", verbose = FALSE )
tdi <- vcfR2genind(tdivcf)


## see individual names
indNames(tdi)


## eliminate monikensis individual (from ForSale)
indNames(tdi)
Dat1 <- tdi[indNames(tdi) != "Tps_04_24"]
Dat1 -> tdi

## set pop tdi
tdi@pop <- as.factor(c("madonna", "madonna", "Philo", "Horseranch", "Horseranch", "Horseranch", 
                        "Horseranch", "Horseranch",
                        "Fort Bragg", "Fort Bragg", "NA", "Fort Bragg", "Fort Bragg", "Fort Bragg", 
                        "Fish", "Fish",
                        "Iverson", "Iverson", "Iverson", "Iverson", "Orr", "Iverson", "Ft Ross",
                        "Ft Ross", "Ft Ross", "Ft Ross",
                        "Ft Ross", "Wpt128", "Wpt128", "Swanton", "Swanton", "Orr", "Swanton", 
                        "Swanton", "Swanton", "Big pullout vista",
                        "Big pullout vista", "Big pullout vista", "Big pullout vista", "Orr", "Orr", 
                        "ORR_11", "Summit rd", "Summit rd",
                        "Summit rd", "Summit rd", "Fish-rock netbag 3", "Fish-rock netbag 1", 
                        "Fish-rock netbag 2", "Fish-rock netbag 9",
                        "Orr", "madonna", "madonna", "madonna", "Manch_5", "Manch_5", "NOSIL", 
                        "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL",
                        "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross", "cross",
                        "Manch_4", "manch_1", "manch_1", "Manch_1", "Manch_3", "Manch_1", "Manch_1", 
                        "Manch_1", "MANCH_5", "ORR_1", "MANCH_1",
                        "ORR_8", "MANCH_2", "ORR_1", "MANCH_1", "ORR_8", "MANCH_1", "ORR_1", "MANCH_6",
                        "MANCH_6", "MANCH_12", "MANCH_3", "MANCH_9",
                        "MANCH_2", "MANCH_6", "MANCH_12", "MANCH_1", "MANCH_1", "MANCH_3", "MANCH_12",
                        "MANCH_1", "MANCH_3", "MANCH_4", "MANCH_6",
                        "ORR_3", "MANCH_1", "MANCH_1", "MANCH_4", "MANCH_7", "MANCH_5", "MANCH_4", 
                        "MANCH_5", "MANCH_2", "MANCH_2", "MANCH_1",
                        "MANCH_1", "MANCH_7", "MANCH_7", "MANCH_2", "MANCH_4", "MANCH_4", "MANCH_5",
                        "MANCH_5", "MANCH_3", "MANCH_2", "MANCH_4",
                        "MANCH_7", "MANCH_7", "MANCH_4", "MANCH_7", "ORR_5", "ORR_8", "MANCH_1", "ORR_1", 
                        "ORR_4", "MANCH_3", "MANCH_1", "MANCH_3",
                        "MANCH_4", "MANCH_10", "cross", "cross", "cross", "pop inconnue", "MANCH_2", "ORR_4"))
tdi@pop

# confirm if we are doing it correctly
chr_unq <- unique(tdi@pop)
length(chr_unq)

##### PCoA - ADEGENET #######
## make PCA
x.cows <- tab(tdi, freq=TRUE, NA.method="mean")
pca.cows <- dudi.pca(x.cows, center=TRUE, scale=TRUE)

# 5  Multivariate analyses
# 5.1  Principal Component Analysis (PCA)
# Principal Component Analysis (PCA) is the amongst the most common multivariate
# analyses used in genetics.  Running a PCA ongenindobject is straightforward. 
# One needs to firstextract allelic data (as frequencies) and replace missing values 
# using the accessortabandthen use the PCA procedure (dudi.pca).  Let us use this
# approach on themicrobovdata.
# Let us first load the data:
# data(microbov)
# x.cows <- tab(microbov, freq=TRUE, NA.method="mean")
# pca.cows <- dudi.pca(x.cows, center=TRUE, scale=FALSE)
#The function dudi.pca displays a barplot of eigenvalues (thescreeplot) and asks 
# for a numberof retained principal components.  In general,  eigenvalues represent
# the amount of geneticdiversity  -  as  measured  by  the  multivariate  method  being
# used  -  represented  by  each principal component (PC). Here,  each eigenvalue is the
# variance of the corresponding PC.A sharp decrease in the eigenvalues is usually
# indicative of the boundaries between relevantstructures and random noise. 
# Here, how many axes would you retain?
# Tge = 10
# tdi = 20
# Tsi = 5

png("tdi_all_PCA.png", width = 580, height = 500)
s.label(pca.cows$li)
s.class(pca.cows$li, fac=pop(tdi), col=funky(41))
dev.off()


# tdi
indNames(tdi)
fac.score <- factor(c("madonna", "madonna", "Philo", "Horseranch", "Horseranch", "Horseranch", 
                        "Horseranch", "Horseranch",
                        "Fort Bragg", "Fort Bragg", "NA", "Fort Bragg", "Fort Bragg", "Fort Bragg", "Fish", "Fish",
                        "Iverson", "Iverson", "Iverson", "Iverson", "Orr", "Iverson", "Ft Ross", "Ft Ross",
                        "Ft Ross", "Ft Ross",
                        "Ft Ross", "Wpt128", "Wpt128", "Swanton", "Swanton", "Orr", "Swanton",
                        "Swanton", "Swanton", "Big pullout vista",
                        "Big pullout vista", "Big pullout vista", "Big pullout vista", "Orr", "Orr", 
                        "ORR_11", "Summit rd", "Summit rd",
                        "Summit rd", "Summit rd", "Fish-rock netbag 3", "Fish-rock netbag 1", 
                        "Fish-rock netbag 2", "Fish-rock netbag 9",
                        "Orr", "madonna", "madonna", "madonna", "Manch_5", "Manch_5", "NOSIL", 
                        "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL",
                        "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross",
                        "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", "cross", 
                        "cross", "cross", "cross", "cross",
                        "Manch_4", "manch_1", "manch_1", "Manch_1", "Manch_3", "Manch_1", "Manch_1", 
                        "Manch_1", "MANCH_5", "ORR_1", "MANCH_1",
                        "ORR_8", "MANCH_2", "ORR_1", "MANCH_1", "ORR_8", "MANCH_1", "ORR_1", "MANCH_6", 
                        "MANCH_6", "MANCH_12", "MANCH_3", "MANCH_9",
                        "MANCH_2", "MANCH_6", "MANCH_12", "MANCH_1", "MANCH_1", "MANCH_3", "MANCH_12", 
                        "MANCH_1", "MANCH_3", "MANCH_4", "MANCH_6",
                        "ORR_3", "MANCH_1", "MANCH_1", "MANCH_4", "MANCH_7", "MANCH_5", "MANCH_4",
                        "MANCH_5", "MANCH_2", "MANCH_2", "MANCH_1",
                        "MANCH_1", "MANCH_7", "MANCH_7", "MANCH_2", "MANCH_4", "MANCH_4", "MANCH_5",
                        "MANCH_5", "MANCH_3", "MANCH_2", "MANCH_4",
                        "MANCH_7", "MANCH_7", "MANCH_4", "MANCH_7", "ORR_5", "ORR_8", "MANCH_1", "ORR_1",
                        "ORR_4", "MANCH_3", "MANCH_1", "MANCH_3",
                        "MANCH_4", "MANCH_10", "cross", "cross", "cross", "pop inconnue", "MANCH_2", "ORR_4"))




png("tdi_all_PCA2.png", width = 580, height = 500)
par(mfrow = c(1,1))
s.class(pca.cows$li, fac.score,
        col=transp(funky(42),1),
       # grid = FALSE,
       possub = "topright",
       sub = "PCA - tdi",
        cellipse= TRUE,
        axesell=TRUE, cstar=1, cpoint=1,
        addaxes=TRUE,
       # ylim = c(-100,100)
        )

dev.off()

png("tdi_all_PCA-eig.png", width = 580, height = 500)
par(mfrow = c(1,1))
s.class(pca.cows$li, fac.score,
        col=transp(funky(41),1),
       # grid = FALSE,
       possub = "topright",
       sub = "PCA - Tsi",
        cellipse= TRUE,
        axesell=TRUE, cstar=1, cpoint=1,
        addaxes=TRUE,
       # ylim = c(-100,100)
        )

add.scatter.eig(pca.cows$eig[1:20], xax=1, yax=2,
                ratio=.2, posi="bottomleft")
dev.off()



### calculate genetic distances
library(poppr)
## I have to use a genind obj
# tdinei <- nei.dist(tdi, warning = TRUE)
## calculate euclidean distance
D <- dist(tab(tdi))
## put them in a tree
library(ape)
tre <- njs(D)
par(xpd=TRUE)

png("tdi_tree_labels.png", width = 580, height = 500)
temp <- as.integer(pop(tdi))
myCol <- transp(funky(41),.9)[temp]
plot(tre, type="unrooted", edge.w=2, font =1, show.tip.label = FALSE)

#edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))
tiplabels(pch = 19, col = myCol, adj = 0, cex = 2)
dev.off()


# Allele presence absencedata are extracted and NAs replaced usingtab:
X <- tab(tdi, NA.method="mean")
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(tdi))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]
## basic plot
png("tdi_all_PCA.png", width = 580, height = 500)
png("tdi_all_PCA_labels.png", width = 580, height = 500)
plot(pca1$li, col=myCol, cex=3, pch=myPch, xlim =c(-130, 130), ylim=c(-80, 90))
#dev.off()

## use wordcloud for non-overlapping labels
library(wordcloud)
textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=0.7, new=FALSE)
dev.off()


### check eigenvalues and other PCA parameters.
## I followed this tutorial: http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
library(factoextra)
png("tdi_all_PCAeigen.png", width = 580, height = 500)
fviz_eig(pca1)
dev.off()
```

### Working with the subset data

Exclude individuals with vcftools:
```bash
vcftools --recode --recode-INFO-all --vcf Tdi_final75.recode.vcf --remove-indv Tps_04_24  --remove-indv Tps_04_17 \
                --remove-indv Tps_04_22  --remove-indv Tps_04_23 --remove-indv Tps_04_25 --remove-indv Tps_04_30 \
                --remove-indv Tps_04_35 --remove-indv Tps_04_42 --remove-indv Tps_04_5 --remove-indv Tps_04_6 \
                --remove-indv Tpsgen.5096112 --remove-indv Tpsgen.5096121 --remove-indv Tpsgen.5096122 \
                --remove-indv Tpsgen.5096123 --remove-indv Tpsgen.5096126 --remove-indv Tpsgen.5096128 \
                --remove-indv Tpsgen.5096164 --remove-indv Tps_LMA_13 --remove-indv Tps_LMA_14 \
                --remove-indv Tps_LMA_15 --remove-indv Tps_LMA_16 --remove-indv Tps_LMA_17 \
                --remove-indv Tps_LMA_18 --remove-indv Tps_LMA_19 --remove-indv Tps_LMA_20 \
                --remove-indv Tps_LMA_21 --remove-indv Tps_LMA_22 --remove-indv Tps_LMA_23 \
                --remove-indv Tps_LMA_24 --remove-indv Tps_LMA_25 --remove-indv Tps_LMA_26 \
                --remove-indv Tps_LMA_27 --remove-indv Tps_LMA_28 --remove-indv Tps_LMA_29 \
                --remove-indv Tps_LMA_30 --remove-indv Tps_LMA_31 --remove-indv Tps_LMA_32 \
                --remove-indv Tps_LMA_33 --remove-indv Tps_LMA_34 --remove-indv Tps_LMA_35 \
                --remove-indv Tps_LMA_36 --remove-indv Tps_LMA_37 --remove-indv Tps_LMA_38 \
                --remove-indv Tps_LMA_39 --remove-indv Tps_LMA_40 --remove-indv Tps_LMA_41 \
                --remove-indv Tps_LMA_42 --remove-indv Tps_LMA_43 --remove-indv Tps_LMA_44 \
                --remove-indv Tps_LMA_45 --remove-indv Tps_LMA_46 --remove-indv Tps_LMA_47 \
                --remove-indv Tps_LMA_48 --remove-indv Tps_LMA_49 --remove-indv Tps_LMA_50 \
                --remove-indv Tps_LMA_51 --remove-indv Tps_LMA_52 --remove-indv Tps_LMA_53 \
                --remove-indv Tps_LMA_54 --remove-indv Tps_LMA_55 --remove-indv Tps_LMA_56 \
                --remove-indv Tps_LMA_57 --remove-indv Tps_LMA_58 --remove-indv Tps_LMA_59 \
                --remove-indv Tps_LMA_60 --remove-indv Tps_LMA_61 --remove-indv Tps_LMA_63 \
                --remove-indv Tps_LMA_64 --remove-indv Tps_LMA_65 --remove-indv Tps_LMA_66 \
                --remove-indv Tps_LMA_67 --remove-indv Tps_LMA_68 --remove-indv Tps_LMA_69 \
                --remove-indv Tps_LMA_70 --remove-indv Tps_LMA_71 --remove-indv Tps_LMA_72 \
                --remove-indv Tps_LMA_73 --remove-indv Tps_LMA_74 --remove-indv Tps_LMA_75 \
                --remove-indv Tps_LMA_76 --remove-indv Tps_LMA_77 --remove-indv Tps_LMA_78 \
                --remove-indv Tps_LMA_79 --remove-indv Tps_LMA_80 --remove-indv Tps_LMA_81 \
                --remove-indv Tps_LMA_82 --remove-indv Tps_LMA_83 --remove-indv Tps_LMA_84 \
                --remove-indv Tps_LMA_85 --remove-indv Tps_LMA_86 --remove-indv Tps_LMA_87 \
                --remove-indv Tps_LMA_88 --remove-indv Tps_LMA_89 --remove-indv Tps_LMA_90 \
                --remove-indv Tps_LMA_91 --remove-indv Tps_LMA_92 --remove-indv Tps_LMA_93 \
                --remove-indv Tps_LMA_94 --remove-indv Tps_LMA_95 --remove-indv Tps_LMA_96 \
                --remove-indv Tps_SRA_75 --remove-indv Tps_SRA_76 --remove-indv Tps_SRA_77 \
                --out Tdi_75_subset



```
back in R:
```R

library(vcfR)
library(adegenet)
#library(adegraphics)
#library(pegas)
#library(StAMPP)
#library(lattice)
#library(gplots)
#library(ape)
#library(ggmap)
library(poppr)

#### set working directory
setwd("/home/cravinhos/Documents/cryptic_gene_flow/tree_guillaume/3_distancetree_R/")
## 75% max missing data - subset

tdivcf <- read.vcfR("Tdi_75_subset.recode.vcf", verbose = FALSE )
tdi <- vcfR2genind(tdivcf)

## see individual names
indNames(tdi)


## set pop tdi
tdi@pop <- as.factor(c("madonna", "madonna", "Philo", "Horseranch", "Horseranch", "Horseranch", 
        "Horseranch", "FortBragg", "FortBragg", "NA", "FortBragg", "FortBragg", "Iverson", 
        "Iverson", "Iverson", "Iverson", "Orr", "FtRoss", "FtRoss", "FtRoss", "FtRoss", "Wpt128", 
        "Wpt128", "Swanton", "Swanton", "Orr", "Swanton", "Swanton", "Bigpulloutvista", 
        "Bigpulloutvista", "Bigpulloutvista", "Bigpulloutvista", "Orr", "O_11", "Summitrd", "Summitrd", 
        "Summitrd", "Summitrd", "Fishrock", "Fishrock", "Fishrock", "Fishrock", "madonna", "madonna", 
        "madonna", "M_5", "M_5", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "cross", "cross", 
        "M_4", "M_1", "M_1", "M_1", "M_3", "M_1", "M_1", "M_1", "M_5", "O_1", "M_1", "O_8", "M_2", 
        "O_1", "M_1", "O_8", "M_1", "O_1", "M_6", "M_6", "M_12", "M_3", "M_9", "M_2", "M_6", "M_12", 
        "M_1", "M_1", "M_3", "M_12", "M_1", "M_3", "M_4", "M_6", "O_3", "M_1", "M_1", "M_4", "M_7", 
        "M_5", "M_4", "M_5", "M_2", "M_2", "M_1", "M_1", "M_7", "M_7", "M_2", "M_4", "M_4", "M_5", "M_5", 
        "M_3", "M_2", "M_4", "M_7", "M_7", "M_4", "M_7", "O_5", "O_8", "M_1", "O_1", "O_4", "M_3", "M_1", 
        "M_3", "M_4", "M_10", "NA", "M_2", "O_4"))


tdi@pop

# confirm if we are doing it correctly
chr_unq <- unique(tdi@pop)
length(chr_unq)

##### PCoA - ADEGENET #######
## make PCA
x.cows <- tab(tdi, freq=TRUE, NA.method="mean")
pca.cows <- dudi.pca(x.cows, center=TRUE, scale=TRUE)

# 5  Multivariate analyses
# 5.1  Principal Component Analysis (PCA)
# Principal Component Analysis (PCA) is the amongst the most common multivariate
# analyses used in genetics.  Running a PCA ongenindobject is straightforward. 
# One needs to firstextract allelic data (as frequencies) and replace missing values 
# using the accessortabandthen use the PCA procedure (dudi.pca).  Let us use this
# approach on themicrobovdata.
# Let us first load the data:
# data(microbov)
# x.cows <- tab(microbov, freq=TRUE, NA.method="mean")
# pca.cows <- dudi.pca(x.cows, center=TRUE, scale=FALSE)
#The function dudi.pca displays a barplot of eigenvalues (thescreeplot) and asks 
# for a numberof retained principal components.  In general,  eigenvalues represent
# the amount of geneticdiversity  -  as  measured  by  the  multivariate  method  being
# used  -  represented  by  each principal component (PC). Here,  each eigenvalue is the
# variance of the corresponding PC.A sharp decrease in the eigenvalues is usually
# indicative of the boundaries between relevantstructures and random noise. 
# Here, how many axes would you retain?
# Tge = 10
# tdi = 20
# Tsi = 5

png("tdi_subset_PCA.png", width = 580, height = 500)
s.label(pca.cows$li)
s.class(pca.cows$li, fac=pop(tdi), col=funky(31))
dev.off()


# tdi
indNames(tdi)
fac.score <- factor(c("madonna", "madonna", "Philo", "Horseranch", "Horseranch", "Horseranch", 
        "Horseranch", "FortBragg", "FortBragg", "NA", "FortBragg", "FortBragg", "Iverson", 
        "Iverson", "Iverson", "Iverson", "Orr", "FtRoss", "FtRoss", "FtRoss", "FtRoss", "Wpt128", 
        "Wpt128", "Swanton", "Swanton", "Orr", "Swanton", "Swanton", "Bigpulloutvista", 
        "Bigpulloutvista", "Bigpulloutvista", "Bigpulloutvista", "Orr", "O_11", "Summitrd", "Summitrd", 
        "Summitrd", "Summitrd", "Fishrock", "Fishrock", "Fishrock", "Fishrock", "madonna", "madonna", 
        "madonna", "M_5", "M_5", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "cross", "cross", 
        "M_4", "M_1", "M_1", "M_1", "M_3", "M_1", "M_1", "M_1", "M_5", "O_1", "M_1", "O_8", "M_2", 
        "O_1", "M_1", "O_8", "M_1", "O_1", "M_6", "M_6", "M_12", "M_3", "M_9", "M_2", "M_6", "M_12", 
        "M_1", "M_1", "M_3", "M_12", "M_1", "M_3", "M_4", "M_6", "O_3", "M_1", "M_1", "M_4", "M_7", 
        "M_5", "M_4", "M_5", "M_2", "M_2", "M_1", "M_1", "M_7", "M_7", "M_2", "M_4", "M_4", "M_5", "M_5", 
        "M_3", "M_2", "M_4", "M_7", "M_7", "M_4", "M_7", "O_5", "O_8", "M_1", "O_1", "O_4", "M_3", "M_1", 
        "M_3", "M_4", "M_10", "NA", "M_2", "O_4"))




png("tdi_subset_PCA2.png", width = 580, height = 500)
par(mfrow = c(1,1))
s.class(pca.cows$li, fac.score,
        col=transp(funky(31),1),
       # grid = FALSE,
       possub = "topright",
       sub = "PCA - tdi",
        cellipse= TRUE,
        axesell=TRUE, cstar=1, cpoint=1,
        addaxes=TRUE,
        ylim = c(-60,100)#,
        #xlim = c(-30,5)
        
        )

dev.off()

png("tdi_subset_PCA-eig.png", width = 580, height = 500)
par(mfrow = c(1,1))
s.class(pca.cows$li, fac.score,
        col=transp(funky(41),1),
       # grid = FALSE,
       possub = "topright",
       sub = "PCA - Tsi",
        cellipse= TRUE,
        axesell=TRUE, cstar=1, cpoint=1,
        addaxes=TRUE,
       # ylim = c(-100,100)
        )

add.scatter.eig(pca.cows$eig[1:20], xax=1, yax=2,
                ratio=.2, posi="bottomleft")
dev.off()



### calculate genetic distances
library(poppr)
## I have to use a genind obj
# tdinei <- nei.dist(tdi, warning = TRUE)
## calculate euclidean distance
D <- dist(tab(tdi))
## put them in a tree
library(ape)
tre <- njs(D)
par(xpd=TRUE)

uColSub <- unique(myCol)
uPopSub <- unique(pop(tdi))

png("tdi_subsettree_labelsnames.png", width = 580, height = 500)
temp <- as.integer(pop(tdi))
myCol <- funky(31)[temp]
plot(tre, type="unrooted", edge.w=2, font =1, show.tip.label = FALSE)

#edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))
plot(tre, type="unrooted", edge.w=2, font =.1, show.tip.label = FALSE, cex = .5)
tiplabels(pch = 19, col = myCol, adj = 0, cex = 2)

dev.off()





## setting Manch and Orr same colour
## set pop tdi
tdi@pop <- as.factor(c("madSummit", "madSummit", "Philo", "Horseranch", "Horseranch", "Horseranch", 
        "Horseranch", "FortBragg", "FortBragg", "NA", "FortBragg", "FortBragg", "Iverson", 
        "Iverson", "Iverson", "Iverson", "Orr", "FtRoss", "FtRoss", "FtRoss", "FtRoss", "Wpt128", 
        "Wpt128", "BPVswanton", "BPVswanton", "Orr", "BPVswanton", "BPVswanton", "BPVswanton", 
        "BPVswanton", "BPVswanton", "BPVswanton", "Orr", "O", "madSummit", "madSummit", 
        "madSummit", "madSummit", "Fishrock", "Fishrock", "Fishrock", "Fishrock", "madSummit", "madSummit", 
        "madSummit", "M", "M", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "NOSIL", "cross", "cross", 
        "M", "M", "M", "M", "M", "M", "M", "M", "M", "O", "M", "O", "M", 
        "O", "M", "O", "M", "O", "M", "M", "M", "M", "M", "M", "M", "M", 
        "M", "M", "M", "M", "M", "M", "M", "M", "O", "M", "M", "M", "M", 
        "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", "M", 
        "M", "M", "M", "M", "M", "M", "M", "O", "O", "M", "O", "O", "M", "M", 
        "M", "M", "M", "NA", "M", "O"))



# Allele presence absencedata are extracted and NAs replaced usingtab:
X <- tab(tdi, NA.method="mean")
pca1 <- dudi.pca(X,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(tdi))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]
## basic plot
png("tdi_all_PCA.png", width = 580, height = 500)
png("tdi_all_PCA_labels.png", width = 580, height = 500)
plot(pca1$li, col=myCol, cex=3, pch=myPch, xlim =c(-130, 130), ylim=c(-80, 90))
#dev.off()

## use wordcloud for non-overlapping labels
library(wordcloud)
textplot(pca1$li[,1], pca1$li[,2], words=rownames(X), cex=0.7, new=FALSE)
dev.off()


### check eigenvalues and other PCA parameters.
## I followed this tutorial: http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
library(factoextra)
png("tdi_all_PCAeigen.png", width = 580, height = 500)
fviz_eig(pca1)
dev.off()
```


