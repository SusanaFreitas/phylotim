The goal of this analysis is to estimate LD across chromosome, using the new genomes and Nosil samples.
The samples have very low no of reads, but lets try anyway.



- Working folder:
/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/phylotim/tree_april_21

1) Map reads
2) Sort bam files
3) Call SNPs - VarScan2 (because I am not managing with FB)
4) Get DP and choose DP filters
5) Estimate LD

To make a vcf file with the sample names (instead of sample1, sample2, etc) with mpileup:
```bash
ls *filtered.bam > order_samples
less order_samples 
sed -i 's/-filtered.bam//g' order_samples
```

```bash
screen -S will_SNPs
module load Bioinformatics/Software/vital-it
module load UHTS/Analysis/samtools/1.10
samtools mpileup -f Tps_LRv5b.fasta *filtered.bam > nosil_tps.mpileup
java -jar VarScan.v2.3.9.jar mpileup2snp nosil_Tps.mpileup  --min-reads2 5 --output-vcf 1 --vcf-sample-list order_samples > nosil_tps.vcf
```
I didnt save the names of the samples (to do so, see [here](https://www.biostars.org/p/94505/). So I need to check the order of bam files, which corresponds to the order in the vcf file.


After getting the vcf file, we will observe the no of reads per position:
```bash
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
GenomeAnalysisTK CreateSequenceDictionary --REFERENCE Tps_LRv5b.fasta
GenomeAnalysisTK VariantsToTable -R Tps_LRv5b.fasta -V nosil_tps.vcf -F CHROM -F POS -GF GT -GF DP --output Tps.nosil.DP.table
```
Got an error:
The provided VCF file is malformed at approximately line number 266519: Duplicate allele added to VariantContext: C, for input source: nosil_tps.vcf
This was because on this line:
```bash
sed -n 266519p nosil_tps.vcf
266519
```
there were more than one genotype (biallelic position)
So I removed all likes with commas:
```bash
awk -F '\t' '/^#/ {print;next;} {OFS="\t";R=$4;n=split($5,a,/[,]/);s="";for(i=1;i<=n;i++) {s=sprintf("%s%s%s%s",s,(i==1?"":","),a[i],a[i]==R?"AAAAAAAAA":"");} $5=s; print;}' nosil_tps.vcf > nosil_tps_withvirg.vcf
GenomeAnalysisTK VariantsToTable -R Tps_LRv5b.fasta -V file.vcf -F CHROM -F POS -GF GT -GF DP --output Tps.nosil.DP.table 
```
The above didn't work, so I will do like last time and just delete the lines that don't work.
```bash
sed -e '266519d' nosil_tps.vcf > file.vcf
sed -e '284207d' file.vcf > file2.vcf
GenomeAnalysisTK VariantsToTable -R Tps_LRv5b.fasta -V file2.vcf -F CHROM -F POS -GF GT -GF DP --output Tps.nosil.DP.table
```

Download to my computer:
```bash
scp sfreitas1@axiom-front1.unil.ch:/scratch/wally/FAC/FBM/DEE/tschwand/default/sfreitas/phylotim/tree_april_21/*table /home/cravinhos/Documents/other_people/will_LD
```

Download to my computer:
```bash
scp sfreitas1@axiom-front1.unil.ch:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/sample_order /home/cravinhos/Dropbox/Timema_cryptic_geneflow/2_from_reads_to_vcf/monikensis/Tms_all_varscan
scp sfreitas1@axiom-front1.unil.ch:/scratch/axiom/FAC/FBM/DEE/tschwand/default/sfreitas/geneflow/05-mapping/Tms/call_old_and_new/all_Tms.vcf /home/cravinhos/Dropbox/Timema_cryptic_geneflow/2_from_reads_to_vcf/monikensis/Tms_all_varscan
```

#### BASH ####
```bash
for ((i=3; i<=96; i +=2)); do cut -f $i,$((i+1)) Tps.nosil.DP.table  | awk '$1 != "./." {print $2}' > $i.DP; done
```

#### R ####
# define the file names list (10 samples here) 
```R
## this will make a list of samples
nameList <- c()
for (i in 3:96) { # 21 - odd number for 10 samples 
  if (i %% 2 ==1) nameList <- append(nameList, paste(i, ".DP", sep=""))
}

## make a new table with 3 columns and the number of plots we want
qlist <- matrix(nrow = 47, ncol = 3) # define number of samples (10 samples here)
## change the rownames to fit the sample order (and the nameList)
qlist <- data.frame(qlist, row.names=nameList)
## name each of the columns with the different bars
colnames(qlist)<-c('5%', '10%', '99%')

png("Tdi_Nosil.fb.DP.png", height=2000, width=1600)
par(mar=c(1, 1, 1, 1), cex=1.5, mfrow=c(12,4)) # define position for our plots
for (i in 1:48) {
  DP <- read.table(nameList[i], header = T) # DP (or number of reads) for each of the SNPs (a list of number of reads per sample)
  qlist[i,] <- quantile(DP[,1], c(.05, .1, .99), na.rm=T) # find in the DP distribution (within sample) the different quantiles
  d <- density(DP[,1], from=0, to=100, bw=1, na.rm =T) # density curve that fits the distribution of DP per sample
  plot(d, xlim = c(0,100), main=nameList[i], col="blue", xlab = dim(DP)[1], lwd=2) # on the x axis: number of reads; on the x axis: number of positions with that DP value
  abline(v=qlist[i,c(1,3)], col='red', lwd=3)
}
dev.off()

write.table(qlist, "Tps.nosildata.DP.percentiles.txt")
```
The coverage is low, I will make cutoof 8 - 80

```bash
## but first we need to create a dictionary for the fasta reference with Picard
module load UHTS/Analysis/picard-tools/2.9.0
picard-tools CreateSequenceDictionary REFERENCE=Tps_LRv5b.fasta OUTPUT=Tps_LRv5b.dict


##### - keep reads DP > 8 and DP < 80
GenomeAnalysisTK VariantFiltration -R Tps_LRv5b.fasta -V file2.vcf --genotype-filter-expression 'DP<8||DP>80' --genotype-filter-name 'DP_8-80' --set-filtered-genotype-to-no-call --output Tps.nosil_DPfilter.vcf

## OPTIONS:
## --set-filtered-genotype-to-no-call : used this to set the filtered genotypes (the ones that less than 8 > DP > 200)
## to ./.
## https://gatkforums.broadinstitute.org/gatk/discussion/7577/removing-variants-based-on-the-format-field

## I had a problem with this script and I was getting constant 'invalid argument' errors.
## I solved it by removing all spaces in between the formal arguments...
# http://selectvariants1.rssing.com/chan-9504484/all_p2.html
```



