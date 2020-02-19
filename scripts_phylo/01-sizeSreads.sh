#!/bin/sh
### script to intersect a vcf file with a bed using vcflib: vcfintersect
## usage: sbatch 04-trimming-Tge.sh
#
#SBATCH --account=tschwand_default ## account name to deduct value of run
#
#SBATCH --partition=long  ## options: normal(24h), long(10d)
#
#SBATCH --job-name=TdiTsiTRIM3
# memory
#SBATCH --ntasks=1 ## no of tasks (or threads)
#SBATCH --cpus-per-task=1 ## no of cores used per thread
#SBATCH --mem-per-cpu=60G ## memory used per cpu (or core) in Mb
# #SBATCH --time=23:00:00 # for this option I will use the default, which is 5 days
# outputs
#SBATCH --output=TgeTRIM3.out
#SBATCH --error=TgeTRIM3.err
#SBATCH --mail-user=susana.freitas@unil.ch
#SBATCH --mail-type=ALL


module add Bioinformatics/Software/vital-it
module load UHTS/Quality_control/cutadapt/2.3


for dir in */; do
echo $dir
cd $dir;
	for file in *.fastq.gz; do
        	cutadapt -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 60 -o $file-size60.fq.gz $file;
	done;
cd ../
done

module rm UHTS/Quality_control/cutadapt/2.3

# To trim a 3’ adapter, the basic command-line for Cutadapt is:
# cutadapt -a AACCGGTT -o output.fastq input.fastq


# the original Illumina adapter is : AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC.
# However, I got this message for ALL reads in the output file of the 
# first Cutadapt run: 'The adapter is preceded by "A" extremely often'
# To fix the problem I added "A" to the beginning of the adapter sequence and re-ran Cutadapt.
# Also, I eliminated all sequences smaller than 80 bp: -m 80

