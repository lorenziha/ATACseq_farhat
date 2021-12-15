#!/usr/bin/bash

set -o errexit    # Used to exit upon error, avoiding cascading errors

Help()
{
   # Display Help
   echo "This script runs STAR on a fasta file to generate STAR index files."
   echo
   echo "Syntax: ${0} [-p|-d|-t|-R|-a|-h]"
   echo "options:"
   echo "-p     seq file prefix. [./samples.prefix]"
   echo "-d     genome dir where BOWTIE2 index files are stored. [./BOWTIE2]"
   echo "-t     Number of threads (up to 16). [16]"
   echo "-R     Reads directory. [./READS]"
   echo "-a     Annotation file in gtf format."
   echo "-h     Prints this help."
   echo
}

Log(){
	echo >> ATACseq.log 
	echo $1 >> ATACseq.log
	echo >> ATACseq.log
}

# Initialize log file
if [[ -e ATACseq.log ]]; then
	rm -f ATACseq.log
fi

# Get options
while getopts "hp:d:t:R:a:" option; do
   case $option in
        h) # display Help
                Help
                exit;;
        p) input=${OPTARG:-./samples.prefix};;
        d) GENOME_PATH=${OPTARG:-./BOWTIE2/genome.fasta};;
        t) CPU=${OPTARG:-16};;
        R) READS=${OPTARG:-./READS};;
	a) GTF_ANNOTATION=${OPTARG};;
        \?) # incorrect option
                echo
                echo "Error, Invalid option"
                echo
                Help
                exit;;
   esac
done



# Job Name
#$ -N ATACseq

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' uder current working directory (cwd)
if [ ! -d "ATACseq_output" ]; then #Create output directory in case it does NOT exist
    mkdir ATACseq_outputfi
fi
#$ -o ATACseq_output/

# Tell the job your cpu and memory requirements
#$ -pe threaded 16 
# -l mem_free=20G,h_vmem=24G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M hernan.lorenzi@nih.gov

# Make temporary directory
if [ ! -d "tmp" ]; then 
    mkdir tmp
fi

Log "Prefix file = ${input}"
Log "Reads directory = ${READS}"
Log "Annotation file = ${GTF_ANNOTATION}"
Log "Genome directory = ${GENOME_DIR}"
Log "CPU usage = ${CPU}"

################
# QC raw reads
################

fun_FASTQC(){
	module load fastqc
	prefixe_file=$1
	cpu=$2
	dir=$3
	while IFS= read -r prefix
	do
		read1=${dir}/${prefix}.R1.*fastq.gz
		read2=${dir}/${prefix}.R2.*fastq.gz
		Log "Running FASTQC on $file1 and $file2"
		Log "Starting time `date`"
		Log "fastqc -t ${cpu} ${read1} ${read2}"
		fastqc -t ${cpu} ${read1} ${read2}
		Log "Done!!! `date`"

	done < "$prefixe_file"
	module purge
}

Log "### QC raw reads ###"
fun_FASTQC ${input} ${CPU} ${READS}

# Run multiQC to merge all fastQC files together
module load multiqc
multiqc ./${READS}/
module purge

##################################
## Trimming reads with trimmomatic
## http://www.usadellab.org/cms/?page=trimmomatic
##################################

if [ ! -d "trimmomatic_output" ]; then #Create output directory in case it does NOT exist
    mkdir trimmomatic_output 
fi


Log "TRIMMING READS"


module load trimmomatic

# Genertae adapter file from trimmomatic 

cat $EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa $EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa > adapters.fasta

while IFS= read -r prefix
do
	file1=${READS}/${prefix}.R1.fastq.gz
	file2=${READS}/${prefix}.R2.fastq.gz
	outP1=./trimmomatic_output/${prefix}.R1.paired.fastq.gz
	outP2=./trimmomatic_output/${prefix}.R2.paired.fastq.gz
	outUP1=trimmomatic_output/${prefix}.R1.unpaired.fastq.gz
	outUP2=trimmomatic_output/${prefix}.R2.unpaired.fastq.gz
	log=trimmomatic_output/trim.log

	Log "Trimming reads" 
	Log "java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${CPU} -trimlog $log $file1 $file2 $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30"
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${CPU} -trimlog $log $file1 $file2 $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:adapters.fasta:2:30:10 LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
	Log "Trimmed $file 1 and $file2"
done <${input}

module purge

#############################
# Run QC on trimmed reads
#############################

fun_FASTQC ${input} ${CPU} "./trimmomatic_output"

##############################
# Run multiQC on trimmed reads
##############################

module load multiqc
multiqc ./trimmomatic_output/
module purge


##################################
# Mapped reads to reference genome
##################################

module load bowtie2/2.3.4.1 picard samtools bedtools

# Check for indexed bowtie2 library
GENOME_DIR=${GENOME_PATH%/*}
REFERENCE=${GENOME_PATH##*/}

if [[ ! -e ${GENOME_DIR}/${REFERENCE} ]]; then
	echo "ERROR, I cannot find reference genome ${GENOME_PATH}"; echo;
	exit 1
elif [[ ! -e ${GENOME_DIR}/${REFERENCE}.rev2.bt2 ]]; then
	bowtie2-build ${GENOME_DIR}/${REFERENCE} ${REFERENCE}
	mv ${REFERENCE}.*bt2 ./${GENOME_DIR}/	
fi

DB="./${GENOME_DIR}/${REFERENCE}"

while IFS= read -r prefix
do
	Log "Running bowtie2 on ${prefix}"
	Log "bowtie2 -p ${CPU} -x ${DB} -1 ${READS}/${prefix}.R1.paired.fastq.gz -2 ${READS}/${prefix}.R2.paired.fastq.gz \| samtools view -hb - \|samtools sort -@ 16 -T tmp -O BAM - \> ATACseq_output/${prefix}.sorted.bam"
	bowtie2 -p ${CPU} -x ${DB} -1 ${READS}/${prefix}.R1.paired.fastq.gz -2 ${READS}/${prefix}.R2.paired.fastq.gz | samtools view -hb - |samtools sort -@ 16 -T tmp -O BAM - > ATACseq_output/${prefix}.sorted.bam 

	Log "Removing duplicated reads with picard on ATACseq_output/${prefix}.sorted.bam"
	java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates I=ATACseq_output/${prefix}.sorted.bam O=ATACseq_output/${prefix}.sorted.dedup.bam M=${prefix}.sorted.dedup.txt READ_NAME_REGEX=null REMOVE_DUPLICATES=true

	Log "Creating bam index for ATACseq_output/${prefix}.sorted.dedup.bam"
	samtools index ATACseq_output/${prefix}.sorted.dedup.bam

	Log "Filter out blacklist regions for ATACseq_output/${prefix}.sorted.dedup.bam"
	bedtools intersect -v -abam ATACseq_output/${prefix}.sorted.dedup.bam  -b blacklist.bed > ATACseq_output/${prefix}.sorted.dedup.filter.bam

	Log "Convert bam file to bed for ATACseq_output/${prefix}.sorted.dedup.filter.bam"
	bedtools bamtobed -i ATACseq_output/${prefix}.sorted.dedup.filter.bam > ATACseq_output/${prefix}.sorted.dedup.filter.bed

	Log "Done!!"

done <$input
module purge

exit 0

# Remove reads that overlap with black list bed file
# using bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam> -v
# bedtools intersect -v -abam FILE.BAM -b BLACKLIST.BED > FILTERED.BAM 

# Running MACS2 to identify peaks
#module load macs2/2.1.0.20150731-goolf-1.7.20-Python-2.7.9

#module purge

