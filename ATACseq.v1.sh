#!/usr/bin/bash

set -o errexit    # Used to exit upon error, avoiding cascading errors

Help()
{
   # Display Help
   echo "This script runs the ATACseq pipeline using the following approach:"
   echo "cutadapt (adapter removal and quality trimming) > bowtie2 (mapping) > bedtools (remove blacklist and mitochondrial reads) > macs2 (peak calling)"
   echo
   echo "Syntax: ${0} [-p|-d|-t|-R|-o|-h]"
   echo "options:"
   echo "-p     seq file prefix. [./samples.prefix]"
   echo "-d     genome dir where BOWTIE2 index files are stored. [./BOWTIE2]"
   echo "-t     Number of threads (up to 16). [16]"
   echo "-R     Reads directory. [./READS]"
   echo "-o     Output directory [./ATACseq_output]"
   echo "-h     Prints this help."
   echo
}

# Get options
while getopts "hp:d:t:R:o:" option; do
   case $option in
        h) # display Help
                Help
                exit;;
        p) input=${OPTARG:-./samples.prefix};;
        d) GENOME_PATH=${OPTARG:-./BOWTIE2/genome.fasta};;
        t) CPU=${OPTARG:-16};;
        R) READS=${OPTARG:-./READS};;
	o) WORKING_DIR=${OPTARG:-./ATACseq_output};;
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
if [ ! -d "${WORKING_DIR}" ]; then #Create output directory in case it does NOT exist
    mkdir ${WORKING_DIR}
fi
#$ -o ${WORKING_DIR}/

# Tell the job your cpu and memory requirements
#$ -pe threaded 16 
# -l mem_free=20G,h_vmem=24G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M hernan.lorenzi@nih.gov

Log(){
        echo >> ${WORKING_DIR}/ATACseq.log
        echo "$1" >> ${WORKING_DIR}/ATACseq.log
        echo >> ${WORKING_DIR}/ATACseq.log
}

# Initialize log file
if [[ -e ATACseq.log ]]; then
	rm -f ATACseq.log
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

Log "echo ### QC raw reads ###"
fun_FASTQC ${input} ${CPU} ${READS}

# Run multiQC to merge all fastQC files together
module load multiqc
multiqc --outdir ${WORKING_DIR}/multiqc ./${READS}/
module purge

##################################
## Trimming reads with cutadapt
## http://www.usadellab.org/cms/?page=cutadapt
##################################
echo WORKING_DIR = ${WORKING_DIR} .
if [ ! -d "${WORKING_DIR}/cutadapt_output" ]; then #Create output directory in case it does NOT exist
    mkdir ${WORKING_DIR}/cutadapt_output 
fi


Log "TRIMMING READS"


module load cutadapt/2.7-Python-3.7.3

while IFS= read -r prefix
do

	file1=${READS}/${prefix}.R1.fastq.gz
	file2=${READS}/${prefix}.R2.fastq.gz
	outP1=${WORKING_DIR}/cutadapt_output/${prefix}.R1.paired.fastq.gz
	outP2=${WORKING_DIR}/cutadapt_output/${prefix}.R2.paired.fastq.gz

	Log "Trimming reads"
        Log "cutadapt -u 15 -U 15 --minimum-length=3 --quality-base=33 --error-rate=0.1 -g CTGTCTCTTATACACATCT -G CTGTCTCTTATACACATCT -o $outP1  -p $outP2 $file1 $file2"	
	cutadapt -u 15 -U 15 --minimum-length=3 --quality-base=33 --error-rate=0.1 \
		-a CTGTCTCTTATACACATCT \
		-A CTGTCTCTTATACACATCT \
		-o $outP1  -p $outP2 $file1 $file2
	Log "Trimmed $file 1 and $file2"
done <${input}

module purge


#############################
# Run QC on trimmed reads
#############################

fun_FASTQC ${input} ${CPU} "${WORKING_DIR}/cutadapt_output"

##############################
# Run multiQC on trimmed reads
##############################

module load multiqc
multiqc --outdir ${WORKING_DIR}/multiqc ${WORKING_DIR}/cutadapt_output/
module purge

##################################
Log "Map reads to reference genome"
##################################

module load bowtie2/2.3.4.1 picard samtools bedtools macs2

# Check for indexed bowtie2 library
GENOME_DIR=${GENOME_PATH%/*}
REFERENCE=${GENOME_PATH##*/}

if [[ ! -e ${GENOME_DIR}/${REFERENCE} ]]; then
	echo "ERROR, I cannot find reference genome ${GENOME_PATH}"; echo;
	exit 1
elif [[ ! -e ${GENOME_DIR}/${REFERENCE}.rev.2.bt2 ]]; then
	bowtie2-build ${GENOME_DIR}/${REFERENCE} ${REFERENCE}
	mv ${REFERENCE}*bt2 ./${GENOME_DIR}/	
fi

DB="./${GENOME_DIR}/${REFERENCE}"

while IFS= read -r prefix
do
	Log "Running bowtie2 on ${prefix}"
	Log "bowtie2 --rg "ID:${prefix}\tPL:Illumina\tLB:${prefix}\tPU:${prefix}\tSM:${prefix}" -p ${CPU} -x ${DB} -1 ${WORKING_DIR}/cutadapt_output/${prefix}.R1.paired.fastq.gz -2 ${WORKING_DIR}/cutadapt_output/${prefix}.R2.paired.fastq.gz \| samtools view -hb - \|samtools sort -@ 16 -T tmp -O BAM - \> ${WORKING_DIR}/${prefix}.sorted.bam"
	bowtie2 --rg "ID:${prefix}	PL:Illumina	SM:${prefix}	LB:${prefix}" -p ${CPU} -x ${DB} -1 ${WORKING_DIR}/cutadapt_output/${prefix}.R1.paired.fastq.gz -2 ${WORKING_DIR}/cutadapt_output/${prefix}.R2.paired.fastq.gz | samtools view -hb - |samtools sort -@ 16 -T tmp -O BAM - > ${WORKING_DIR}/${prefix}.sorted.bam 

	Log "java -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics R=./GRCh38/GRCh38.primary_assembly.genome.fa I=./${WORKING_DIR}/${prefix}.sorted.bam O=${prefix}.metrics.txt"
	java -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics R=./GRCh38/GRCh38.primary_assembly.genome.fa I=./${WORKING_DIR}/${prefix}.sorted.bam O=${WORKING_DIR}/${prefix}.metrics.txt

	Log "Removing duplicated reads with picard on ${WORKING_DIR}/${prefix}.sorted.bam"
	java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates I=${WORKING_DIR}/${prefix}.sorted.bam O=${WORKING_DIR}/${prefix}.sorted.dedup.bam M=${prefix}.sorted.dedup.txt READ_NAME_REGEX=null REMOVE_DUPLICATES=true

	Log "Creating bam index for ${WORKING_DIR}/${prefix}.sorted.dedup.bam"
	samtools index ${WORKING_DIR}/${prefix}.sorted.dedup.bam
	mv ${prefix}.sorted.dedup.txt ${WORKING_DIR}/

	Log "Filter out blacklist and mitochondrial egions from ${WORKING_DIR}/${prefix}.sorted.dedup.bam"
	bedtools intersect -v -abam ${WORKING_DIR}/${prefix}.sorted.dedup.bam  -b ./blacklist/blacklist_and_mito.bed > ${WORKING_DIR}/${prefix}.sorted.dedup.filter.bam

	Log "Convert bam file to bed for ${WORKING_DIR}/${prefix}.sorted.dedup.filter.bam"
	bedtools bamtobed -i ${WORKING_DIR}/${prefix}.sorted.dedup.filter.bam > ${WORKING_DIR}/${prefix}.sorted.dedup.filter.bed

	Log "Calling peaks with macs2 on ${WORKING_DIR}/${prefix}.sorted.dedup.filter.bed"
	macs2 callpeak --treatment ${WORKING_DIR}/${prefix}.sorted.dedup.filter.bed \
		        --name ${prefix} --format BED \
			--gsize hs --shift -75 --extsize 150 \
			--nomodel --nolambda -B --SPMR --keep-dup all \
			-p 0.05 --verbose 3 --outdir ${WORKING_DIR}/macs2_output

	Log "Done!!"

done <$input
module purge

# Generate metrics summary
grep -h '^CATEGORY' ${WORKING_DIR}/*.metrics.txt |head -1 | sed -e 's/CATEGORY/SAMPLE_ID	CATEGORY/'> ${WORKING_DIR}/read_mapping_metrics_summary.txt
while IFS= read -r prefix
do
	grep -hv '^#\|^CATEGORY' ${WORKING_DIR}/${prefix}.metrics.txt| grep . |sed 's/^/'${prefix}'	/' >> ${WORKING_DIR}/read_mapping_metrics_summary.txt

done <$input

exit 0

