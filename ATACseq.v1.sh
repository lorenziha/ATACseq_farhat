#!/usr/bin/bash

Help()
{
   # Display Help
   echo "This script runs STAR on a fasta file to generate STAR index files."
   echo
   echo "Syntax: ${0} [-p|-d|-t|-R|-a|-h]"
   echo "options:"
   echo "-p     seq file prefix. [./samples.prefix]"
   echo "-d     genome dir where STAR index files are stored. [./STAR_DB]"
   echo "-t     Number of threads (up to 16). [16]"
   echo "-R     Reads directory. [./READS]"
   echo "-a     Annotation file in gtf format. [./gencode.v38.primary_assembly.annotation.gtf.gz]"
   echo "-h     Prints this help."
   echo
}

while getopts "hp:d:t:R:a:" option; do
   case $option in
        h) # display Help
                Help
                exit;;
        p) input=${OPTARG};;
        d) GENOME_DIR=${OPTARG};;
        t) CPU=${OPTARG};;
        R) READS=${OPTARG};;
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
#$ -N RNAseq

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' uder current working directory (cwd)
# - if [ ! -d "RNASEQ_output" ]; then #Create output directory in case it does NOT exist
# -     mkdir RNASEQ_output
# - fi
# -$ -o RNASEQ_output/

# Tell the job your cpu and memory requirements
#$ -pe threaded 16 
# -l mem_free=20G,h_vmem=24G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M hernan.lorenzi@nih.gov


# Add in the line below the path to the file containing raw read file prefixes; default "./samples.prefix".
if [ -z ${input+x} ]; 
	then 
	echo "Using prefix file default ./samples.prefix"
	input=./samples.prefix
	else 
	echo Using prefix file ${input}
fi

# Modify the path below to point to the directory where yoiu store your raw read files. Subsequent processing read files and QC files will be stored in the same directory; default "./READS"
if [ -z ${READS+x} ];
        then
        echo "Using reads directory default ./READS"
        READS=./READS
        else
        echo Using prefix file ${READS}
fi

# Modify the line below to point to the Reference annotation file in GTF format; default human chromosome 1 ONLY "./REFERENCE/grch38/Chr1.GRCh38.103.gtf"
if [ -z ${GTF_ANNOTATION+x} ]; 
        then 
        echo "Using annotation file default ./gencode.v38.primary_assembly.annotation.gtf.gz"
        GTF_ANNOTATION=./gencode.v38.primary_assembly.annotation.gtf.gz
        else 
        echo Using prefix file ${GTF_ANNOTATION}
fi

# Modify the path below to point at the directory where you store the reference genome already indexed for STAR (you can get the indexed reference for human genome from STAR website)
if [ -z ${GENOME_DIR+x} ]; 
        then 
        echo "Using STAR DB directory default ./STAR_DB"
        GENOME_DIR=./STAR_DB
        else 
        echo Using prefix file ${GENOME_DIR}
fi

# Enter CPU usage, default = 8
if [ -z ${CPU+x} ];
        then
        echo "Using CPU usage default = 16"
        CPU=16
        else
        echo Using CPU usage = ${CPU}
fi



# QC raw reads

module load fastqc
while IFS= read -r prefix
do
	file1=${READS}/${prefix}.R1.fastq.gz
	file2=${READS}/${prefix}.R2.fastq.gz
	echo
	echo Running FASTQC on $file1 and $file2
	echo Starting time `date`
	echo fastqc -t ${CPU} $file1 $file2
	
	fastqc -t ${CPU} $file1 $file2
	echo
	echo Done!!! `date`
	echo
done < "$input"
module unload fastqc

# Run multiQC to merge all fastQC files together

module load multiqc
multiqc ${READS}/
module unload multiqc

## Trimming reads with trimmomatic
## http://www.usadellab.org/cms/?page=trimmomatic

echo
echo TRIMMING READS
echo

module load trimmomatic

while IFS= read -r prefix
do
	file1=${READS}/${prefix}.R1.fastq.gz
	file2=${READS}/${prefix}.R2.fastq.gz
	outP1=${READS}/${prefix}.R1.paired.fastq.gz
	outP2=${READS}/${prefix}.R2.paired.fastq.gz
	outUP1=${READS}/${prefix}.R1.unpaired.fastq.gz
	outUP2=${READS}/${prefix}.R2.unpaired.fastq.gz
	log=${READS}/${line}.trim.log

	echo Trimming reads 
	echo java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${CPU} -trimlog $log $file1 $file2 $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
	echo
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads ${CPU} -trimlog $log $file1 $file2 $outP1 $outUP1 $outP2 $outUP2 ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10 LEADING:10 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
	echo
done <$input
module unload trimmomatic
module unload Java

module load fastqc
while IFS= read -r prefix
do
        file1=${READS}/${prefix}.R1.paired.fastq.gz
        file2=${READS}/${prefix}.R2.paired.fastq.gz
        echo
        echo Running FASTQC on $file1 and $file2
        echo Starting time `date`
	echo fastqc -t ${CPU} $file1 $file2

        fastqc -t ${CPU} $file1 $file2
        echo
        echo Done!!! `date`
        echo
done < "$input"
module purge


# Run STAR to map RNAseq reads to reference
module load STAR

while IFS= read -r prefix
do
	echo
	echo running STAR on ${prefix}
	echo
	echo STAR --runMode alignReads --runThreadN ${CPU} --genomeDir ${GENOME_DIR} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ${READS}/${prefix}.R1.paired.fastq.gz ${READS}/${prefix}.R2.paired.fastq.gz
echo 

	STAR --runMode alignReads --runThreadN ${CPU} --genomeDir ${GENOME_DIR} --readFilesCommand zcat --outFileNamePrefix ${prefix}. --outSAMtype BAM SortedByCoordinate --readFilesIn ${READS}/${prefix}.R1.paired.fastq.gz ${READS}/${prefix}.R2.paired.fastq.gz
	mv ${prefix}.Aligned.sortedByCoord.out.bam ${prefix}.sorted.bam
	echo
	echo Done!!
	echo
done <$input
module purge

# Sort the resulting sam file by chromosome position with samtools; remove duplicated reads with picard; 

module load samtools
module load picard

files=()
while IFS= read -r prefix 
do
	# No sorting needed. Sorting by coordinate is done during read mapping with STAR

	echo REMOVE DUPLICATES
	echo java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.sorted.dedup.bam M=${prefix}.sorted.dedup.txt READ_NAME_REGEX=null REMOVE_DUPLICATES=true
	java -jar ${EBROOTPICARD}/picard.jar MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.sorted.dedup.bam M=${prefix}.sorted.dedup.txt READ_NAME_REGEX=null REMOVE_DUPLICATES=true
	
	echo GENERATING BAM INDEX
	echo samtools index ${prefix}.sorted.dedup.bam
	samtools index ${prefix}.sorted.dedup.bam
	
	# Storing BAM file names for counting reads in bulk below
	files+=(${prefix}.sorted.dedup.bam)
        echo
        echo Done!!
        echo
done <$input
module purge 

# Counting reads per gene per sample with featureCounts  and store results in file "all.counts"
# Duplicated reads are removed during picard MarkDuplicates
module load subread
echo
echo COUNTING READS FROM THE FOLLOWING sorted.dedup.bam FILES:
echo ${files[@]}
echo featureCounts -p -a $GTF_ANNOTATION -T ${CPU} -s 0 -F GTF -t exon -g gene_id -o all_counts.txt -R BAM --extraAttributes gene_name ${files[@]}
# Counting PE reads as a single fragment (-p option)
featureCounts -p -a $GTF_ANNOTATION -T ${CPU} -s 0 -F GTF -t exon -g gene_id -o all_counts.txt -R BAM --extraAttributes gene_name -O ${files[@]} 

module purge
