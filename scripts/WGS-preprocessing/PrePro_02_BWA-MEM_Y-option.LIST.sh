#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V

list=$1
sampleinfo=$(cat ${list} | head -n ${SGE_TASK_ID} | tail -n 1)
sampleid=$(echo $sampleinfo | awk '{print $1}')
fq1=$(echo $sampleinfo | awk '{print $2}')
fq2=$(echo $sampleinfo | awk '{print $3}')
fqs=$(echo $sampleinfo | awk '{print $4}')

cpu=32

mkdir -p ${sampleid}
cd ${sampleid}


CPU=${cpu}

REF=~/data/reference/hg38/Homo_sapiens_assembly38.fasta


# -pe OpenMP 8~16 -l s_vmem=6G

### processing-for-variant-discovery-gatk4.wdl ###################################
## This WDL pipeline implements data pre-processing according to the GATK Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile 
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
## Software version requirements 
## - GATK 4 or later
## - BWA 0.7.15-r1140
## - Picard 2.16.0-SNAPSHOT
## - Samtools 1.3.1 (using htslib 1.3.1)
## - Python 2.7




source ~/conda-pack/gatk/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB


### Get version of BWA
version=$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //')


### Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
if [ ! -f "Aligned.bam" ]; then
	mkdir -p log
	/usr/bin/time -f "Memory:%M KB time:%E" -o log/SamToFastq.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 SamToFastq\
	 -I Unmapped.rg.bam\
	 -F /dev/stdout\
	 --INTERLEAVE true\
	 --INCLUDE_NON_PF_READS true\
	 2> log/SamToFastq.log\
	 |\
	 /usr/bin/time -f "Memory:%M KB time:%E" -o log/bwamem.txt \
	 bwa mem -K 100000000 -p -v 3 -t ${CPU} -Y ${REF} /dev/stdin -  2> >(tee log/bwa.stderr.log >&2)\
	 |\
	 samtools view -1 - > Aligned.bam
fi	

