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

platform=ILLUMINA

mkdir -p ${sampleid}
cd ${sampleid}


ID=${sampleid}
# N or P
FASTQ1=${fq1}
FASTQ2=${fq2}
PLATFORM=${platform}
# ILLUMINA or DNBSEQ

REF=~/data/reference/hg38/Homo_sapiens_assembly38.fasta


# -pe OpenMP 4 -l s_vmem=6G

### paired-fastq-to-unmapped-bam.wdl ###################################
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - The following metada descriptors per sample:
##  - readgroup
##  - sample_name
##  - library_name
##  - platform_unit
##  - run_date
##  - platform_name
##  - sequecing_center
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs


source ~/conda-pack/gatk4.1.9.0/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB

### Convert a pair of FASTQs to uBAM
# ReadGroup name: default A. To change, "-RG"
if [ ! -f "Unmapped.rg.bam" ]; then
	mkdir -p log
	mkdir -p FTSTmp
	/usr/bin/time -f "Memory:%M KB time:%E" -o log/FastqToSam.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 FastqToSam\
	 -F1 ${FASTQ1}\
	 -F2 ${FASTQ2}\
	 -O Unmapped.rg.bam\
	 -SM ${ID} -PL ${PLATFORM} -PU unit1 -CN OSAKA -LB Macrogen\
	 --MAX_RECORDS_IN_RAM 1000000\
	 --TMP_DIR FTSTmp\
	 > log/FastqToSam.log 2>&1
	rm -rf FTSTmp
fi


### Validates the unmapped BAM file
if [ ! -f "ValidateSamFile_uBAM.txt" ]; then
	mkdir -p log
	mkdir -p VSFTmp
	/usr/bin/time -f "Memory:%M KB time:%E" -o log/ValidateSamFile_uBAM.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 ValidateSamFile\
	 -I Unmapped.rg.bam\
	 --MODE SUMMARY\
	 -O ValidateSamFile_uBAM.txt\
	 -R ${REF}\
	 --TMP_DIR VSFTmp\
	 --MAX_RECORDS_IN_RAM 1000000\
	 > log/ValidateSamFile_uBAM.log 2>&1
	rm -rf VSFTmp
fi

