#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V

list=$1
sampleinfo=$(cat ${list} | head -n ${SGE_TASK_ID} | tail -n 1)
sampleid=$(echo $sampleinfo | awk '{print $1}')
fq1=$(echo $sampleinfo | awk '{print $2}')
fq2=$(echo $sampleinfo | awk '{print $3}')
fqs=$(echo $sampleinfo | awk '{print $4}')

postqc_dir=$2
# ~/data/PCAWG_other/GATK4_pre-processing_bwamem-y/PostQC

mkdir -p ${sampleid}
cd ${sampleid}


ID=${sampleid}
#ID+N/P
INPUT_BAM=${ID}-hg38.bam
POSTQC_DIR=${postqc_dir}
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis/WGS_Data/PostQC

#-pe OpenMP 4 -l s_vmem=6G

source ~/conda-pack/gatk/bin/activate

mkdir -p PostAlignmentQC.log

#Samtools flagstat
mkdir -p ${POSTQC_DIR}/Samtools_flagstat
/usr/bin/time -f "Memory:%M KB time:%E" -o PostAlignmentQC.log/Samtools_flagstat.txt \
samtools flagstat -@ 4 ${INPUT_BAM} > ${POSTQC_DIR}/Samtools_flagstat/${ID}_flagstat.txt



#FastQC
source ~/conda-pack/RNA-seq/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB


mkdir -p FastQCTmp
mkdir -p ${POSTQC_DIR}/FastQC
/usr/bin/time -f "Memory:%M KB time:%E" -o PostAlignmentQC.log/FastQC.txt \
fastqc --nogroup -t 4 -d FastQCTmp -o ${POSTQC_DIR}/FastQC ${INPUT_BAM} > PostAlignmentQC.log/FastQC.log 2>&1

rm -rf FastQCTmp

#MultiQC
# cd ${FASTQC_DIR}
# multiqc .




