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
read_length=$(echo $sampleinfo | awk '{print $5}')

postqc_dir=$2
# ~/data/PCAWG_other/GATK4_pre-processing_bwamem-y/PostQC
#read_length=$3
# 151

mkdir -p ${sampleid}
cd ${sampleid}


ID=${sampleid}
#ID+N/P
INPUT_BAM=${ID}-hg38.bam
POSTQC_DIR=${postqc_dir}
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis/WGS_Data/PostQC
READ_LENGTH=${read_length}
REF=~/data/reference/hg38/Homo_sapiens_assembly38.fasta
INTERVALS=~/data/reference/hg38/Homo_sapiens_assembly38_Autosomal.interval_list
DB_SNP=~/data/reference/SNP/GCF_000001405.38.dbSNP153.GRCh38p12.GATK.vcf.gz

#-pe OpenMP 4 -l s_vmem=6G

source ~/conda-pack/gatk/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB

mkdir -p PostAlignmentQC.log
mkdir -p ${POSTQC_DIR}/gatk4

#Varidates a SAM or BAM file
if [ ! -f "PostAlignmentQC.log/ValidateSamFile.txt" ]; then
	mkdir -p VSFTmp
	/usr/bin/time -f "Memory:%M KB time:%E" -o PostAlignmentQC.log/ValidateSamFile.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 ValidateSamFile\
	 -I ${INPUT_BAM}\
	 --MODE SUMMARY\
	 -O ${POSTQC_DIR}/gatk4/${ID}_ValidateSamFile.txt\
	 -R ${REF}\
	 --TMP_DIR VSFTmp\
	 --MAX_RECORDS_IN_RAM 1000000\
	 > PostAlignmentQC.log/ValidateSamFile.log 2>&1
	echo "ValidateSamFile Done"
	rm -rf VSFTmp
else
	echo "ValidateSamFile Already Done"
fi


#gatk4 CollectMultipleMetrics
if [ ! -f "PostAlignmentQC.log/CollectMultipleMetrics.txt" ]; then
	mkdir -p CMMTmp
	/usr/bin/time -f "Memory:%M KB time:%E" -o PostAlignmentQC.log/CollectMultipleMetrics.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 CollectMultipleMetrics\
	 -R ${REF}\
	 -I ${INPUT_BAM}\
	 -O ${POSTQC_DIR}/gatk4/${ID}_multiple_metrics\
	 --DB_SNP ${DB_SNP}\
	 --TMP_DIR CMMTmp\
	 --MAX_RECORDS_IN_RAM 1000000\
	 > PostAlignmentQC.log/CollectMultipleMetrics.log 2>&1
	echo "CollectMultipleMetrics Done"
	rm -rf CMMTmp
else
	echo "CollectMultipleMetrics Already Done"
fi

#gatk4 CollectWgsMetrics
#The CollectRawWgsMetrics have base and mapping quality score thresholds set '3' and '0' respectively (less stringent),
# while the CollectWgsMetrics tool has the default threshold values set to '20'.
if [ ! -f "PostAlignmentQC.log/CollectWgsMetrics.txt" ]; then
	mkdir -p CWMTmp
	/usr/bin/time -f "Memory:%M KB time:%E" -o PostAlignmentQC.log/CollectWgsMetrics.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 CollectWgsMetrics\
	 -R ${REF}\
	 -I ${INPUT_BAM}\
	 -O ${POSTQC_DIR}/gatk4/${ID}_WgsMetrics\
	 --TMP_DIR CWMTmp\
	 --MAX_RECORDS_IN_RAM 1000000\
	 --INTERVALS ${INTERVALS}\
	 --READ_LENGTH ${READ_LENGTH}\
	 > PostAlignmentQC.log/CollectWgsMetrics.log 2>&1
	echo "CollectWgsMetrics Done"
	rm -rf CWMTmp
else
	echo "CollectWgsMetrics Already Done"
fi

#gatk4 CollectGcBiasMetrics
if [ ! -f "PostAlignmentQC.log/CollectGcBiasMetrics.txt" ]; then
	mkdir -p CGCBMTmp
	/usr/bin/time -f "Memory:%M KB time:%E" -o PostAlignmentQC.log/CollectGcBiasMetrics.txt \
	gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
	 CollectGcBiasMetrics\
	 -R ${REF}\
	 -I ${INPUT_BAM}\
	 -O ${POSTQC_DIR}/gatk4/${ID}_gc_bias_metrics\
	 -CHART ${POSTQC_DIR}/gatk4/${ID}_gc_bias_metrics.pdf\
	 -S ${POSTQC_DIR}/gatk4/${ID}_gc_bias_summary.txt\
	 --TMP_DIR CGCBMTmp\
	 --MAX_RECORDS_IN_RAM 1000000\
	 > PostAlignmentQC.log/CollectGcBiasMetrics.log 2>&1
	echo "CollectGcBiasMetrics Done"
	rm -rf CGCBMTmp
else
	echo "CollectGcBiasMetrics Already Done"
fi


