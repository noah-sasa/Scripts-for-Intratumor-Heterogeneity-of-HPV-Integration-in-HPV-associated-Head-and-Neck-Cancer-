#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V



ID=$1

#-pe OpenMP 2 -l s_vmem=6G

source /work23/home/nsasa/conda-pack/gatk/bin/activate

mkdir -p PrePro_00_log

#convert MGI header to ILLIMINA header
/usr/bin/time -f "Memory:%M KB time:%E" -o PrePro_00_log/${ID}_R1_convertHeader.txt \
python /work23/home/nsasa/tools/convertHeaders_5elements_LessThan90Letters.py \
	-i ${ID}_R1.fastq.gz \
	-o ${ID}_ILLUMINA_R1.fastq.gz \
	> PrePro_00_log/${ID}_R1_convertHeader.log 2>&1


source /work23/home/nsasa/conda-pack/gatk/bin/deactivate

