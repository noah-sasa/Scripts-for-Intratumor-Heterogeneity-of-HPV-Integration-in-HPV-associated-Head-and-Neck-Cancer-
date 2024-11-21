#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V



INPUT_BAM=$1
ID=$2
#ID+N/P

#-pe OpenMP 2 -l s_vmem=6G

source ~/conda-pack/gatk/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB


#Varidates a SAM or BAM file
/usr/bin/time -f "Memory:%M KB time:%E" -o ${ID}_ValidateSamFile.txt \
gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
 ValidateSamFile\
 -I ${INPUT_BAM}\
 --MODE SUMMARY\
 > ${ID}_validatesamfile.log 2>&1



