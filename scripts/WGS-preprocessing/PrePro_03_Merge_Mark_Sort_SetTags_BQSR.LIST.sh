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



ID=${sampleid}
# N or P

REF=~/data/reference/hg38/Homo_sapiens_assembly38.fasta


#-pe OpenMP 4 -l s_vmem=6-10G

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
command="bwa mem -K 100000000 -p -v 3 -Y -M ${REF}"


source ~/conda-pack/gatk4.1.9.0/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB

# #BWA-MEM Aligned BAM query-grouped -> queryname-sorted (To avoid "Too many open files @ MergeBamAlignment")
# if [ ! -f "Aligned.qs.bam" ]; then
#   mkdir -p SortqsTmp
#   /usr/bin/time -f "Memory:%M KB time:%E" -o SortSam-queryname.txt \
#   gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
#    SortSam\
#    -I Aligned.bam\
#    -O Aligned.qs.bam\
#    --SORT_ORDER "queryname"\
#    --TMP_DIR ./SortqsTmp\
#    --MAX_RECORDS_IN_RAM 1000000\
#    > sortsam-queryname.log 2>&1
# else
#   echo "SortSam-queryname Already Done"
# fi


### Merge original uBAM and BWA-aligned BAM 
if [ ! -f "${ID}-hg38.bam" ]; then
if [ ! -f "recal_data.table" ]; then
if [ ! -f "Aligned.rg.marked.cs.fixed.bam" ]; then
if [ ! -f "Aligned.rg.marked.bam" ]; then
if [ ! -f "Aligned.rg.bam" ]; then
  mkdir -p log
  mkdir -p MergeTmp
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/MergeBamAlignment.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   MergeBamAlignment\
   --VALIDATION_STRINGENCY SILENT\
   --EXPECTED_ORIENTATIONS FR\
   --ATTRIBUTES_TO_RETAIN X0\
   --ALIGNED_BAM Aligned.bam\
   --UNMAPPED_BAM Unmapped.rg.bam\
   --OUTPUT Aligned.rg.bam\
   --REFERENCE_SEQUENCE ${REF}\
   --SORT_ORDER "unsorted"\
   --IS_BISULFITE_SEQUENCE false\
   --ALIGNED_READS_ONLY false\
   --CLIP_ADAPTERS false\
   --MAX_RECORDS_IN_RAM 2000000\
   --ADD_MATE_CIGAR true\
   --MAX_INSERTIONS_OR_DELETIONS -1\
   --PRIMARY_ALIGNMENT_STRATEGY MostDistant\
   --PROGRAM_RECORD_ID "bwamem"\
   --PROGRAM_GROUP_VERSION "${version}"\
   --PROGRAM_GROUP_COMMAND_LINE "${command}"\
   --PROGRAM_GROUP_NAME "bwamem"\
   --UNMAPPED_READ_STRATEGY COPY_TO_TAG\
   --ALIGNER_PROPER_PAIR_FLAGS true\
   --UNMAP_CONTAMINANT_READS true\
   --TMP_DIR MergeTmp\
   > log/MergeBamAlignment.log 2>&1
  rm -rf MergeTmp
fi

if [ -f "Aligned.rg.bam" ]; then
  if [ -f "Unmapped.rg.bam" ]; then
    rm Unmapped.rg.bam
  fi
  if [ -f "Aligned.bam" ]; then
    rm Aligned.bam
  fi
fi
fi
fi
fi
fi



### Mark duplicate reads to avoid counting non-independent observations
# Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
# This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
# While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
if [ ! -f "${ID}-hg38.bam" ]; then
if [ ! -f "recal_data.table" ]; then
if [ ! -f "Aligned.rg.marked.cs.fixed.bam" ]; then
if [ ! -f "Aligned.rg.marked.bam" ]; then
  mkdir -p log
  mkdir -p MarkTmp
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/MarkDuplicates.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   MarkDuplicates\
   -I Aligned.rg.bam\
   -O Aligned.rg.marked.bam\
   -M marked_dup_metrics.txt\
   --ASSUME_SORT_ORDER "queryname"\
   --VALIDATION_STRINGENCY SILENT\
   --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500\
   --CREATE_MD5_FILE true\
   --TMP_DIR ./MarkTmp\
   --MAX_RECORDS_IN_RAM 1000000\
   > log/MarkDuplicates.log 2>&1
  rm -rf MarkTmp
fi

if [ -f "Aligned.rg.marked.bam" ]; then
  if [ -f "Aligned.rg.bam" ]; then
    rm Aligned.rg.bam
  fi
fi
fi
fi
fi


### Sort BAM file by coordinate order and fix tag values for NM and UQ
if [ ! -f "${ID}-hg38.bam" ]; then
if [ ! -f "recal_data.table" ]; then
if [ ! -f "Aligned.rg.marked.cs.fixed.bam" ]; then
  mkdir -p log
  mkdir SortcsTmp
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/SortSam-coordinate.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   SortSam\
   -I Aligned.rg.marked.bam\
   -O /dev/stdout\
   --SORT_ORDER "coordinate"\
   --CREATE_INDEX false\
   --CREATE_MD5_FILE false\
   --TMP_DIR SortcsTmp\
   --MAX_RECORDS_IN_RAM 1000000\
   2> log/SortSam-coordinate.log\
   |\
   /usr/bin/time -f "Memory:%M KB time:%E" -o log/SetNmMdAndUqTags.txt \
   gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
    SetNmMdAndUqTags\
    -R ${REF}\
    -I /dev/stdin\
    -O Aligned.rg.marked.cs.fixed.bam\
    --CREATE_INDEX true\
    --CREATE_MD5_FILE true\
    > log/SetNmMdAndUqTags.log 2>&1
  rm -rf SortcsTmp
fi

if [ -f "Aligned.rg.marked.cs.fixed.bam" ]; then
  if [ -f "Aligned.rg.marked.bam" ]; then
    rm Aligned.rg.marked.bam
  fi
  if [ -f "Aligned.rg.marked.bam.md5" ]; then
    rm Aligned.rg.marked.bam.md5
  fi
fi
fi
fi



### Generate Base Quality Score Recalibration (BQSR) model
if [ ! -f "${ID}-hg38.bam" ]; then
if [ ! -f "recal_data.table" ]; then
  mkdir -p log
  mkdir -p BRTmp
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/BaseRecalibrator.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   BaseRecalibrator\
   -I Aligned.rg.marked.cs.fixed.bam\
   -R ${REF}\
   --use-original-qualities\
   --known-sites ~/data/reference/SNP/GCF_000001405.38.dbSNP153.GRCh38p12.GATK.vcf.gz\
   --known-sites ~/data/reference/INDEL/Homo_sapiens_assembly38.known_indels.vcf\
   --known-sites ~/data/reference/INDEL/Mills_and_1000G_gold_standard.indels.hg38.vcf\
   -O recal_data.table\
   --tmp-dir BRTmp\
   > log/BaseRecalibrator.log 2>&1
  rm -rf BRTmp
fi
fi

### Apply the recalibration model
if [ ! -f "${ID}-hg38.bam" ]; then
  mkdir -p log
  mkdir -p ApplyTmp
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/ApplyBQSR.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   ApplyBQSR\
   -I Aligned.rg.marked.cs.fixed.bam\
   -R ${REF}\
   --bqsr-recal-file recal_data.table\
   -O ${ID}-hg38.bam\
   --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30\
   --add-output-sam-program-record\
   --create-output-bam-md5\
   --use-original-qualities\
   --tmp-dir ApplyTmp\
   > log/ApplyBQSR.log 2>&1
  rm -rf ApplyTmp
fi

if [ -f "${ID}-hg38.bam" ]; then
  if [ -f "Aligned.rg.marked.cs.fixed.bam" ]; then
    rm Aligned.rg.marked.cs.fixed.bam
  fi
  if [ -f "Aligned.rg.marked.cs.fixed.bai" ]; then
    rm Aligned.rg.marked.cs.fixed.bai
  fi
  if [ -f "Aligned.rg.marked.cs.fixed.bam.md5" ]; then
    rm Aligned.rg.marked.cs.fixed.bam.md5
  fi
fi
