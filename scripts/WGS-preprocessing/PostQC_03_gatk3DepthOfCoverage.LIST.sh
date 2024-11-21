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

postqc_dir=~/data/PCAWG_other/GATK4_pre-processing_bwamem-y/PostQC

mkdir -p ${sampleid}
cd ${sampleid}


ID=${sampleid}
#ID+N/P
INPUT_BAM=${ID}-hg38.bam
POSTQC_DIR=${postqc_dir}
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis/WGS_Data/PostQC
REF=~/data/reference/hg38/Homo_sapiens_assembly38.fasta
INTERVALS=~/data/reference/hg38/Homo_sapiens_assembly38_Autosomal.interval_list
DB_SNP=~/data/reference/SNP/GCF_000001405.38.dbSNP153.GRCh38p12.GATK.vcf.gz
AUTOSOMAL_PAR=~/data/reference/hg38/Homo_sapiens_assembly38_Autosomal_PAR.bed

#-pe OpenMP 4 -l s_vmem=6G

source ~/conda-pack/GATK3.8/bin/activate
unset JAVA_TOOL_OPTIONS
unset PERL5LIB

mkdir -p DepthOfCoverage

#gatk3.8 DepthOfCoverage <Autosomal_PAR>
/usr/bin/time -f "Memory:%M KB time:%E" -o DepthOfCoverage/DepthOfCoverage_A_PAR.txt \
GenomeAnalysisTK -Xmx6g -Xms6g -XX:ParallelGCThreads=4\
 -T DepthOfCoverage\
 -R ${REF}\
 -I ${INPUT_BAM}\
 -nt 4\
 -o DepthOfCoverage/${ID}_Autosomal_PAR\
 -L ${AUTOSOMAL_PAR}\
 -omitBaseOutput\
 -omitIntervals\
 > DepthOfCoverage/DepthOfCoverage_A_PAR.log 2>&1

#gatk3.8 DepthOfCoverage <per Chrromosome>
for CHROM in $(seq 1 22) X Y;
do
	/usr/bin/time -f "Memory:%M KB time:%E" -o DepthOfCoverage/DepthOfCoverage_chr$CHROM.txt \
	GenomeAnalysisTK -Xmx6g -Xms6g -XX:ParallelGCThreads=4\
	 -T DepthOfCoverage\
	 -R ${REF}\
	 -I ${INPUT_BAM}\
	 -nt 4\
     -o DepthOfCoverage/${ID}_chr$CHROM\
     -L chr$CHROM\
     -omitBaseOutput\
     -omitIntervals\
     > DepthOfCoverage/DepthOfCoverage_chr$CHROM.log 2>&1;
done

mkdir -p ${POSTQC_DIR}/DepthOfCoverage

head -n 1 DepthOfCoverage/${ID}_Autosomal_PAR.sample_summary > DepthOfCoverage/header.tsv
cat DepthOfCoverage/${ID}_Autosomal_PAR.sample_summary DepthOfCoverage/${ID}_chr*_summary | grep "${ID}" >> DepthOfCoverage/header.tsv
echo -e chr\\nAutosomal_PAR\\nchr1\\nchr2\\nchr3\\nchr4\\nchr5\\nchr6\\nchr7\\nchr8\\nchr9\\nchr10\\nchr11\\nchr12\\nchr13\\nchr14\\nchr15\\nchr16\\nchr17\\nchr18\\nchr19\\nchr20\\nchr21\\nchr22\\nchrX\\nchrY > DepthOfCoverage/chr.tsv
paste DepthOfCoverage/chr.tsv DepthOfCoverage/header.tsv > ${POSTQC_DIR}/DepthOfCoverage/${ID}_DepthOfCoverage.tsv


