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
#~/data/PCAWG_other/GATK4_pre-processing_bwamem-y/PostQC

mkdir -p ${sampleid}
cd ${sampleid}


ID=${sampleid}
#ID+N/P
INPUT_BAM=${ID}-hg38.bam
REF=Homo_sapiens_assembly38.fasta
REF_DIR=~/data/reference/hg38
CONTAMINATION_SITES_UD=Homo_sapiens_assembly38.contam.UD
CONTAMINATION_SITES_MU=Homo_sapiens_assembly38.contam.mu
CONTAMINATION_SITES_BED=Homo_sapiens_assembly38.contam.bed
contamination_underestimation_factor=0.75
POSTQC_DIR=${postqc_dir}
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis/WGS_Data/PostQC

#-pe OpenMP 4 -l s_vmem=6G

. ~/miniconda3/etc/profile.d/conda.sh
#conda activate udocker_1.3.16
conda activate verifybamid2
unset PERL5LIB
#unset UDOCKER_DIR
#export UDOCKER_DIR=~/udocker/.udocker$UDOCKER_DIR

mkdir -p VerifyBamID2

### VerifyBamID2 Estimate level of cross-sample contamination
# creates a ${ID}.selfSM file, a TSV file with 2 rows, 19 colums.
# First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
#/usr/bin/time -f "Memory:%M KB time:%E" -o VerifyBamID2/VerifyBamID2.txt \
#udocker run -v ${REF_DIR}:/mnt/mnt0 -v `pwd`:/mnt/mnt1 --workdir=/mnt/mnt1 verifybamid2 \
#VerifyBamID\
# --Verbose\
# --NumPC 4\
# --Output VerifyBamID2/${ID}\
# --BamFile /mnt/mnt1/${INPUT_BAM}\
# --Reference /mnt/mnt0/${REF}\
# --UDPath /mnt/mnt0/contam/${CONTAMINATION_SITES_UD}\
# --MeanPath /mnt/mnt0/contam/${CONTAMINATION_SITES_MU}\
# --BedPath /mnt/mnt0/contam/${CONTAMINATION_SITES_BED}\
# > VerifyBamID2/VerifyBamID2.log 2>&1

/usr/bin/time -f "Memory:%M KB time:%E" -o VerifyBamID2/VerifyBamID2.txt \
verifybamid2\
 --SVDPrefix ~/data/reference/hg38/contam/Homo_sapiens_assembly38.contam\
 --Verbose\
 --NumPC 4\
 --Output VerifyBamID2/${ID}\
 --BamFile ${INPUT_BAM}\
 --Reference ~/data/reference/hg38/${REF}\
 > VerifyBamID2/VerifyBamID2.log 2>&1

cd VerifyBamID2



# used to read from the selfSM file and calculate contamination, which gets printed out
#source ~/conda-pack/RNA-seq/bin/activate
#unset PERL5LIB

python3 <<CODE
import csv
import sys
with open('${ID}.selfSM') as selfSM:
    reader = csv.DictReader(selfSM, delimiter='\t')
    i = 0
    for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
    	    # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
            # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
            # vcf and bam.
            sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
            sys.exit(1)
        print(float(row["FREEMIX"])/${contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
      	    sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
        sys.exit(2)
CODE


mkdir -p $POSTQC_DIR/VerifyBamID2
cp ${ID}.selfSM $POSTQC_DIR/VerifyBamID2
