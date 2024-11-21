#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V

ID=$1
# with P
Purple=$2
# PURPLE_SV_GRIPSS/DO14282/PURPLE_SV_plus_medcov/PURPLE/DO14282P.purple.cnv.somatic.tsv
Sclust_iCN=$3
# step4_SV_mod (or step5_SV_mod) /DO14282P_iCN.seg #WGDの確認を！
Sclust_allelic=$4
# step4_SV_mod/DO14282P_allelic_states.txt
Sclust_subclonal=$5
# step4_SV_mod/DO14282P_subclonal_cn.txt
FACETS_cval400=$6
# depth15_cval400/DO14282P.vcf.gz
FACETS_cval50=$7
# depth15_cval50/DO14282P.vcf.gz

mkdir -p ${ID}
cd ${ID}

# -pe OpenMP 1, -l s_vmem=4G, -l m_mem_free=4G


mkdir -p Purple
mkdir -p Sclust
mkdir -p FACETS

## chr1-22 TCN: Purple vs Sclust_iCN vs FACETS
## chr1-22 LCN: Purple vs Sclust_allelic vs FACETS
## chrX TCN: Purple
## chrX LCN: Purple
## chrY TCN: Purple
## chrY LCN: Purple



# Purple CNV
if [ ! -f "Purple/${ID}.Purple_CNV.tsv" ]; then
  less ${Purple} | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$15}' > Purple/${ID}.Purple_CNV.tsv
fi

# Sclust iCN: chr1-22 TCN
if [ ! -f "Sclust/${ID}.Sclust_CNV.iCN.tsv" ]; then
  echo -e "chr\tstart\tend\tiCN" > Sclust/${ID}.Sclust_CNV.iCN.tsv
  cat ${Sclust_iCN} | awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$6}' >> Sclust/${ID}.Sclust_CNV.iCN.tsv
fi


# Sclust _allelic_states.txt + _subclonal_cn.txt: LCN chr1-22,Y   ***なぜかchrXの情報なし！
if [ ! -f "Sclust/${ID}.Sclust_CNV.allelic.tsv" ]; then
  echo -e "chr\tstart\tend\tCopyNr_Raw\tB_CopyNr" > Sclust/${ID}.Sclust_CNV.allelic.tsv
  cat <(tail -n +2 ${Sclust_allelic} | awk -F "\t" '$13==0{print $2"\t"$3"\t"$4"\t"$5"\t"$8}')\
   <(tail -n +2 ${Sclust_subclonal} | awk -F "\t" '{print $2"\t"$3"\t"$4"\t"$5"\t"$7*$8+$10*$11}') | sort -V >> Sclust/${ID}.Sclust_CNV.allelic.tsv
fi

# FACETS_cval400:  Purple_seg > 1Mb: Use cval 400 results
if [ ! -f "FACETS/${ID}.FACETS_CNV.cval400.tsv" ]; then
  purity=`less ${FACETS_cval400} | grep "##purity=" | awk -F "##purity=" '{print $2}'`
  dipLogR=`less ${FACETS_cval400} | grep "##dipLogR=" | awk -F "##dipLogR=" '{print $2}'`
  echo -e "chr\tstart\tend\tTCN\tLCN" > FACETS/${ID}.FACETS_CNV.cval400.tsv
  less ${FACETS_cval400} | grep -v ^# | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
    CNLR_MEDIAN=`echo $line | awk '{print $8}' | awk -F "CNLR_MEDIAN=" '{print $2}' | awk -F ";" '{print $1}'`
    Raw_TCN=`echo "$CNLR_MEDIAN $dipLogR $purity" | awk '{print ((2*2^($1-$2)-2*(1-$3))/$3);}'`
    TCN=`printf "%.2f" ${Raw_TCN}`
    LCN=`echo $line | awk '{print $8}' | awk -F "LCN_EM=" '{print $2}' | awk -F ";" '{print $1}'`
    echo -e "${chr}\t${start}\t${end}\t${TCN}\t${LCN}"
  done >> FACETS/${ID}.FACETS_CNV.cval400.tsv
fi

# FACETS_cval50:  Purple_seg <= 1Mb: Use cval 50 results
if [ ! -f "FACETS/${ID}.FACETS_CNV.cval50.tsv" ]; then
  purity=`less ${FACETS_cval50} | grep "##purity=" | awk -F "##purity=" '{print $2}'`
  dipLogR=`less ${FACETS_cval50} | grep "##dipLogR=" | awk -F "##dipLogR=" '{print $2}'`
  echo -e "chr\tstart\tend\tTCN\tLCN" > FACETS/${ID}.FACETS_CNV.cval50.tsv
  less ${FACETS_cval50} | grep -v ^# | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
    CNLR_MEDIAN=`echo $line | awk '{print $8}' | awk -F "CNLR_MEDIAN=" '{print $2}' | awk -F ";" '{print $1}'`
    Raw_TCN=`echo "$CNLR_MEDIAN $dipLogR $purity" | awk '{print ((2*2^($1-$2)-2*(1-$3))/$3);}'`
    TCN=`printf "%.2f" ${Raw_TCN}`
    LCN=`echo $line | awk '{print $8}' | awk -F "LCN_EM=" '{print $2}' | awk -F ";" '{print $1}'`
    echo -e "${chr}\t${start}\t${end}\t${TCN}\t${LCN}"
  done >> FACETS/${ID}.FACETS_CNV.cval50.tsv
fi




### segmentation bed
if [ ! -f "FACETS/${ID}.FACETS_seg.cval50.sorted.bed" ]; then
  tail -n +2 Purple/${ID}.Purple_CNV.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3}' | sort -k1,1 -k2,2n > Purple/${ID}.Purple_seg.sorted.bed
  tail -n +2 Sclust/${ID}.Sclust_CNV.iCN.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3}' | sort -k1,1 -k2,2n > Sclust/${ID}.Sclust_seg.iCN.sorted.bed
  tail -n +2 Sclust/${ID}.Sclust_CNV.allelic.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3}' | sort -k1,1 -k2,2n > Sclust/${ID}.Sclust_seg.allelic.sorted.bed
  tail -n +2 FACETS/${ID}.FACETS_CNV.cval400.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3}' | sort -k1,1 -k2,2n > FACETS/${ID}.FACETS_seg.cval400.sorted.bed
  tail -n +2 FACETS/${ID}.FACETS_CNV.cval50.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3}' | sort -k1,1 -k2,2n > FACETS/${ID}.FACETS_seg.cval50.sorted.bed
fi


# segmentation intersect
if [ ! -f "FACETS/${ID}.Purple_FACETS_seg.cval50.sorted.bed" ]; then
  bedtools intersect -a Purple/${ID}.Purple_seg.sorted.bed -b Sclust/${ID}.Sclust_seg.iCN.sorted.bed -sorted > Sclust/${ID}.Purple_Sclust_seg.iCN.sorted.bed
  bedtools intersect -a Purple/${ID}.Purple_seg.sorted.bed -b Sclust/${ID}.Sclust_seg.allelic.sorted.bed -sorted > Sclust/${ID}.Purple_Sclust_seg.allelic.sorted.bed
  bedtools intersect -a Purple/${ID}.Purple_seg.sorted.bed -b FACETS/${ID}.FACETS_seg.cval400.sorted.bed -sorted > FACETS/${ID}.Purple_FACETS_seg.cval400.sorted.bed
  bedtools intersect -a Purple/${ID}.Purple_seg.sorted.bed -b FACETS/${ID}.FACETS_seg.cval50.sorted.bed -sorted > FACETS/${ID}.Purple_FACETS_seg.cval50.sorted.bed
fi

# Purple total minor cn bedgraph 5+1
if [ ! -f "Purple/${ID}.Purple_CNV.sorted.bedgraph" ]; then
  tail -n +2 Purple/${ID}.Purple_CNV.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\tPurple"}' | sort -k1,1 -k2,2n > Purple/${ID}.Purple_CNV.sorted.bedgraph
fi

# Sclust total cn bedgraph 4+1
if [ ! -f "Sclust/${ID}.Sclust_CNV.iCN.sorted.bedgraph" ]; then
  tail -n +2 Sclust/${ID}.Sclust_CNV.iCN.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$4"\tSclust"}' | sort -k1,1 -k2,2n > Sclust/${ID}.Sclust_CNV.iCN.sorted.bedgraph
fi

# Sclust total minor cn bedgraph 5+1
if [ ! -f "Sclust/${ID}.Sclust_CNV.allelic.sorted.bedgraph" ]; then
  tail -n +2 Sclust/${ID}.Sclust_CNV.allelic.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\tSclust"}' | sort -k1,1 -k2,2n > Sclust/${ID}.Sclust_CNV.allelic.sorted.bedgraph
fi

# FACETS total minor cn bedgraph 5+1
if [ ! -f "FACETS/${ID}.FACETS_CNV.cval50.sorted.bedgraph" ]; then
  tail -n +2 FACETS/${ID}.FACETS_CNV.cval400.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\tFACETS"}' | sort -k1,1 -k2,2n > FACETS/${ID}.FACETS_CNV.cval400.sorted.bedgraph
  tail -n +2 FACETS/${ID}.FACETS_CNV.cval50.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\tFACETS"}' | sort -k1,1 -k2,2n > FACETS/${ID}.FACETS_CNV.cval50.sorted.bedgraph
fi



# Sclust, FACETS bedgraph splited by 100
mkdir -p tmp

if [ ! -f "Sclust/${ID}.Sclust_CNV.iCN.sorted.100bp.bedgraph" ]; then
  cat Sclust/${ID}.Sclust_CNV.iCN.sorted.bedgraph | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN=`echo $line | awk '{print $4}'`
    seq $start 100 $end > tmp/seq1
    cat <(tail -n +2 tmp/seq1) <(echo $end) > tmp/seq2
    paste tmp/seq1 tmp/seq2 | awk -F "\t" -v chr=${chr} -v TCN=${TCN} '{print chr"\t"$1"\t"$2"\t"TCN"\tSclust"}'
  done > Sclust/${ID}.Sclust_CNV.iCN.sorted.100bp.bedgraph
fi

if [ ! -f "Sclust/${ID}.Sclust_CNV.allelic.sorted.100bp.bedgraph" ]; then
  cat Sclust/${ID}.Sclust_CNV.allelic.sorted.bedgraph | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN=`echo $line | awk '{print $4}'`
    LCN=`echo $line | awk '{print $5}'`
    seq $start 100 $end > tmp/seq1
    cat <(tail -n +2 tmp/seq1) <(echo $end) > tmp/seq2
    paste tmp/seq1 tmp/seq2 | awk -F "\t" -v chr=${chr} -v TCN=${TCN} -v LCN=${LCN} '{print chr"\t"$1"\t"$2"\t"TCN"\t"LCN"\tSclust"}'
  done > Sclust/${ID}.Sclust_CNV.allelic.sorted.100bp.bedgraph
fi

if [ ! -f "FACETS/${ID}.FACETS_CNV.cval50.sorted.100bp.bedgraph" ]; then
  cat FACETS/${ID}.FACETS_CNV.cval50.sorted.bedgraph | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN=`echo $line | awk '{print $4}'`
    LCN=`echo $line | awk '{print $5}'`
    seq $start 100 $end > tmp/seq1
    cat <(tail -n +2 tmp/seq1) <(echo $end) > tmp/seq2
    paste tmp/seq1 tmp/seq2 | awk -F "\t" -v chr=${chr} -v TCN=${TCN} -v LCN=${LCN} '{print chr"\t"$1"\t"$2"\t"TCN"\t"LCN"\tFACETS"}'
  done > FACETS/${ID}.FACETS_CNV.cval50.sorted.100bp.bedgraph
fi

if [ ! -f "FACETS/${ID}.FACETS_CNV.cval50.sorted.100bp.LCN.bedgraph" ]; then
  cat FACETS/${ID}.FACETS_CNV.cval50.sorted.bedgraph | awk -F "\t" '$5!="."' | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN=`echo $line | awk '{print $4}'`
    LCN=`echo $line | awk '{print $5}'`
    seq $start 100 $end > tmp/seq1
    cat <(tail -n +2 tmp/seq1) <(echo $end) > tmp/seq2
    paste tmp/seq1 tmp/seq2 | awk -F "\t" -v chr=${chr} -v TCN=${TCN} -v LCN=${LCN} '{print chr"\t"$1"\t"$2"\t"TCN"\t"LCN"\tFACETS"}'
  done > FACETS/${ID}.FACETS_CNV.cval50.sorted.100bp.LCN.bedgraph
fi

if [ ! -f "FACETS/${ID}.FACETS_CNV.cval400.sorted.100bp.bedgraph" ]; then
  cat FACETS/${ID}.FACETS_CNV.cval400.sorted.bedgraph | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN=`echo $line | awk '{print $4}'`
    LCN=`echo $line | awk '{print $5}'`
    seq $start 100 $end > tmp/seq1
    cat <(tail -n +2 tmp/seq1) <(echo $end) > tmp/seq2
    paste tmp/seq1 tmp/seq2 | awk -F "\t" -v chr=${chr} -v TCN=${TCN} -v LCN=${LCN} '{print chr"\t"$1"\t"$2"\t"TCN"\t"LCN"\tFACETS"}'
  done > FACETS/${ID}.FACETS_CNV.cval400.sorted.100bp.bedgraph
fi

if [ ! -f "FACETS/${ID}.FACETS_CNV.cval400.sorted.100bp.LCN.bedgraph" ]; then
  cat FACETS/${ID}.FACETS_CNV.cval400.sorted.bedgraph | awk -F "\t" '$5!="."' | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN=`echo $line | awk '{print $4}'`
    LCN=`echo $line | awk '{print $5}'`
    seq $start 100 $end > tmp/seq1
    cat <(tail -n +2 tmp/seq1) <(echo $end) > tmp/seq2
    paste tmp/seq1 tmp/seq2 | awk -F "\t" -v chr=${chr} -v TCN=${TCN} -v LCN=${LCN} '{print chr"\t"$1"\t"$2"\t"TCN"\t"LCN"\tFACETS"}'
  done > FACETS/${ID}.FACETS_CNV.cval400.sorted.100bp.LCN.bedgraph
fi



# map Purpleとのintersect segmentの中での平均（結果としてsegment長で重みづけされた平均） ※最頻値(mode)でやると0.01刻みなので一番segmentが長いところになってずれることがある。
if [ ! -f "FACETS/${ID}.FACETS_Minor.Puple_seg.cval400.sorted.bedgraph" ]; then
  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b Sclust/${ID}.Sclust_CNV.iCN.sorted.100bp.bedgraph -c 4 -o mean -sorted > Sclust/${ID}.Sclust_Total.Puple_seg.iCN.sorted.bedgraph

  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b Sclust/${ID}.Sclust_CNV.allelic.sorted.100bp.bedgraph -c 4 -o mean -sorted > Sclust/${ID}.Sclust_Total.Puple_seg.allelic.sorted.bedgraph
  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b Sclust/${ID}.Sclust_CNV.allelic.sorted.100bp.bedgraph -c 5 -o mean -sorted > Sclust/${ID}.Sclust_Minor.Puple_seg.allelic.sorted.bedgraph

  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b FACETS/${ID}.FACETS_CNV.cval50.sorted.100bp.bedgraph -c 4 -o mean -sorted > FACETS/${ID}.FACETS_Total.Puple_seg.cval50.sorted.bedgraph
  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b FACETS/${ID}.FACETS_CNV.cval50.sorted.100bp.LCN.bedgraph -c 5 -o mean -sorted > FACETS/${ID}.FACETS_Minor.Puple_seg.cval50.sorted.bedgraph

  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b FACETS/${ID}.FACETS_CNV.cval400.sorted.100bp.bedgraph -c 4 -o mean -sorted > FACETS/${ID}.FACETS_Total.Puple_seg.cval400.sorted.bedgraph
  bedtools map -a Purple/${ID}.Purple_seg.sorted.bed -b FACETS/${ID}.FACETS_CNV.cval400.sorted.100bp.LCN.bedgraph -c 5 -o mean -sorted > FACETS/${ID}.FACETS_Minor.Puple_seg.cval400.sorted.bedgraph
fi



# bedgraph
if [ ! -f "FACETS/${ID}.FACETS_CNV.Puple_seg.cval400.sorted.bedgraph" ]; then
  cat Sclust/${ID}.Sclust_Total.Puple_seg.iCN.sorted.bedgraph | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\tSclust"}' > Sclust/${ID}.Sclust_CNV.Puple_seg.iCN.sorted.bedgraph

  paste Sclust/${ID}.Sclust_Total.Puple_seg.allelic.sorted.bedgraph Sclust/${ID}.Sclust_Minor.Puple_seg.allelic.sorted.bedgraph | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\tSclust"}' > Sclust/${ID}.Sclust_CNV.Puple_seg.allelic.sorted.bedgraph

  paste FACETS/${ID}.FACETS_Total.Puple_seg.cval50.sorted.bedgraph FACETS/${ID}.FACETS_Minor.Puple_seg.cval50.sorted.bedgraph | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\tFACETS"}' > FACETS/${ID}.FACETS_CNV.Puple_seg.cval50.sorted.bedgraph
  paste FACETS/${ID}.FACETS_Total.Puple_seg.cval400.sorted.bedgraph FACETS/${ID}.FACETS_Minor.Puple_seg.cval400.sorted.bedgraph | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\tFACETS"}' > FACETS/${ID}.FACETS_CNV.Puple_seg.cval400.sorted.bedgraph
fi


# bedgraph merge
if [ ! -f "${ID}.ALL_CNV.Purple_seg.sorted.bedgraph" ]; then
  echo -e "chromosome\tstart\tend\tTCN_Purple\tLCN_Purple\tTCN_Sclust_iCN\tTCN_Sclust_allelic\tLCN_Sclust_allelic\tTCN_FACETS_cval50\tLCN_FACETS_cval50\tTCN_FACETS_cval400\tLCN_FACETS_cval400" > ${ID}.ALL_CNV.Purple_seg.sorted.bedgraph
  paste Purple/${ID}.Purple_CNV.sorted.bedgraph Sclust/${ID}.Sclust_CNV.Puple_seg.iCN.sorted.bedgraph Sclust/${ID}.Sclust_CNV.Puple_seg.allelic.sorted.bedgraph FACETS/${ID}.FACETS_CNV.Puple_seg.cval50.sorted.bedgraph FACETS/${ID}.FACETS_CNV.Puple_seg.cval400.sorted.bedgraph |\
   awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$10"\t"$15"\t"$16"\t"$21"\t"$22"\t"$27"\t"$28}' | sort -k1,1V -k2,2n >> ${ID}.ALL_CNV.Purple_seg.sorted.bedgraph
fi

# select Sclust, FACETS
if [ ! -f "${ID}.ALL_FINAL_CNV.Purple_seg.bedgraph" ]; then
  echo -e "chromosome\tstart\tend\tTCN_Purple\tTCN_Sclust\tTCN_FACETS\tLCN_Purple\tLCN_Sclust\tLCN_FACETS" > ${ID}.ALL_FINAL_CNV.Purple_seg.bedgraph
  tail -n +2 ${ID}.ALL_CNV.Purple_seg.sorted.bedgraph | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    seg_length=`echo ${start} ${end} | awk '{print $2-$1}'`
    TCN_Purple=`echo $line | awk '{print $4}'`
    LCN_Purple=`echo $line | awk '{print $5}'`
    TCN_Sclust=`echo $line | awk '{print $6}'`
    LCN_Sclust=`echo $line | awk '{print $8}'`
    if [[ "${seg_length}" -gt 1000000 ]]; then
      TCN_FACETS=`echo $line | awk '{print $11}'`
      LCN_FACETS=`echo $line | awk '{print $12}'`
    else
      TCN_FACETS=`echo $line | awk '{print $9}'`
      LCN_FACETS=`echo $line | awk '{print $10}'`
    fi
    echo -e "${chr}\t${start}\t${end}\t${TCN_Purple}\t${TCN_Sclust}\t${TCN_FACETS}\t${LCN_Purple}\t${LCN_Sclust}\t${LCN_FACETS}"
  done | sort -V >>  ${ID}.ALL_FINAL_CNV.Purple_seg.bedgraph
fi

if [ ! -f "${ID}.ALL_FINAL_CNV.Purple_seg.tsv" ]; then
  echo -e "chromosome\tstart\tend\tTCN_Purple\tTCN_Sclust\tTCN_FACETS\tMin_TCN\tMedian_TCN\tMean_TCN\tMax_TCN\tLCN_Purple\tLCN_Sclust\tLCN_FACETS\tMin_LCN\tMedian_LCN\tMean_LCN\tMax_LCN" > ${ID}.ALL_FINAL_CNV.Purple_seg.tsv
  tail -n +2 ${ID}.ALL_FINAL_CNV.Purple_seg.bedgraph | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2+1}'`
    end=`echo $line | awk '{print $3}'`
    T_P=`echo $line | awk '{if($4=="."){print "NA"}else{print $4}}'`
    T_S=`echo $line | awk '{if($5=="."){print "NA"}else{print $5}}'`
    T_F=`echo $line | awk '{if($6=="."){print "NA"}else{print $6}}'`
    Rscript_TCN=`Rscript -e "summary(c(${T_P},${T_S},${T_F}))" | tail -n +2`
    Min_TCN=`echo $Rscript_TCN | awk '{print $1}'`
    Median_TCN=`echo $Rscript_TCN | awk '{print $3}'`
    Mean_TCN=`echo $Rscript_TCN | awk '{print $4}'`
    Max_TCN=`echo $Rscript_TCN | awk '{print $6}'`
    L_P=`echo $line | awk '{if($7=="."){print "NA"}else{print $7}}'`
    L_S=`echo $line | awk '{if($8=="."){print "NA"}else{print $8}}'`
    L_F=`echo $line | awk '{if($9=="."){print "NA"}else{print $9}}'`
    Rscript_LCN=`Rscript -e "summary(c(${L_P},${L_S},${L_F}))" | tail -n +2`
    Min_LCN=`echo $Rscript_LCN | awk '{print $1}'`
    Median_LCN=`echo $Rscript_LCN | awk '{print $3}'`
    Mean_LCN=`echo $Rscript_LCN | awk '{print $4}'`
    Max_LCN=`echo $Rscript_LCN | awk '{print $6}'`
    echo -e "${chr}\t${start}\t${end}\t${T_P}\t${T_S}\t${T_F}\t${Min_TCN}\t${Median_TCN}\t${Mean_TCN}\t${Max_TCN}\t${L_P}\t${L_S}\t${L_F}\t${Min_LCN}\t${Median_LCN}\t${Mean_LCN}\t${Max_LCN}"
  done >> ${ID}.ALL_FINAL_CNV.Purple_seg.tsv
fi


if [ ! -f "${ID}.ALL_FINAL_CNV.merge.tsv" ]; then
  echo -e "chromosome\tstart\tend\tTCN_Purple\tMedian_TCN\tMean_TCN\tLCN_Purple\tMedian_LCN\tMean_LCN" > ${ID}.ALL_FINAL_CNV.merge.tsv
  tail -n +2 ${ID}.ALL_FINAL_CNV.Purple_seg.tsv | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN_Purple=`echo $line | awk '{print $4}'`
    LCN_Purple=`echo $line | awk '{print $11}'`
    if [[ `echo $line | awk '{print $4,$5,$6}' | grep "NA"` ]] ; then
      Median_TCN=`echo $TCN_Purple`
      Mean_TCN=`echo $TCN_Purple`
    else
      Median_TCN=`echo $line | awk '{print $8}'`
      Mean_TCN=`echo $line | awk '{print $9}'`
    fi
    if [[ `echo $line | awk '{print $11,$12,$13}' | grep "NA"` ]] ; then
      Median_LCN=`echo $LCN_Purple`
      Mean_LCN=`echo $LCN_Purple`
    else
      Median_LCN=`echo $line | awk '{print $15}'`
      Mean_LCN=`echo $line | awk '{print $16}'`
    fi
    echo -e "$chr\t$start\t$end\t$TCN_Purple\t$Median_TCN\t$Mean_TCN\t$LCN_Purple\t$Median_LCN\t$Mean_LCN"
  done >> ${ID}.ALL_FINAL_CNV.merge.tsv
fi




### major_cnを用いた修正案
if [ ! -f "${ID}.ALL_FINAL_CNV.Purple_seg.mod.tsv" ]; then
  echo -e "chromosome\tstart\tend\tTotal_Purple\tTotal_Sclust\tTotal_FACETS\tMin_Total\tMedian_Total\tMean_Total\tMax_Total\tminor_Purple\tminor_Sclust\tminor_FACETS\tMin_minor\tMedian_minor\tMean_minor\tMax_minor" > ${ID}.ALL_FINAL_CNV.Purple_seg.mod.tsv
  tail -n +2 ${ID}.ALL_FINAL_CNV.Purple_seg.tsv | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    T_P=`echo $line | awk '{if($4=="."){print "NA"}else{print $4}}'`
    T_S=`echo $line | awk '{if($5=="."){print "NA"}else{print $5}}'`
    T_F=`echo $line | awk '{if($6=="."){print "NA"}else{print $6}}'`
    L_P=`echo $line | awk '{if($7=="."){print "NA"}else{print $11}}'`
    L_S=`echo $line | awk '{if($8=="."){print "NA"}else{print $12}}'`
    L_F=`echo $line | awk '{if($9=="."){print "NA"}else{print $13}}'`
    T_P_bc=`echo $T_P | sed -e 's/e/*10^/g'`
    T_S_bc=`echo $T_S | sed -e 's/e/*10^/g'`
    T_F_bc=`echo $T_F | sed -e 's/e/*10^/g'`
    L_P_bc=`echo $L_P | sed -e 's/e/*10^/g'`
    L_S_bc=`echo $L_S | sed -e 's/e/*10^/g'`
    L_F_bc=`echo $L_F | sed -e 's/e/*10^/g'`
    if [[ ! `echo ${T_P} ${L_P} | grep "NA"` ]] ; then
      # TCN <= LCNならばTotal=TCN minor=0
      result=`echo "${T_P_bc} <= ${L_P_bc}" | bc`
      if [[ ${result} -eq 1 ]]; then
        Total_Purple=`echo ${T_P}`
        minor_Purple="0"
      else
        M_P=`echo ${T_P} ${L_P} | awk '{print $1-$2}'`
        M_P_bc=`echo $M_P | sed -e 's/e/*10^/g'`
        # TCN-LCN <= LCNならばminor=TCN-LCN
        result=`echo "${M_P_bc} <= ${L_P_bc}" | bc`
        if [[ ${result} -eq 1 ]]; then
          Total_Purple=`echo ${T_P}`
          minor_Purple=`echo ${M_P}`
        else
          Total_Purple=`echo ${T_P}`
          minor_Purple=`echo ${L_P}`
        fi
      fi
    else
      Total_Purple=`echo ${T_P}`
      minor_Purple=`echo ${L_P}`
    fi
    if [[ ! `echo ${T_S} ${L_S} | grep "NA"` ]] ; then
      # TCN <= LCNならばTotal=TCN minor=0
      result=`echo "${T_S_bc} <= ${L_S_bc}" | bc`
      if [[ ${result} -eq 1 ]]; then
        TotaL_Sclust=`echo ${T_S}`
        minor_Sclust="0"
      else
        M_S=`echo ${T_S} ${L_S} | awk '{print $1-$2}'`
        M_S_bc=`echo $M_S | sed -e 's/e/*10^/g'`
        # TCN-LCN <= LCNならばminor=TCN-LCN
        result=`echo "${M_S_bc} <= ${L_S_bc}" | bc`
        if [[ ${result} -eq 1 ]]; then
          TotaL_Sclust=`echo ${T_S}`
          minor_Sclust=`echo ${M_S}`
        else
          TotaL_Sclust=`echo ${T_S}`
          minor_Sclust=`echo ${L_S}`
        fi
      fi
    else
      TotaL_Sclust=`echo ${T_S}`
      minor_Sclust=`echo ${L_S}`
    fi
    if [[ ! `echo ${T_F} ${L_F} | grep "NA"` ]] ; then
      # TCN <= LCNならばTotal=TCN minor=0
      result=`echo "${T_F_bc} <= ${L_F_bc}" | bc`
      if [[ ${result} -eq 1 ]]; then
        TotaL_FACETS=`echo ${T_F}`
        minor_FACETS="0"
      else
        M_F=`echo ${T_F} ${L_F} | awk '{print $1-$2}'`
        M_F_bc=`echo $M_F | sed -e 's/e/*10^/g'`
        # TCN-LCN <= LCNならばminor=TCN-LCN
        result=`echo "${M_F_bc} <= ${L_F_bc}" | bc`
        if [[ ${result} -eq 1 ]]; then
          TotaL_FACETS=`echo ${T_F}`
          minor_FACETS=`echo ${M_F}`
        else
          TotaL_FACETS=`echo ${T_F}`
          minor_FACETS=`echo ${L_F}`
        fi
      fi
    else
      TotaL_FACETS=`echo ${T_F}`
      minor_FACETS=`echo ${L_F}`
    fi
    Rscript_Total=`Rscript -e "summary(c(${Total_Purple},${TotaL_Sclust},${TotaL_FACETS}))" | tail -n +2`
    Min_Total=`echo $Rscript_Total | awk '{print $1}'`
    Median_Total=`echo $Rscript_Total | awk '{print $3}'`
    Mean_Total=`echo $Rscript_Total | awk '{print $4}'`
    Max_Total=`echo $Rscript_Total | awk '{print $6}'`
    Rscript_minor=`Rscript -e "summary(c(${minor_Purple},${minor_Sclust},${minor_FACETS}))" | tail -n +2`
    Min_minor=`echo $Rscript_minor | awk '{print $1}'`
    Median_minor=`echo $Rscript_minor | awk '{print $3}'`
    Mean_minor=`echo $Rscript_minor | awk '{print $4}'`
    Max_minor=`echo $Rscript_minor | awk '{print $6}'`
    echo -e "${chr}\t${start}\t${end}\t${Total_Purple}\t${TotaL_Sclust}\t${TotaL_FACETS}\t${Min_Total}\t${Median_Total}\t${Mean_Total}\t${Max_Total}\t${minor_Purple}\t${minor_Sclust}\t${minor_FACETS}\t${Min_minor}\t${Median_minor}\t${Mean_minor}\t${Max_minor}"
  done >> ${ID}.ALL_FINAL_CNV.Purple_seg.mod.tsv
fi

if [ ! -f "${ID}.ALL_FINAL_CNV.merge.mod.tsv" ]; then
  echo -e "chromosome\tstart\tend\tTotal_Purple\tMedian_Total\tMean_Total\tminor_Purple\tMedian_minor\tMean_minor" > ${ID}.ALL_FINAL_CNV.merge.mod.tsv
  tail -n +2 ${ID}.ALL_FINAL_CNV.Purple_seg.mod.tsv | while read line
  do
    chr=`echo $line | awk '{print $1}'`
    start=`echo $line | awk '{print $2}'`
    end=`echo $line | awk '{print $3}'`
    TCN_Purple=`echo $line | awk '{if($4<0){print "0"}else{print $4}}'`
    LCN_Purple=`echo $line | awk '{if($11<0){print "0"}else{print $11}}'`
    TCN_Purple_bc=`echo $TCN_Purple | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
    LCN_Purple_bc=`echo $LCN_Purple | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
    if [[ `echo $line | awk '{print $4,$5,$6}' | grep "NA"` ]] ; then
      Median_TCN=`echo $TCN_Purple`
      Mean_TCN=`echo $TCN_Purple`
      Median_TCN_bc=`echo $Median_TCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
      Mean_TCN_bc=`echo $Mean_TCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
    else
      Median_TCN=`echo $line | awk '{if($8<0){print "0"}else{print $8}}'`
      Mean_TCN=`echo $line | awk '{if($9<0){print "0"}else{print $9}}'`
      Median_TCN_bc=`echo $Median_TCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
      Mean_TCN_bc=`echo $Mean_TCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
    fi
    if [[ `echo $line | awk '{print $11,$12,$13}' | grep "NA"` ]] ; then
      Median_LCN=`echo $LCN_Purple`
      Mean_LCN=`echo $LCN_Purple`
      Median_LCN_bc=`echo $Median_LCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
      Mean_LCN_bc=`echo $Mean_LCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
    else
      Median_LCN=`echo $line | awk '{if($15<0){print "0"}else{print $15}}'`
      Mean_LCN=`echo $line | awk '{if($16<0){print "0"}else{print $16}}'`
      Median_LCN_bc=`echo $Median_LCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
      Mean_LCN_bc=`echo $Mean_LCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
    fi
    # TCN <= LCNならばTotal=TCN minor=0
    result=`echo "${Median_TCN_bc} <= ${Median_LCN_bc}" | bc`
    if [[ ${result} -eq 1 ]]; then
      Median_Total=`echo ${Median_TCN}`
      Median_minor="0"
    else
      Median_MCN=`echo ${Median_TCN} ${Median_LCN} | awk '{print $1-$2}'`
      Median_MCN_bc=`echo $Median_MCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
      # TCN-LCN <= LCNならばminor=TCN-LCN
      result=`echo "${Median_MCN_bc} <= ${Median_LCN_bc}" | bc`
      if [[ ${result} -eq 1 ]]; then
        Median_Total=`echo ${Median_TCN}`
        Median_minor=`echo ${Median_MCN}`
      else
        Median_Total=`echo ${Median_TCN}`
        Median_minor=`echo ${Median_LCN}`
      fi
    fi
    # TCN <= LCNならばTotal=TCN minor=0
    result=`echo "${Mean_TCN_bc} <= ${Mean_LCN_bc}" | bc`
    if [[ ${result} -eq 1 ]]; then
      Mean_Total=`echo ${Mean_TCN}`
      Mean_minor="0"
    else
      Mean_MCN=`echo ${Mean_TCN} ${Mean_LCN} | awk '{print $1-$2}'`
      Mean_MCN_bc=`echo $Mean_MCN | sed -e 's/e/*10^/g' | sed -e 's/+//g'`
      # TCN-LCN <= LCNならばminor=TCN-LCN
      result=`echo "${Mean_MCN_bc} <= ${Mean_LCN_bc}" | bc`
      if [[ ${result} -eq 1 ]]; then
        Mean_Total=`echo ${Mean_TCN}`
        Mean_minor=`echo ${Mean_MCN}`
      else
        Mean_Total=`echo ${Mean_TCN}`
        Mean_minor=`echo ${Mean_LCN}`
      fi
    fi
    FINAL_Median_Total=`Rscript -e "$Median_Total" | awk '{print $2}'`
    FINAL_Median_minor=`Rscript -e "$Median_minor" | awk '{print $2}'`
    FINAL_Mean_Total=`Rscript -e "$Mean_Total" | awk '{print $2}'`
    FINAL_Mean_minor=`Rscript -e "$Mean_minor" | awk '{print $2}'`
    echo -e "$chr\t$start\t$end\t$TCN_Purple\t$FINAL_Median_Total\t$FINAL_Mean_Total\t$LCN_Purple\t$FINAL_Median_minor\t$FINAL_Mean_minor"
  done >> ${ID}.ALL_FINAL_CNV.merge.mod.tsv
fi


### PURPLE only
if [ ! -f "${ID}.ALL_FINAL_CNV.merge.mod.PURPLE_ONLY.tsv" ]; then
  cat ${ID}.ALL_FINAL_CNV.merge.mod.tsv \
  | awk -F "\t" '{print $1"\t"$2"\t"$3"\ttotal_cn->\t"$4"\t"$5"\tminor_cn->\t"$7"\t"$8}' > ${ID}.ALL_FINAL_CNV.merge.mod.PURPLE_ONLY.tsv
fi  



### CNA binning
source ~/conda-pack/gatk4.1.9.0/bin/activate
unset PERL5LIB
unset JAVA_TOOL_OPTIONS

if [ ! -f "100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz" ]; then
  mkdir -p 100bp_bin
  tail -n +2 ${ID}.ALL_FINAL_CNV.merge.mod.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$5"\t"$8}' > 100bp_bin/${ID}.ALL_FINAL_CNV.merge.mod.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100bp.window.bed.gz -b 100bp_bin/${ID}.ALL_FINAL_CNV.merge.mod.bedgraph -c 4 -o mean | bgzip > 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz
  tabix -p bed 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100bp.window.bed.gz -b 100bp_bin/${ID}.ALL_FINAL_CNV.merge.mod.bedgraph -c 5 -o mean | bgzip > 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz
  tabix -p bed 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz
fi

if [ ! -f "1Mb_bin/${ID}.ALL_FINAL_CNV.1Mb_bin.MinorMedian.bedgraph" ]; then
  mkdir -p 1Mb_bin
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1Mb.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > 1Mb_bin/${ID}.ALL_FINAL_CNV.1Mb_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1Mb.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > 1Mb_bin/${ID}.ALL_FINAL_CNV.1Mb_bin.MinorMedian.bedgraph
fi

if [ ! -f "100kb_bin/${ID}.ALL_FINAL_CNV.100kb_bin.MinorMedian.bedgraph" ]; then
  mkdir -p 100kb_bin
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100kb.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > 100kb_bin/${ID}.ALL_FINAL_CNV.100kb_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100kb.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > 100kb_bin/${ID}.ALL_FINAL_CNV.100kb_bin.MinorMedian.bedgraph
fi

if [ ! -f "1kb_bin/${ID}.ALL_FINAL_CNV.1kb_bin.MinorMedian.bedgraph" ]; then
  mkdir -p 1kb_bin
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1kb.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > 1kb_bin/${ID}.ALL_FINAL_CNV.1kb_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1kb.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > 1kb_bin/${ID}.ALL_FINAL_CNV.1kb_bin.MinorMedian.bedgraph
fi

if [ ! -f "pq_bin/${ID}.ALL_FINAL_CNV.pq_bin.MinorMedian.bedgraph" ]; then
  mkdir -p pq_bin
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.pq.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > pq_bin/${ID}.ALL_FINAL_CNV.pq_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.pq.window.bed -b 100bp_bin/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > pq_bin/${ID}.ALL_FINAL_CNV.pq_bin.MinorMedian.bedgraph
fi




### CNA binning PURPLE_ONLY
#source /work23/home/nsasa/conda-pack/gatk4.1.9.0/bin/activate

if [ ! -f "100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz" ]; then
  mkdir -p 100bp_bin_PURPLE
  tail -n +2 ${ID}.ALL_FINAL_CNV.merge.mod.PURPLE_ONLY.tsv | awk -F "\t" '{print $1"\t"$2-1"\t"$3"\t"$5"\t"$8}' > 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.merge.mod.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100bp.window.bed.gz -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.merge.mod.bedgraph -c 4 -o mean | bgzip > 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz
  tabix -p bed 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100bp.window.bed.gz -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.merge.mod.bedgraph -c 5 -o mean | bgzip > 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz
  tabix -p bed 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz
fi

if [ ! -f "1Mb_bin_PURPLE/${ID}.ALL_FINAL_CNV.1Mb_bin.MinorMedian.bedgraph" ]; then
  mkdir -p 1Mb_bin_PURPLE
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1Mb.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > 1Mb_bin_PURPLE/${ID}.ALL_FINAL_CNV.1Mb_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1Mb.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > 1Mb_bin_PURPLE/${ID}.ALL_FINAL_CNV.1Mb_bin.MinorMedian.bedgraph
fi

if [ ! -f "100kb_bin_PURPLE/${ID}.ALL_FINAL_CNV.100kb_bin.MinorMedian.bedgraph" ]; then
  mkdir -p 100kb_bin_PURPLE
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100kb.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > 100kb_bin_PURPLE/${ID}.ALL_FINAL_CNV.100kb_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.100kb.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > 100kb_bin_PURPLE/${ID}.ALL_FINAL_CNV.100kb_bin.MinorMedian.bedgraph
fi

if [ ! -f "1kb_bin_PURPLE/${ID}.ALL_FINAL_CNV.1kb_bin.MinorMedian.bedgraph" ]; then
  mkdir -p 1kb_bin_PURPLE
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1kb.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > 1kb_bin_PURPLE/${ID}.ALL_FINAL_CNV.1kb_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.1kb.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > 1kb_bin_PURPLE/${ID}.ALL_FINAL_CNV.1kb_bin.MinorMedian.bedgraph
fi

if [ ! -f "pq_bin_PURPLE/${ID}.ALL_FINAL_CNV.pq_bin.MinorMedian.bedgraph" ]; then
  mkdir -p pq_bin_PURPLE
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.pq.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.TotalMedian.bedgraph.gz -c 4 -o mean > pq_bin_PURPLE/${ID}.ALL_FINAL_CNV.pq_bin.TotalMedian.bedgraph
  bedtools map -a /home/n_sasa/data/reference/hg38/bedtools/hg38.1-22.X.Y.pq.window.bed -b 100bp_bin_PURPLE/${ID}.ALL_FINAL_CNV.100bp_bin.MinorMedian.bedgraph.gz -c 4 -o mean > pq_bin_PURPLE/${ID}.ALL_FINAL_CNV.pq_bin.MinorMedian.bedgraph
fi
