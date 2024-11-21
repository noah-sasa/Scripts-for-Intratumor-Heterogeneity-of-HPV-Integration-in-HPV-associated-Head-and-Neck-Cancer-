#!/bin/bash

#$ -S /bin/bash
#$ -cwd
#$ -V

ID=$1
# with P or E/L
M2_VCF=$2
S2_SNV_VCF=$3
S2_INDEL_VCF=$4
GRIDSS_VCF=$5
Manta_VCF=$6
Delly_VCF=$7
REF_FASTA=/home/n_sasa/data/reference/hg38/Homo_sapiens_assembly38.fasta

#-pe OpenMP 2 -l s_vmem=6G

ID_N=$(echo ${ID/%?/}N)

mkdir -p ${ID}
cd ${ID}


source ~/conda-pack/gatk4.1.9.0/bin/activate
unset PERL5LIB
unset JAVA_TOOL_OPTIONS


if [ ! -f "M2/M2.PASS.Norm.noMNP.indel.vcf" ]; then
  ###M2
  mkdir -p M2
  cd M2
  mkdir -p log
  
  ### M2 PASS
  less ${M2_VCF} | grep "^#" > M2.PASS.vcf
  bcftools view -H ${M2_VCF} | awk '$7=="PASS"{print}' >> M2.PASS.vcf
  bgzip M2.PASS.vcf
  tabix -p vcf M2.PASS.vcf.gz
  
  ### M2 Normalization: Split multiallelics into biallelics, left align and trim alleles
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/LeftAlignAndTrimVariants_M2.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   LeftAlignAndTrimVariants\
   -R ${REF_FASTA}\
   -V M2.PASS.vcf.gz\
   -O M2.PASS.Norm.vcf\
   --split-multi-allelics\
   --keep-original-ac\
   > log/LeftAlignAndTrimVariants_M2.log 2>&1
  
  ### extract MNP in M2
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/SelectVariants_extract-MNP.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   SelectVariants\
   -R ${REF_FASTA}\
   -V M2.PASS.Norm.vcf\
   -O M.PASS.Norm.extract-MNP.vcf\
   --select-type-to-include MNP\
   > log/SelectVariants_extract-MNP.log 2>&1
  
  
  ### MNPtoSNP GATK3  INPUT:vcf.gzでWARN
  source ~/conda-pack/GATK3.8/bin/activate
  unset PERL5LIB
  unset JAVA_TOOL_OPTIONS

  /usr/bin/time -f "Memory:%M KB time:%E" -o log/VariantsToAllelicPrimitives.txt \
  GenomeAnalysisTK -Xmx6g -Xms6g -XX:ParallelGCThreads=4\
   -T VariantsToAllelicPrimitives\
   -R ${REF_FASTA}\
   -V M2.PASS.Norm.vcf\
   -o M2.PASS.Norm.noMNP.vcf.gz\
   > log/VariantsToAllelicPrimitives.log 2>&1
  source ~/conda-pack/gatk4.1.9.0/bin/activate
  unset PERL5LIB
  unset JAVA_TOOL_OPTIONS
  
  ### SNV or INDEL
  less M2.PASS.Norm.noMNP.vcf.gz | grep "^#" > M2.PASS.Norm.noMNP.snv.vcf
  less M2.PASS.Norm.noMNP.vcf.gz | grep "^#" > M2.PASS.Norm.noMNP.indel.vcf
  bcftools view -H M2.PASS.Norm.noMNP.vcf.gz | awk 'length($4)==1 && length($5)==1{print}' >> M2.PASS.Norm.noMNP.snv.vcf
  bcftools view -H M2.PASS.Norm.noMNP.vcf.gz | awk 'length($4)>1 || length($5)>1{print}' >> M2.PASS.Norm.noMNP.indel.vcf
  # CombineVcfs vcf.gzでWARN

  cd ..
  echo "M2 Done"
else
  echo "M2 Already Done"
fi


if [ ! -f "S2/S2.INDEL.PASS.Norm.name-changed.vcf" ]; then
  ### S2
  mkdir -p S2
  cd S2
  mkdir -p log
  
  ### S2 PASS
  less ${S2_SNV_VCF} | grep "^#" > S2.SNV.PASS.vcf
  bcftools view -H ${S2_SNV_VCF} | awk '$7=="PASS"{print}' >> S2.SNV.PASS.vcf
  bgzip S2.SNV.PASS.vcf
  tabix -p vcf S2.SNV.PASS.vcf.gz
  
  less ${S2_INDEL_VCF} | grep "^#" > S2.INDEL.PASS.vcf
  bcftools view -H ${S2_INDEL_VCF} | awk '$7=="PASS"{print}' >> S2.INDEL.PASS.vcf
  bgzip S2.INDEL.PASS.vcf
  tabix -p vcf S2.INDEL.PASS.vcf.gz
  
  ### S2 Normalization
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/LeftAlignAndTrimVariants_SNV.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   LeftAlignAndTrimVariants\
   -R ${REF_FASTA}\
   -V S2.SNV.PASS.vcf.gz\
   -O S2.SNV.PASS.Norm.vcf.gz\
   --split-multi-allelics\
   --keep-original-ac\
   > log/LeftAlignAndTrimVariants_SNV.log 2>&1
  
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/LeftAlignAndTrimVariants_INDEL.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   LeftAlignAndTrimVariants\
   -R ${REF_FASTA}\
   -V S2.INDEL.PASS.vcf.gz\
   -O S2.INDEL.PASS.Norm.vcf.gz\
   --split-multi-allelics\
   --keep-original-ac\
   > log/LeftAlignAndTrimVariants_INDEL.log 2>&1
  
  ### change ID name
  less S2.SNV.PASS.Norm.vcf.gz | sed -e "s/NORMAL/${ID_N}/g" | sed -e "s/TUMOR/${ID}/g" > S2.SNV.PASS.Norm.name-changed.vcf
  less S2.INDEL.PASS.Norm.vcf.gz | sed -e "s/NORMAL/${ID_N}/g" | sed -e "s/TUMOR/${ID}/g" > S2.INDEL.PASS.Norm.name-changed.vcf
  # CombineVcfs vcf.gzでWARN

  cd ..
  echo "S2 Done"
else
  echo "S2 Already Done"
fi


if [ ! -f "SNV/M2.S2.SNV.merged.withMNP.vcf" ]; then
  ### CombineVariants GATK3
  mkdir -p SNV
  cd SNV
  mkdir -p log
  mkdir -p CVTmp
  source ~/conda-pack/GATK3.8/bin/activate
  unset PERL5LIB
  unset JAVA_TOOL_OPTIONS
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/CombineVariants.txt \
  GenomeAnalysisTK -Xmx6g -Xms6g -XX:ParallelGCThreads=4 -Djava.io.tmpdir=CVTmp\
   -T CombineVariants\
   -R ${REF_FASTA}\
   -genotypeMergeOptions PRIORITIZE\
   --rod_priority_list mutect,strelka\
   --variant:mutect ../M2/M2.PASS.Norm.noMNP.snv.vcf\
   --variant:strelka ../S2/S2.SNV.PASS.Norm.name-changed.vcf\
   -o M2.S2.SNV.combined.vcf.gz\
   -nt 1\
   > log/CombineVariants.log 2>&1
  
  
  ### extract Mutect2 MNP to SNP -> TARGET
  /usr/bin/time -f "Memory:%M KB time:%E" -o log/VariantsToAllelicPrimitives.txt \
  GenomeAnalysisTK -Xmx6g -Xms6g -XX:ParallelGCThreads=4\
   -T VariantsToAllelicPrimitives\
   -R ${REF_FASTA}\
   -V ../M2/M.PASS.Norm.extract-MNP.vcf\
   -o M2.PASS.Norm.extract-MNP.toSNP.vcf.gz\
   > log/VariantsToAllelicPrimitives.log 2>&1
  
  source ~/conda-pack/gatk4.1.9.0/bin/activate
  unset PERL5LIB
  unset JAVA_TOOL_OPTIONS
  
  # extract multiallelic SNV
  bcftools view -H M2.S2.SNV.combined.vcf.gz | awk '$5~/,/{print}' > M2.S2.SNV.extract-multiallelic
  ## if include INDEL
  #bcftools view -H M2.S2.combined.vcf | awk '{print $1,$2}' |uniq -d > M2.S2.combined.extract-samePOS
  
  ### set=Intersection
  /usr/bin/time -f "Memory:%M KB time:%E" -o ExtractSelectVariants.txt \
  gatk --java-options "-Xmx6g -Xms6g -XX:+UseSerialGC"\
   SelectVariants\
   -R ${REF_FASTA}\
   -V M2.S2.SNV.combined.vcf.gz\
   -select 'set=="Intersection"'\
   -O M2.S2.SNV.combined.intersection.vcf.gz\
   > ExtractSelectVariants.log 2>&1
  
  ### bcftools TARGET
  # Reagionsで絞る場合にはSTARTがRegionsに含まれていなくても、変異自体がoverlapしていたら抽出されてしまう。
  bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' M2.S2.SNV.combined.intersection.vcf.gz | bgzip -c > M2.S2.SNV.tsv.gz && tabix -s1 -b2 -e2 M2.S2.SNV.tsv.gz
  
  less ../M2/M2.PASS.Norm.noMNP.snv.vcf | grep "^#" > M2.S2.SNV.merged.noMNP.vcf
  bcftools view -H ../M2/M2.PASS.Norm.noMNP.snv.vcf -T M2.S2.SNV.tsv.gz >> M2.S2.SNV.merged.noMNP.vcf
  
  
  ### MNPtoSNP -> TARGET
  bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' M2.PASS.Norm.extract-MNP.toSNP.vcf.gz | bgzip -c > M2.MNV.tsv.gz && tabix -s1 -b2 -e2 M2.MNV.tsv.gz
  bcftools view -H M2.S2.SNV.merged.noMNP.vcf -T M2.MNV.tsv.gz > SNPtoMNP
  ### SNPtoMNPのPOSをBEDに直すことでbedtools mergeを使ってSNPtoMNPの連続するPOSをmergeできる。
  cat SNPtoMNP | awk '{print $1"\t"$2-1"\t"$2}' | bedtools merge > SNPtoMNP.bed
  
  ### BEDを1行ずつ読み込み、その範囲をTARGETとしてnoMNP.vcfからMNPに直すべきSNPを抜き出しMNPに直す。この時SNPtoMNPの連続しないものは単なるSNPとして書かれる。
  cat SNPtoMNP.bed | while read line
  do
  CHR=`echo $line | awk '{print $1}'`
  START=`echo $line | awk '{print $2}'`
  END=`echo $line | awk '{print $3}'`
  TARGET_S=`expr $START + 1`
  ID=`bcftools view -H M2.S2.SNV.merged.noMNP.vcf -t ${CHR}:${TARGET_S}-${END} | head -n 1 | awk '{print $3}'`
  REF=`bcftools view -H M2.S2.SNV.merged.noMNP.vcf -t ${CHR}:${TARGET_S}-${END} | awk '{printf $4}'`
  ALT=`bcftools view -H M2.S2.SNV.merged.noMNP.vcf -t ${CHR}:${TARGET_S}-${END} | awk '{printf $5}'`
  LATTER=`bcftools view -H M2.S2.SNV.merged.noMNP.vcf -t ${CHR}:${TARGET_S}-${END} | head -n 1 | awk '{c="";for(i=6;i<=NF;i++) c=c $i"\t"; print c}' | sed 's/\t*$//'`
  echo -e "${CHR}\t${TARGET_S}\t${ID}\t${REF}\t${ALT}\t${LATTER}"
  done > MNP
  # LATTERの最後でなぜか\tが入るのでsedで除去

  ### noMNP.vcfから予め除外するSNPを作ろうと思って下記をしてみたがよくよく考えるとSNPtoMNPと同じものを作ってるだけ、、、
  # cat SNPtoMNP.bed | while read line
  # do
  # CHR=`echo $line | awk '{print $1}'`
  # START=`echo $line | awk '{print $2}'`
  # END=`echo $line | awk '{print $3}'`
  # TARGET_S=`expr $START + 1`
  # bcftools view -H M2.S2.SNV.merged.noMNP.vcf -t ${CHR}:${TARGET_S}-${END}
  # done > SNP
  
  cat <(cat M2.S2.SNV.merged.noMNP.vcf | grep "^#") <(cat <(cat <(bcftools view -H M2.S2.SNV.merged.noMNP.vcf) SNPtoMNP | sort | uniq -u) MNP | sort -V) > M2.S2.SNV.merged.withMNP.vcf

  echo -e "<<<multiallelicのcheck>>>\n<SNV/M2.S2.SNV.extract-multiallelic>\n\n" > ../ToDo

  cd ..
  echo "SNV merge Done"
else
  echo "SNV merge Already Done"
fi


if [ ! -f "SV_VCF/Delly.PASS.vcf" ]; then
  ### SV
  mkdir -p SV_VCF
  cd SV_VCF
  
  # GRIDSS
  # extract o,h PASS
  less ${GRIDSS_VCF} | grep "^#" > GRIDSS.oh.PASS.vcf
  # sortだとhが先に来るので奇数行と偶数行を入れ替えるsed (Yahoo知恵袋のchange.shは上手く動かず)
  bcftools view -H ${GRIDSS_VCF} | awk '{if($3 ~ /o$/ || $3 ~ /h$/) print $0}' | sort -k3 | sed "N; s/\(.*\)\n\(.*\)/\2\\n\1/g" | awk '$7=="PASS"{print}' >> GRIDSS.oh.PASS.vcf
  # extract b PASS
  bcftools view -H ${GRIDSS_VCF} | awk '$3 ~ /b$/' | awk '$7=="PASS"' > gridss.b.pass
### 必ずしも単独BreakpointがINSというわけではない。大きすぎるINSは単独BreakpointとしかShortReadでは検出できない。
#  cat gridss.b.pass | while read line
#  do
#  	CHR=`echo $line | awk '{print $1}'`
#  	POS=`echo $line | awk '{print $2}'`
#  	NAME=`echo $line | awk '{print $3}'`
#  	REF=`echo $line | awk '{print $4}'`
#  	ALT=`echo $line | awk '{print $5}'`
#  	QUAL=`echo $line | awk '{print $6}'`
#  	PASS=`echo $line | awk '{print $7}'`
#  	INFO=`echo $line | awk '{print $8}'`
#  	FORMAT=`echo $line | awk '{print $9}'`
#  	FORMAT_1=`echo $line | awk '{print $10}'`
#  	FORMAT_2=`echo $line | awk '{print $11}'`
#  	ALT_head=`echo $ALT | cut -c 1-1`
#  	if [[ "${ALT_head}" == "." ]]; then
#  		POS_2=`expr "${POS}" + 1`
#  		INS=`echo $ALT | cut -c 2-`
#  		echo -e "$CHR\t$POS\t$NAME-1\t$REF\t$REF$INS[$CHR:$POS_2[\t$QUAL\t$PASS\tMATEID=$NAME-2;$INFO\t$FORMAT\t$FORMAT_1\t$FORMAT_2"
#  		echo -e "$CHR\t$POS_2\t$NAME-2\tN\t]$CHR:$POS]${INS}N\t$QUAL\t$PASS\tMATEID=$NAME-1;$INFO\t$FORMAT\t$FORMAT_1\t$FORMAT_2"
#  	else
#  		POS_2=`expr "${POS}" - 1`
#  		INS=`echo $ALT | rev | cut -c 2- | rev`
#  		echo -e "$CHR\t$POS_2\t$NAME-1\tN\tN$INS[$CHR:$POS[\t$QUAL\t$PASS\tMATEID=$NAME-2;$INFO\t$FORMAT\t$FORMAT_1\t$FORMAT_2"
#  		echo -e "$CHR\t$POS\t$NAME-2\t$REF\t]$CHR:$POS_2]$INS$REF\t$QUAL\t$PASS\tMATEID=$NAME-1;$INFO\t$FORMAT\t$FORMAT_1\t$FORMAT_2"
#  	fi
#  done > gridss.b.pass.vcf
#  cat <(cat GRIDSS.oh.PASS.vcf | grep ^#) <(cat <(cat GRIDSS.oh.PASS.vcf | grep -v ^#) gridss.b.pass.vcf | sort -V) > GRIDSS.ohb.PASS.vcf
  cat <(cat GRIDSS.oh.PASS.vcf | grep ^#) <(cat <(cat GRIDSS.oh.PASS.vcf | grep -v ^#) gridss.b.pass | sort -V) > GRIDSS.ohb.PASS.vcf
  
  # Manta
  less ${Manta_VCF} | grep "^#" > Manta.PASS.vcf
  bcftools view -H ${Manta_VCF} | awk '$7=="PASS"{print}' | while read line
  do
    CHR1=`echo $line | awk '{print $1}'`
    POS1=`echo $line | awk '{print $2}'`
    ID1=`echo $line | awk '{print $3}'`
    REF=`echo $line | awk '{print $4}'`
    REF1=`echo $REF | cut -c 1-1`
    CHR2=`echo $line | awk '{print $1}'`
    INFO=`echo $line | awk '{print $8}'`
    QUAL=`echo $line | awk '{print $6}'`
    FILTER=`echo $line | awk '{print $7}'`
    FORMAT=`echo $line | awk '{print $9}'`
    NORMAL=`echo $line | awk '{print $10}'`
    TUMOR=`echo $line | awk '{print $11}'`
    ALT=`echo $line | awk '{print $5}'`
    if [[ "`echo $line | grep SVTYPE=INS`" ]]; then
      if [[ "${ALT}" == "<INS>" ]]; then
        if [ ! -f "manta.b.pass" ]; then
          > manta.b.pass
        fi
        END=`echo $INFO | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        LEFT_SVINSSEQ=`echo $INFO | awk -F "LEFT_SVINSSEQ=" '{print $2}' | awk -F ";" '{print $1}'`
        RIGHT_SVINSSEQ=`echo $INFO | awk -F "RIGHT_SVINSSEQ=" '{print $2}' | awk -F ";" '{print $1}'`
        ALT1="${REF1}${LEFT_SVINSSEQ}."
        ALT2=".${RIGHT_SVINSSEQ}${REF2}"
        id=`echo ${ID1} | rev | cut -c 2- | rev`
        NEW_ID=`echo "${id}" | sed -e 's/INS/BND:fromINS/g'`
        NEW_ID1=${NEW_ID}0
        NEW_ID2=${NEW_ID}1

        BND1="${CHR1}\t${POS1}\t${NEW_ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${NEW_ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1} | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/END=.*;SVTYPE=/SVTYPE=/g' 
        echo -e ${BND2} | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/END=.*;SVTYPE=/SVTYPE=/g'
        echo -e ${BND1} | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/END=.*;SVTYPE=/SVTYPE=/g' >> manta.b.pass
        echo -e ${BND2} | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/END=.*;SVTYPE=/SVTYPE=/g' >> manta.b.pass
      else
        echo -e "${CHR1}\t${POS1}\t${ID1}\t${REF}\t${ALT}\t${QUAL}\t${FILTER}\t${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      fi
    else
      echo -e "${CHR1}\t${POS1}\t${ID1}\t${REF}\t${ALT}\t${QUAL}\t${FILTER}\t${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
    fi
  done | sort -V >> Manta.PASS.vcf
  
  if [ -f "manta.b.pass" ]; then
    > ../MantaにSingleBreakEndあり
    echo -e "<<<GRIDSSとMantaのSingleBreakEndをcheck!!!!!>>>\n<SV_VCF/gridss.pass.b & manta.pass.b>\n\n" > ../MantaにSingleBreakEndあり
  fi

  # Delly
  # TRAはBNDで表現されているがpairが記述されず、後ろの染色体の表記のみ
  less ${Delly_VCF} | grep "^#" > Delly.BND.PASS.vcf
  bcftools view -H ${Delly_VCF} | awk '$3~/BND/{print}' | awk '$7=="PASS"{print}' >> Delly.BND.PASS.vcf
  
  ## pairを足す。INFOなどはcpして手直しせざるを得ない。
  ## INFOに##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">を追加。INFOにMATEID=を書く。＃これしないとRでpairと認識されず。
  ## IDに.0,.1を追加　INFOのENDもPOS+1　CHR2,POS2も書く
  #if N[[		]]N
  #if N]]		N]]
  #if [[N		[[N
  #if ]]N		N[[
  
  
  bcftools view -H Delly.BND.PASS.vcf | while read line
  do
    CHR0=`echo $line | awk '{print $5}' | awk -F '\\\[|\\\]' '{print $2}' | awk -F ':' '{print $1}'`
    POS0=`echo $line | awk '{print $5}' | awk -F '\\\[|\\\]' '{print $2}' | awk -F ':' '{print $2}'`
    ID0=`echo $line | awk '{print $3".0"}'`
    ID1=`echo $line | awk '{print $3}'`
    MATEID0="MATEID=${ID1}"
    MATEID1="MATEID=${ID0}"
    CHR1=`echo $line | awk '{print $1}'`
    POS1=`echo $line | awk '{print $2}'`
    INFO0_1=`echo $line | awk '{print $8}' | awk -F ";END=" '{print $1}'`
    END0=`expr ${POS0} + 1`
    START0=`expr ${POS0} - 1`
    echo -e "${CHR0}\t${START0}\t${POS0}" > tmp.bed
    REF0=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
    INFO0_2=`echo $line | awk '{print $8}' | awk -F ";PE=" '{print $2}'`
    INFO0="${MATEID0};${INFO0_1};END=${END0};CHR2=${CHR1};POS2=${POS1};PE=${INFO0_2}"
    INFO=`echo $line | awk '{print $8}'`
    INFO1="${MATEID1};${INFO}"
    LATTER=`echo $line | awk '{c="";for(i=9;i<=NF;i++) c=c $i"\t"; print c}'`
    QUAL=`echo $line | awk '{print $6}'`
    FILTER=`echo $line | awk '{print $7}'`
    ALT1=`echo $line | awk '{print $5}'`
    ALT_head=`echo ${ALT1} | cut -c 1-1`
    ALT_tail=`echo ${ALT1} | rev | cut -c 1-1`
    REF1=`echo $line | awk '{print $4}' | cut -c 1-1`
    if [[ "${ALT_head}" == "[" ]]; then
      MATE0="${CHR0}\t${POS0}\t${ID0}\t${REF0}\t[${CHR1}:${POS1}[${REF0}\t${QUAL}\t${FILTER}\t${INFO0}\t${LATTER}"
    elif [[ "${ALT_head}" == "]" ]]; then
      MATE0="${CHR0}\t${POS0}\t${ID0}\t${REF0}\t${REF0}[${CHR1}:${POS1}[\t${QUAL}\t${FILTER}\t${INFO0}\t${LATTER}"
    elif [[ "${ALT_tail}" == "[" ]]; then
      MATE0="${CHR0}\t${POS0}\t${ID0}\t${REF0}\t]${CHR1}:${POS1}]${REF0}\t${QUAL}\t${FILTER}\t${INFO0}\t${LATTER}"
    elif [[ "${ALT_tail}" == "]" ]]; then
      MATE0="${CHR0}\t${POS0}\t${ID0}\t${REF0}\t${REF0}]${CHR1}:${POS1}]\t${QUAL}\t${FILTER}\t${INFO0}\t${LATTER}"
    fi
    echo -e ${MATE0}
    MATE1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${INFO1}\t${LATTER}"
    echo -e ${MATE1}
  done | tr ' ' '\t' > TRA
  
  cat <(head -n 9 Delly.BND.PASS.vcf) <(echo '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">') <(tail -n +10 Delly.BND.PASS.vcf | grep "^#") TRA > Delly.TRA.PASS.vcf
  
  ###TRAを手直しした全体のvcfを作る。
  less Delly.TRA.PASS.vcf | grep "^#" > Delly.PASS.vcf
  bcftools view -H ${Delly_VCF} | awk '$3!~/BND/{print}' | awk '$7=="PASS"{print}' >> Delly.PASS.vcf
  bcftools view -H Delly.TRA.PASS.vcf >> Delly.PASS.vcf

  cd ..
  echo "SV_VCF Done"
else
  echo "SV_VCF Already Done"
fi



if [ ! -f "gr/delly_gr.SV.tsv" ]; then
  ### R gr
  mkdir -p gr
  cd gr
  
  if [ ! -f "delly_gr.tsv" ]; then
    ### INDEL match条件1  maxgap=0, sizemargin=NULL, restrictMarginToSizeMultiple=NULL
    ### INDEL match条件2  maxgap=5, sizemargin=0.25, restrictMarginToSizeMultiple=0.05
    ### SV match条件  maxgap=100, sizemargin=0.25, restrictMarginToSizeMultiple=0.5
    source ~/conda-pack/gridss/bin/activate
    unset PERL5LIB
    unset JAVA_TOOL_OPTIONS
    Rscript /home/n_sasa/data/script/EGA_PCAWG/Merge/Merge_gr.R
    source ~/conda-pack/gatk4.1.9.0/bin/activate
    unset PERL5LIB
    unset JAVA_TOOL_OPTIONS
    
    
    # 1,1はずれているのでN\tをsedで挿入し上書き保存（-i）。
    sed -i -e "1s/^/N\t/g" m2_gr.tsv
    sed -i -e "1s/^/N\t/g" s2_gr.tsv
    sed -i -e "1s/^/N\t/g" gridss_gr.tsv
    sed -i -e "1s/^/N\t/g" gridss_single-brealend_gr.tsv
    sed -i -e "1s/^/N\t/g" manta_gr.tsv
    sed -i -e "1s/^/N\t/g" manta_single-breakend_gr.tsv
    sed -i -e "1s/^/N\t/g" delly_gr.tsv
    sed -i -e "1s/^/N\t/g" delly_single-breakend_gr.tsv

    echo "R Done"
  else
    echo "R Already Done"
  fi
  
  if [ ! -f "delly_gr.SV.tsv" ]; then
    ### INDEL条件1	$20-23
    ### INDEL条件2	$24-27
    ### SV条件	$28-31
    
    ### INDEL条件2とSV条件を使う
    
    
    ### INDELとそれ以外に分ける。matchの条件をINDELとSVで変える。
    # m2,s2はINDELのみだがsize>=50bpも含まれる。
    # m2_gr,s2_grのうちsizeで場合分け
    cat <(head -n 1 m2_gr.tsv) <(tail -n +2 m2_gr.tsv | awk -F "\t" '$15 >-50 && $15<50{print}') > m2_gr.INDEL.tsv
    cat <(head -n 1 m2_gr.tsv) <(tail -n +2 m2_gr.tsv | awk -F "\t" '$15 <=-50 || $15>=50{print}') > m2_gr.SV.tsv
    cat <(head -n 1 s2_gr.tsv) <(tail -n +2 s2_gr.tsv | awk -F "\t" '$15 >-50 && $15<50{print}') > s2_gr.INDEL.tsv
    cat <(head -n 1 s2_gr.tsv) <(tail -n +2 s2_gr.tsv | awk -F "\t" '$15 <=-50 || $15>=50{print}') > s2_gr.SV.tsv
    
    #gridss_gr
    #DELでは$svLen=-{(END-1)-START}+$insLen
    # {(END-1)-START}=$insLen-$svLen <50bp
    #DUP,INVでは$svLen={(END-1)-START}+$insLen
    # {(END-1)-START}=$svLen-$insLen <50bp
    #TRAでは$svLen = NA
    
    #TRAでは$svLen-$insLen=0になってしまう$15=NA（初段もそうなる。）
    #INDELは$svLen-$insLen <50bp && >-50bp　（INSの場合にinsLen>=50でないことも調べておく）
    cat <(head -n 1 gridss_gr.tsv) <(tail -n +2 gridss_gr.tsv | awk -F "\t" '$15!="NA" && $15>-50 && $15<50{print}') > gridss_gr.INDEL.tsv
    cat <(head -n 1 gridss_gr.tsv) <(tail -n +2 gridss_gr.tsv | awk -F "\t" '$15=="NA" || $15<=-50 || $15>=50{print}') > gridss_gr.SV.tsv
    
    #manta_gr
    #DELでは$svLen=-{(END-1)-START}+$insLen
    # {(END-1)-START}=$insLen-$svLen <50bp
    #DUP,INVでは$svLen={(END-1)-START}+$insLen
    # {(END-1)-START}=$svLen-$insLen <50bp
    #TRAでは$svLen = NA
    
    #TRAでは$svLen-$insLe=0になってしまう（初段もそうなる。）
    #INDELは$svLen-$insLen <50bp && >-50bp
    cat <(head -n 1 manta_gr.tsv) <(tail -n +2 manta_gr.tsv | awk -F "\t" '$15!="NA" && $15>-50 && $15<50{print}') > manta_gr.INDEL.tsv
    cat <(head -n 1 manta_gr.tsv) <(tail -n +2 manta_gr.tsv | awk -F "\t" '$15=="NA" || $15<=-50 || $15>=50{print}') > manta_gr.SV.tsv
    
    #delly_gr
    #DEL,DUP,INVでは$svLen={(END-1)-START}+$insLen
    # {(END-1)-START}=$svLen-$insLen <50bp
    #TRAでは$svLen = NA
    
    #TRAでは$svLen-$insLe=0になってしまう（初段もそうなる。）
    #INDELは$svLen-$insLen <50bp
    cat <(head -n 1 delly_gr.tsv) <(tail -n +2 delly_gr.tsv | awk -F "\t" '$15!="NA" && $15<50{print}') > delly_gr.INDEL.tsv
    cat <(head -n 1 delly_gr.tsv) <(tail -n +2 delly_gr.tsv | awk -F "\t" '$15=="NA" || $15>=50{print}') > delly_gr.SV.tsv

    echo "INDEL/SV.tsv Done"
  else
    echo "INDEL/SV.tsv Already Done"
  fi
  
  cd ..
else
  echo "R, INDEL/SV.tsv Already Done"
fi

if [ ! -f "INDEL/m2.s2.indel.vcf" ]; then
    ### INDEL
    #INDEL match条件2は$24-$27
    #2つ以上で共通
    mkdir -p INDEL
    cd INDEL

    echo "###count>1になる変異がないか(overlapしてmatchしていないか)check" >> ../ToDo
    echo "tail -n +2 ../gr/m2_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l" >> ../ToDo
    echo `tail -n +2 ../gr/m2_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l` >> ../ToDo
    echo "tail -n +2 ../gr/s2_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l" >> ../ToDo
    echo `tail -n +2 ../gr/s2_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l` >> ../ToDo
    echo "tail -n +2 ../gr/gridss_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l" >> ../ToDo
    echo `tail -n +2 ../gr/gridss_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l` >> ../ToDo
    echo "tail -n +2 ../gr/manta_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l" >> ../ToDo
    echo `tail -n +2 ../gr/manta_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l` >> ../ToDo
    echo "tail -n +2 ../gr/delly_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l" >> ../ToDo
    echo `tail -n +2 ../gr/delly_gr.INDEL.tsv | awk -F "\t" '"$24">1 || "$25">1 || "$26">1 || "$27">1{print}' | wc -l` >> ../ToDo
    
    #m2基準
    #ID,bp1,2でsortして奇数行のPOSやREF/ALTを抜き出してTARGETにする。
    tail -n +2 ../gr/m2_gr.INDEL.tsv | awk -F "\t" '$24!=0 || $25!=0 || $26!=0 || $27!=0{print}' | sort -k 12,12 -V -k 2,3 | awk 'NR%2!=0' | awk '{print $2"\t"$3"\t"$8","$9}' | bgzip -c > m2.indel.tsv.gz && tabix -s1 -b2 -e2 m2.indel.tsv.gz
    
    #s2基準
    tail -n +2 ../gr/s2_gr.INDEL.tsv | awk -F "\t" '$24==0 && ($25!=0 || $26!=0 || $27!=0){print}' | sort -k 12,12 -V -k 2,3 | awk 'NR%2!=0' | awk '{print $2"\t"$3"\t"$8","$9}' | bgzip -c > s2.indel.tsv.gz && tabix -s1 -b2 -e2 s2.indel.tsv.gz
    
    
    ###m2,s2足し合わせたheaderを作る
    grep "^##" ../M2/M2.PASS.Norm.noMNP.indel.vcf | grep -v "^##contig" > m2.header
    grep "^##" ../S2/S2.INDEL.PASS.Norm.name-changed.vcf | grep -v "^##contig" > s2.header
    contig="grep "^##contig" ../M2/M2.PASS.Norm.noMNP.indel.vcf"
    CHROM="grep "^#CHROM" ../M2/M2.PASS.Norm.noMNP.indel.vcf"
    
    cat <(head -n 1 m2.header) <(cat <(tail -n +2 m2.header) <(tail -n +2 s2.header) | grep "^##[A-Z]" | sort | uniq) <(cat <(tail -n +2 m2.header) <(tail -n +2 s2.header) | grep "^##[a-z]" | sort | uniq) <($contig) <($CHROM) > m2.s2.header
    
    
    ###m2,s2から重複変異を抜き出してまとまる。
    M2_INDEL="bcftools view ../M2/M2.PASS.Norm.noMNP.indel.vcf -H -T m2.indel.tsv.gz"
    S2_INDEL="bcftools view ../S2/S2.INDEL.PASS.Norm.name-changed.vcf -H -T s2.indel.tsv.gz"
    
    cat m2.s2.header <($M2_INDEL) <($S2_INDEL) | bcftools sort > m2.s2.indel.vcf
    
    
    ### if SV callerでのみのINTERSECTIONがある場合
    echo "SV callerでのみのINTERSECTIONがあるかcheck" >> ../ToDo
    echo "gridss_gr.INDEL.tsvでM2,S2は0だがほか(Delly)で0以外" >> ../ToDo
    echo `tail -n +2 ../gr/gridss_gr.INDEL.tsv | awk -F "\t" '$24==0 && $25==0 && ($26!=0 || $27!=0){print}' | sort -k 12,12 -V | sed "N; s/\(.*\)\n\(.*\)/\2\\n\1/g" | wc -l` >> ../ToDo
    #grep "^##" ../SV_VCF/GRIDSS.oh.PASS.vcf | grep -v "^##contig" > gridss.header
    #grep "^##" ../SV_VCF/Manta.PASS.vcf | grep -v "^##contig" > manta.header
    #cat <(head -n 1 m2.s2.header) <(cat <(tail -n +2 m2.s2.header) <(tail -n +2 gridss.header) <(tail -n +2 manta.header) | grep "^##[A-Z]" | sort | uniq) <(cat <(tail -n +2 m2.s2.header) <(tail -n +2 gridss.header) <(tail -n +2 manta.header) | grep "^##[a-z]" | sort | uniq) <($contig) <($CHROM) > m2.s2.g.m.header
    #cat m2.s2.g.m.header <($M2_INDEL) <($S2_INDEL) gridss.indel manta.indel > m2.s2.g.m.indel.vcf

    cd ..  
    echo "INDEL merge Done"
else
  echo "INDEL merge Already Done"
fi
  
  
if [ ! -f "SV/g.m.d.sv.vcf" ]; then
  ### SV merge
  
  #SV match条件は$28-$31
  #2つ以上で共通
  mkdir -p SV
  cd SV
  
  echo "count>1になる変異がないか(overlapしてmatchしていないか)check" >> ../ToDo
  echo "tail -n +2 ../gr/gridss_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l" >> ../ToDo
  echo `tail -n +2 ../gr/gridss_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l` >> ../ToDo
  echo "tail -n +2 ../gr/manta_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l" >> ../ToDo
  echo `tail -n +2 ../gr/manta_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l` >> ../ToDo
  echo "tail -n +2 ../gr/delly_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l" >> ../ToDo
  echo `tail -n +2 ../gr/delly_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l` >> ../ToDo
  echo "tail -n +2 ../gr/m2_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l" >> ../ToDo
  echo `tail -n +2 ../gr/m2_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l` >> ../ToDo
  echo "tail -n +2 ../gr/s2_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l" >> ../ToDo
  echo `tail -n +2 ../gr/s2_gr.SV.tsv | awk -F "\t" '"$28">1 || "$29">1 || "$30">1 || "$31">1{print}' | wc -l` >> ../ToDo
  
  
  if [ ! -f "gridss.sv" ]; then
    # gridss基準
    # ID抽出してVCFから検索する。awk内で変数を使うには-vによる引き渡しが必要。
    # pairの片方だけがHitしてしまうパターンもあるっぽい、、、そのためIDからoh抜いて検索してpair/breakendを採用。
    tail -n +2 ../gr/gridss_gr.SV.tsv | awk -F "\t" '$28!=0 || $29!=0 || $30!=0 || $31!=0{print $12}' | rev | cut -c 2- | rev | sort -u | while read line
    do
    bcftools view -H ../SV_VCF/GRIDSS.ohb.PASS.vcf | awk -v ID="$line" '$3==ID"o" || $3==ID"h" || $3==ID"b"'
    done > gridss.breakpoint

    tail -n +2 ../gr/gridss_single-brealend_gr.tsv | awk -F "\t" '$19!=0 || $20!=0{print $12}' | rev | cut -c 2- | rev | sort -u | while read line
    do
    bcftools view -H ../SV_VCF/GRIDSS.ohb.PASS.vcf | awk -v ID="$line" '$3==ID"o" || $3==ID"h" || $3==ID"b"'
    done > gridss.breakend

    cat gridss.breakpoint gridss.breakend | sort -u | sort -V > gridss.sv
  fi
  
  if [ ! -f "manta.sv" ]; then
    #manta基準
    tail -n +2 ../gr/manta_gr.SV.tsv | awk -F "\t" '$30==0 && ($28!=0 || $29!=0 || $31!=0){print $12,$14}' | while read line
    do
      SVTYPE=$(echo $line | awk '{print $2}')
      sourceID=$(echo $line | awk '{print $1}')
      if [[ "`echo $sourceID | grep DEL`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT1="${REF1}[${CHR2}:${POS2}["
        ALT2="]${CHR1}:${POS1}]${REF2}"
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        NORMAL=`echo $VCF | awk '{print $10}'`
        TUMOR=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $sourceID | grep DUP`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT2="${REF2}[${CHR1}:${POS1}["
        ALT1="]${CHR2}:${POS2}]${REF1}"
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        NORMAL=`echo $VCF | awk '{print $10}'`
        TUMOR=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $sourceID | grep INS`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT=`echo $VCF | awk '{print $5}'`
        ALT2="]${CHR1}:${POS1}]${ALT}${REF2}"
        ALT1="${REF1}${ALT}[${CHR2}:${POS2}["
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        NORMAL=`echo $VCF | awk '{print $10}'`
        TUMOR=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      else
        bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'
      fi
    done | sort -u | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/SVTYPE=DUP/SVTYPE=BND/g' | sed -e 's/SVTYPE=DEL/SVTYPE=BND/g' | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' > manta.breakpoint

    tail -n +2 ../gr/manta_single-breakend_gr.tsv | awk -F "\t" '$19==0 && $20!=0{print $12,$14}' | while read line
    do
      SVTYPE=$(echo $line | awk '{print $2}')
      sourceID=$(echo $line | awk '{print $1}')
      bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'
    done | sort -u > manta.breakend

    cat manta.breakpoint manta.breakend | sort -u | sort -V > manta.sv
  fi

#      LATTER=`echo $VCF | awk '{c="";for(i=9;i<=NF;i++) c=c $i"\t"; print c}'`

  
  if [ ! -f "delly.sv" ]; then
    #delly基準 dellyはSampleの並びが多く、順番も他と違う。
    #INVはstrandを見てALTの書き方を決定する。
    tail -n +2 ../gr/delly_gr.SV.tsv | awk -F "\t" '$30==0 && $31==0 && ($28!=0 || $29!=0){print}' | while read line
    do
    NAME=`echo $line | awk '{print $12}' | awk -F "." '{print $1}'`
    STRAND=`echo $line | awk '{print $6}'`
    if [[ "`echo $NAME | grep DEL`" ]]; then
      VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
      CHR1=`echo $VCF | awk '{print $1}'`
      POS1=`echo $VCF | awk '{print $2}'`
      ID1=`echo $VCF | awk '{print $3}'`
      REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
      CHR2=`echo $VCF | awk '{print $1}'`
      END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
      POS2=`expr ${END} + 1`
      echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
      REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
      ALT1="${REF1}[${CHR2}:${POS2}["
      ALT2="]${CHR1}:${POS1}]${REF2}"
      ID2="${ID1}.1"
      MATEID1="MATEID=${ID2}"
      MATEID2="MATEID=${ID1}"
      QUAL=`echo $VCF | awk '{print $6}'`
      FILTER=`echo $VCF | awk '{print $7}'`
      INFO=`echo $VCF | awk '{print $8}'`
      FORMAT=`echo $VCF | awk '{print $9}'`
      TUMOR=`echo $VCF | awk '{print $10}'`
      NORMAL=`echo $VCF | awk '{print $11}'`
      BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      echo -e ${BND1}
      echo -e ${BND2}
    elif [[ "`echo $NAME | grep DUP`" ]]; then
      VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
      CHR1=`echo $VCF | awk '{print $1}'`
      POS1=`echo $VCF | awk '{print $2}'`
      ID1=`echo $VCF | awk '{print $3}'`
      REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
      CHR2=`echo $VCF | awk '{print $1}'`
      END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
      POS2=`expr ${END} + 1`
      echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
      REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
      ALT2="${REF2}[${CHR1}:${POS1}["
      ALT1="]${CHR2}:${POS2}]${REF1}"
      ID2="${ID1}.1"
      MATEID1="MATEID=${ID2}"
      MATEID2="MATEID=${ID1}"
      QUAL=`echo $VCF | awk '{print $6}'`
      FILTER=`echo $VCF | awk '{print $7}'`
      INFO=`echo $VCF | awk '{print $8}'`
      FORMAT=`echo $VCF | awk '{print $9}'`
      TUMOR=`echo $VCF | awk '{print $10}'`
      NORMAL=`echo $VCF | awk '{print $11}'`
      BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      echo -e ${BND1}
      echo -e ${BND2}
    elif [[ "`echo $NAME | grep INS`" ]]; then
      VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
      CHR1=`echo $VCF | awk '{print $1}'`
      POS1=`echo $VCF | awk '{print $2}'`
      ID1=`echo $VCF | awk '{print $3}'`
#      REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
      REF1=`echo $VCF | awk '{print $4}'`
      CHR2=`echo $VCF | awk '{print $1}'`
      END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
      POS2=`expr ${END} + 1`
      echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
      REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
      ALT=`echo $VCF | awk '{print $8}' | awk -F "CONSENSUS==" '{print $2}' | awk -F ";" '{print $1}'`
      ALT2="]${CHR1}:${POS1}]${ALT}${REF2}"
      ALT1="${REF1}${ALT}[${CHR2}:${POS2}["
      ID2="${ID1}.1"
      MATEID1="MATEID=${ID2}"
      MATEID2="MATEID=${ID1}"
      QUAL=`echo $VCF | awk '{print $6}'`
      FILTER=`echo $VCF | awk '{print $7}'`
      INFO=`echo $VCF | awk '{print $8}'`
      FORMAT=`echo $VCF | awk '{print $9}'`
      TUMOR=`echo $VCF | awk '{print $10}'`
      NORMAL=`echo $VCF | awk '{print $11}'`
      BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      echo -e ${BND1}
      echo -e ${BND2}
    elif [[ "`echo $NAME | grep INV`" ]]; then
      VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
      CHR1=`echo $VCF | awk '{print $1}'`
      POS1=`echo $VCF | awk '{print $2}'`
      REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
      CHR2=`echo $VCF | awk '{print $1}'`
      END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
      POS2=`expr ${END} + 1`
      echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
      REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
      if [[ "${STRAND}" == "+" ]]; then
        ID1=`echo $VCF | awk '{print $3"+"}'`
        ID2="${ID1}.1"
        ALT1="${REF1}]${CHR2}:${POS2}]"
        ALT2="${REF2}]${CHR1}:${POS1}]"
      else
        ID1=`echo $VCF | awk '{print $3"-"}'`
        ID2="${ID1}.1"
        ALT1="[${CHR2}:${POS2}[${REF1}"
        ALT2="[${CHR1}:${POS1}[${REF2}"
      fi
      MATEID1="MATEID=${ID2}"
      MATEID2="MATEID=${ID1}"
      QUAL=`echo $VCF | awk '{print $6}'`
      FILTER=`echo $VCF | awk '{print $7}'`
      INFO=`echo $VCF | awk '{print $8}'`
      FORMAT=`echo $VCF | awk '{print $9}'`
      TUMOR=`echo $VCF | awk '{print $10}'`
      NORMAL=`echo $VCF | awk '{print $11}'`
      BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
      echo -e ${BND1}
      echo -e ${BND2}
    else
      bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID || $3==ID".0"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$10}'
    fi
    done | sort | uniq | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/SVTYPE=DUP/SVTYPE=BND/g' | sed -e 's/SVTYPE=DEL/SVTYPE=BND/g' | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' > delly.sv
  fi
  
  
  #m2基準　基本的にはStrelka2はMantaの<50bpのindelから結果を受けているのでSVには絡んでこない。m2-s2のみがcallするSVはないと思われる。
  
  
  
  ###gridss,manta,dellyのheaderを作る
  grep "^##" ../SV_VCF/GRIDSS.oh.PASS.vcf | grep -v "^##contig" > gridss.header
  grep "^##" ../SV_VCF/Manta.PASS.vcf | grep -v "^##contig" > manta.header
  grep "^##" ../SV_VCF/Delly.PASS.vcf | grep -v "^##contig" > delly.header
  contig="grep "^##contig" ../SV_VCF/GRIDSS.oh.PASS.vcf"
  CHROM="grep "^#CHROM" ../SV_VCF/GRIDSS.oh.PASS.vcf"
  
  cat <(head -n 1 gridss.header)\
   <(cat <(tail -n +2 gridss.header) <(tail -n +2 manta.header) <(tail -n +2 delly.header) | grep "^##[A-Z]" | sort | uniq)\
   <(cat <(tail -n +2 gridss.header) <(tail -n +2 manta.header) <(tail -n +2 delly.header) | grep "^##[a-z]" | sort | uniq)\
   <($contig) <($CHROM) > g.m.d.header
  
  cat g.m.d.header <(cat gridss.sv manta.sv delly.sv | sort -V) > g.m.d.sv.vcf

  cd ..
  echo "SV merge Done"
else
  echo "SV merge Already Done"
fi
  
if [ ! -f "SV-UNION/g.m.d.sv.union.vcf" ]; then  
  ### SV UNION
  # g.m.d.sv.vcfにそれぞれでONLYの変異を足す。Manta、Dellyは変換が必要。DellyのINVについては他と被っていたINVのもう片方のPairも含まれるので注意。
  mkdir -p SV-UNION
  cd SV-UNION

  if [ ! -f "gridss.only.sv" ]; then
    # gridss only
    tail -n +2 ../gr/gridss_gr.SV.tsv | awk -F "\t" '$28==0 && $29==0 && $30==0 && $31==0{print $12}' | rev | cut -c 2- | rev | sort -u | while read line
    do
      bcftools view -H ../SV_VCF/GRIDSS.ohb.PASS.vcf | awk -v ID="$line" '$3==ID"o" || $3==ID"h" || $3==ID"b"{print}'
    done > gridss.only.breakpoint
    cat gridss.only.breakpoint ../SV/gridss.breakend | sort -V | uniq -d > gridss.not-only.breakpoint
    cat gridss.only.breakpoint  gridss.not-only.breakpoint | sort -V | uniq -u > gridss.only.sv
  fi
  
  if [ ! -f "manta.only.sv" ]; then
    #manta only
    tail -n +2 ../gr/manta_gr.SV.tsv | awk -F "\t" '$28==0 && $29==0 && $30==0 && $31==0{print $12,$14}' | while read line
    do
      SVTYPE=$(echo $line | awk '{print $2}')
      sourceID=$(echo $line | awk '{print $1}')
      if [[ "`echo $SVTYPE | grep DEL`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT1="${REF1}[${CHR2}:${POS2}["
        ALT2="]${CHR1}:${POS1}]${REF2}"
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        NORMAL=`echo $VCF | awk '{print $10}'`
        TUMOR=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $SVTYPE | grep DUP`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT2="${REF2}[${CHR1}:${POS1}["
        ALT1="]${CHR2}:${POS2}]${REF1}"
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        NORMAL=`echo $VCF | awk '{print $10}'`
        TUMOR=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $SVTYPE | grep INS`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F "END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT=`echo $VCF | awk '{print $5}'`
        ALT2="]${CHR1}:${POS1}]${ALT}${REF2}"
        ALT1="${REF1}${ALT}[${CHR2}:${POS2}["
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        NORMAL=`echo $VCF | awk '{print $10}'`
        TUMOR=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      else
        bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'
      fi
    done | sort | uniq | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/SVTYPE=DUP/SVTYPE=BND/g' | sed -e 's/SVTYPE=DEL/SVTYPE=BND/g' | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' > manta.only.breakpoint

    tail -n +2 ../gr/manta_single-breakend_gr.tsv | awk -F "\t" '$19!=0 || $20!=0{print $12,$14}' | while read line
    do
      SVTYPE=$(echo $line | awk '{print $2}')
      sourceID=$(echo $line | awk '{print $1}')
      bcftools view -H ../SV_VCF/Manta.PASS.vcf | awk -v ID="$sourceID" '$3==ID{print}'
    done | sort -u > manta.not-only.breakend
    cat manta.only.breakpoint manta.not-only.breakend | sort -V | uniq -d > manta.not-only.breakpoint
    cat manta.only.breakpoint manta.not-only.breakpoint | sort -V | uniq -u > manta.only.sv
  fi
  
  if [ ! -f "delly.only.sv" ]; then
    # delly only
    # 重複が分かっているINVのIDをlistにしておき、それと被っていたらpass
    tail -n +2 ../gr/delly_gr.SV.tsv | awk -F "\t" '$28!=0 || $29!=0 || $30!=0 || $31!=0{print $12}' | grep INV | sort | uniq > delly.dupINV.list
    
    tail -n +2 ../gr/delly_gr.SV.tsv | awk -F "\t" '$28==0 && $29==0 && $30==0 && $31==0{print}' | while read line
    do
    NAME=`echo $line | awk '{print $12}' | awk -F "." '{print $1}'`
    STRAND=`echo $line | awk '{print $6}'`
    if [[ ! "`grep "$NAME" delly.dupINV.list`" ]]; then
      if [[ "`echo $NAME | grep DEL`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT1="${REF1}[${CHR2}:${POS2}["
        ALT2="]${CHR1}:${POS1}]${REF2}"
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        TUMOR=`echo $VCF | awk '{print $10}'`
        NORMAL=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $NAME | grep DUP`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT2="${REF2}[${CHR1}:${POS1}["
        ALT1="]${CHR2}:${POS2}]${REF1}"
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        TUMOR=`echo $VCF | awk '{print $10}'`
        NORMAL=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $NAME | grep INS`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        ID1=`echo $VCF | awk '{print $3}'`
  #      REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        REF1=`echo $VCF | awk '{print $4}'`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        ALT=`echo $VCF | awk '{print $8}' | awk -F "CONSENSUS==" '{print $2}' | awk -F ";" '{print $1}'`
        ALT2="]${CHR1}:${POS1}]${ALT}${REF2}"
        ALT1="${REF1}${ALT}[${CHR2}:${POS2}["
        ID2="${ID1}.1"
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        TUMOR=`echo $VCF | awk '{print $10}'`
        NORMAL=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      elif [[ "`echo $NAME | grep INV`" ]]; then
        VCF=`bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID{print}'`
        CHR1=`echo $VCF | awk '{print $1}'`
        POS1=`echo $VCF | awk '{print $2}'`
        REF1=`echo $VCF | awk '{print $4}' | cut -c 1-1`
        CHR2=`echo $VCF | awk '{print $1}'`
        END=`echo $VCF | awk '{print $8}' | awk -F ";END=" '{print $2}' | awk -F ";" '{print $1}'`
        POS2=`expr ${END} + 1`
        echo -e "${CHR2}\t${END}\t${POS2}" > tmp.bed
        REF2=`cat ${REF_FASTA} | seqkit subseq --bed tmp.bed --quiet | tail -n +2`
        if [[ "${STRAND}" == "+" ]]; then
          ID1=`echo $VCF | awk '{print $3"+"}'`
          ID2="${ID1}.1"
          ALT1="${REF1}]${CHR2}:${POS2}]"
          ALT2="${REF2}]${CHR1}:${POS1}]"
        else
          ID1=`echo $VCF | awk '{print $3"-"}'`
          ID2="${ID1}.1"
          ALT1="[${CHR2}:${POS2}[${REF1}"
          ALT2="[${CHR1}:${POS1}[${REF2}"
        fi
        MATEID1="MATEID=${ID2}"
        MATEID2="MATEID=${ID1}"
        QUAL=`echo $VCF | awk '{print $6}'`
        FILTER=`echo $VCF | awk '{print $7}'`
        INFO=`echo $VCF | awk '{print $8}'`
        FORMAT=`echo $VCF | awk '{print $9}'`
        TUMOR=`echo $VCF | awk '{print $10}'`
        NORMAL=`echo $VCF | awk '{print $11}'`
        BND1="${CHR1}\t${POS1}\t${ID1}\t${REF1}\t${ALT1}\t${QUAL}\t${FILTER}\t${MATEID1};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        BND2="${CHR2}\t${POS2}\t${ID2}\t${REF2}\t${ALT2}\t${QUAL}\t${FILTER}\t${MATEID2};${INFO}\t${FORMAT}\t${NORMAL}\t${TUMOR}"
        echo -e ${BND1}
        echo -e ${BND2}
      else
        bcftools view -H ../SV_VCF/Delly.PASS.vcf | awk -v ID="$NAME" '$3==ID || $3==ID".0"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$10}'
      fi
    fi
    done | sort | uniq | sed -e 's/SVTYPE=INS/SVTYPE=BND/g' | sed -e 's/SVTYPE=DUP/SVTYPE=BND/g' | sed -e 's/SVTYPE=DEL/SVTYPE=BND/g' | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' > delly.only.sv
  fi
  

  # SVで抽出したものと被っているときがあるのでuniq（pairの片割れだけがIntercept判定された場合）
  cat ../SV/g.m.d.header <(cat ../SV/gridss.sv ../SV/manta.sv ../SV/delly.sv gridss.only.sv manta.only.sv delly.only.sv | sort -V | uniq) > g.m.d.sv.union.vcf
  cat ../SV/g.m.d.header <(cat g.m.d.sv.union.vcf ../SV/g.m.d.sv.vcf | sort -V | uniq -u) > g.m.d.sv.low-confidence.vcf
  gatk IndexFeatureFile -I g.m.d.sv.low-confidence.vcf
  
  # <(cat g.m.d.sv.union.vcf ../SV/g.m.d.sv.vcf | sort -V | uniq -u) = <(cat gridss.only.sv manta.only.sv delly.only.sv | sort -V)

  #echo "g.m.d.sv.low-confidence.vcfにALT含まれていないかcheck" >> ../ToDo
  #echo "含まれていたら除いてからgatk IndexFeatureFile" >> ../ToDo

  cd ..
  echo "SV UNION Done"
else
  echo "SV UNION Already Done"
fi



