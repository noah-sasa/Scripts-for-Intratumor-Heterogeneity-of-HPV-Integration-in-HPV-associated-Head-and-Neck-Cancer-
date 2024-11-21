#!/bin/bash

ReadLength=$1
# 101 or 150
InsertSize=$2
# /work23/home/nsasa/data/HPVrOPC_WGS_Analysis/WGS_Data/PostQC/gatk4/DO14282P_multiple_metrics.insert_size_metrics
# MEDIAN_INSERT_SIZE

# host+virus.fa
HOST_VIRUS_FASTA=/home/n_sasa/data/reference/SurVirus/hg38+hpv.fa

survirus_fastq=fastq/results.remapped.t1.txt
survirus_fastq_4000=fastq-4000/results.remapped.t1.txt

survirus_virus_fasta=/home/n_sasa/data/reference/SurVirus/hpv.unaligned.fas




mkdir -p survirus_result
> survirus_result/fastq.vcf

cat $survirus_fastq | while read line
do
	human=$(echo $line | awk '{print $2}')
	human_CHR=$(echo $human | awk -F ":" '{print $1}')
	if [[ "$human" =~ "+" ]]; then
		human_POS=$(echo $human | awk -F ":" '{print $2}' | tr -d "+")
		END5=$(echo $human_POS | awk '{print $1+1}')
		START5=$(echo $human_POS | awk '{print $1+1-5}')
		START_5=$(echo $human_POS | awk '{print $1+1}')
		END_5=$(echo $human_POS | awk '{print $1+6}')
		human_Direction="+"
	else
		human_POS=$(echo $human | awk -F ":" '{print $2}' | tr -d "-" | awk '{print $1+1}')
		START5=$(echo $human_POS | awk '{print $1}')
		END5=$(echo $human_POS | awk '{print $1+5}')
		END_5=$(echo $human_POS | awk '{print $1}')
		START_5=$(echo $human_POS | awk '{print $1-5}')
		human_Direction="-"
	fi
	source ~/conda-pack/survirus/bin/activate
	human_BP_BASE=$(samtools faidx $HOST_VIRUS_FASTA ${human_CHR}:${human_POS}-${human_POS} | tail -n +2 | tr '[:lower:]' '[:upper:]')

	hpv=$(echo $line | awk '{print $3}')
	hpv_CHR=$(echo $hpv | awk -F ":" '{print $1}')
	if [[ "$hpv" =~ "+" ]]; then
		hpv_POS=$(echo $hpv | awk -F ":" '{print $2}' | tr -d "+")
	else
		hpv_POS=$(echo $hpv | awk -F ":" '{print $2}' | tr -d "-" | awk '{print $1+1}')
	fi
	hpv_BP_BASE=$(samtools faidx $survirus_virus_fasta ${hpv_CHR}:${hpv_POS}-${hpv_POS} | tail -n +2 | tr '[:lower:]' '[:upper:]')

	ID=$(echo -e $human"~"$hpv)
	if [[ "$human" =~ "+" ]]; then
		if [[ "$hpv" =~ "+" ]]; then
			human_ALT=$(echo -e $human_BP_BASE"]"$hpv_CHR":"$hpv_POS"]")
			hpv_ALT=$(echo -e $hpv_BP_BASE"]"$human_CHR":"$human_POS"]")
		else
			human_ALT=$(echo -e $human_BP_BASE"["$hpv_CHR":"$hpv_POS"[")
			hpv_ALT=$(echo -e "]"$human_CHR":"$human_POS"]"$hpv_BP_BASE)
		fi
	else
		if [[ "$hpv" =~ "+" ]]; then
			human_ALT=$(echo -e "]"$hpv_CHR":"$hpv_POS"]"$human_BP_BASE)
			hpv_ALT=$(echo -e $hpv_BP_BASE"["$human_CHR":"$human_POS"[")
		else
			human_ALT=$(echo -e "["$hpv_CHR":"$hpv_POS"["$human_BP_BASE)
			hpv_ALT=$(echo -e "["$human_CHR":"$human_POS"["$hpv_BP_BASE)
		fi
	fi
	source ~/conda-pack/pysamstats/bin/activate
	Coverage5bp_Median=$(pysamstats --type coverage fastq/bam_0/retained-pairs.remapped.cs.bam --truncate --no-del --pad --chromosome $human_CHR --start $START5 --end $END5\
             | tail -n +2 | awk '{print $3}' | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}' )
	Split5bp_Max=$(pysamstats --type coverage_ext fastq/bam_0/retained-pairs.remapped.cs.bam --truncate --no-del --pad --chromosome $human_CHR --start $START5 --end $END5\
             | tail -n +2 | awk '{print $9}' | awk 'NR==1 {max=$1} {if($1>max) max=$1} END{print max}')
	Split_5bp_Min=$(pysamstats --type coverage_ext fastq/bam_0/retained-pairs.remapped.cs.bam --truncate --no-del --pad --chromosome $human_CHR --start $START_5 --end $END_5\
             | tail -n +2 | awk '{print $9}' | awk 'BEGIN{min=1000000} {if($1 !="" && min>$1) min=$1} END{print min}')
	SUPPORTING_PAIRS=$(echo $line | awk '{print $4}' | awk -F "=" '{print $2}')
	SPLIT_READS=$(echo $line | awk '{print $5}' | awk -F "=" '{print $2}')
	DISCORDANT_PAIRS=$(echo $SUPPORTING_PAIRS $SPLIT_READS | awk '{print $1-$2}')
	SPLIT_READS_RATIO=$(echo $SPLIT_READS $Coverage5bp_Median $Split5bp_Max $Split_5bp_Min | awk '{print $1/($2-($3-$4)+$1)}')
	DISCORDANT_PAIRS_RATIO=$(echo $DISCORDANT_PAIRS $InsertSize $ReadLength $Coverage5bp_Median $Split5bp_Max $Split_5bp_Min | awk '{print $1*(2*$3/($2-2*$3))/($4-($5-$6)+$1*(2*$3/($2-2*$3)))}')
	SUPPORTING_PAIRS_RATIO=$(echo $SUPPORTING_PAIRS $InsertSize $ReadLength $Coverage5bp_Median $Split5bp_Max $Split_5bp_Min | awk '{print $1*(2*$3/$2)/($4-($5-$6)+$1*(2*$3/$2))}')
	Split5bp_Max_RATIO=$(echo $Split5bp_Max $Split_5bp_Min $Coverage5bp_Median | awk '{print ($1-$2)/$3}')
	VALUE=$(echo $SPLIT_READS_RATIO $Split5bp_Max_RATIO | awk '{if($2 >= $1*1.5){print 0}else{print 1}}')
	MEDIAN=$(echo -e $DISCORDANT_PAIRS_RATIO"\n"$SUPPORTING_PAIRS_RATIO"\n"$Split5bp_Max_RATIO | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}')
	VALUE2=$(echo $SPLIT_READS_RATIO $MEDIAN | awk '{if($2 >= $1*1.5){print 0}else{print 1}}')
	if [[ $SPLIT_READS == 0 ]]; then
		SCORE=0
		if [[ $SUPPORTING_PAIRS == 0 ]]; then
			Predict_VAF=$(echo $Split5bp_Max_RATIO)
		else
			# median
			Predict_VAF=$(echo -e $DISCORDANT_PAIRS_RATIO"\n"$SUPPORTING_PAIRS_RATIO"\n"$Split5bp_Max_RATIO | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}')
		fi
	elif [[ $VALUE == 0 && $VALUE2 == 0 ]]; then
		SCORE=10
		Predict_VAF=$(echo -e $DISCORDANT_PAIRS_RATIO"\n"$SUPPORTING_PAIRS_RATIO"\n"$Split5bp_Max_RATIO | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}')
	else
		SCORE=100
		Predict_VAF=$(echo $SPLIT_READS_RATIO)
	fi

	echo -e $human_CHR"\t"$human_POS"\t"$ID"-h\t"$human_BP_BASE"\t"$human_ALT"\t"$SCORE"\tPASS\tSVTYPE=BND\tSPLIT_READS:SUPPORTING_PAIRS:Coverage5bp_Median:Split5bp_Max:Split_5bp_Min:InsertSize:SPLIT_READS_RATIO:DISCORDANT_PAIRS_RATIO:SUPPORTING_PAIRS_RATIO:Split5bp_Max_RATIO:Predict_VAF\t"$SPLIT_READS":"$SUPPORTING_PAIRS":"$Coverage5bp_Median":"$Split5bp_Max":"$Split_5bp_Min":"$InsertSize":"$SPLIT_READS_RATIO":"$DISCORDANT_PAIRS_RATIO":"$SUPPORTING_PAIRS_RATIO":"$Split5bp_Max_RATIO":"$Predict_VAF\
	 >> survirus_result/fastq.vcf
	echo -e $hpv_CHR"\t"$hpv_POS"\t"$ID"-v\t"$hpv_BP_BASE"\t"$hpv_ALT"\t"$SCORE"\tPASS\tSVTYPE=BND\tSPLIT_READS:SUPPORTING_PAIRS:Coverage5bp_Median:Split5bp_Max:Split_5bp_Min:InsertSize:SPLIT_READS_RATIO:DISCORDANT_PAIRS_RATIO:SUPPORTING_PAIRS_RATIO:Split5bp_Max_RATIO:Predict_VAF\t"$SPLIT_READS":"$SUPPORTING_PAIRS":"$Coverage5bp_Median":"$Split5bp_Max":"$Split_5bp_Min":"$InsertSize":"$SPLIT_READS_RATIO":"$DISCORDANT_PAIRS_RATIO":"$SUPPORTING_PAIRS_RATIO":"$Split5bp_Max_RATIO":"$Predict_VAF\
	 >> survirus_result/fastq.vcf
done




# -4000を直す
> survirus_result/fastq-4000.vcf

cat $survirus_fastq_4000 | while read line
do
	human=$(echo $line | awk '{print $2}')
	human_CHR=$(echo $human | awk -F ":" '{print $1}')
	if [[ "$human" =~ "+" ]]; then
		human_POS=$(echo $human | awk -F ":" '{print $2}' | tr -d "+")
		END5=$(echo $human_POS | awk '{print $1+1}')
		START5=$(echo $human_POS | awk '{print $1+1-5}')
		START_5=$(echo $human_POS | awk '{print $1+1}')
		END_5=$(echo $human_POS | awk '{print $1+6}')
		human_Direction="+"
	else
		human_POS=$(echo $human | awk -F ":" '{print $2}' | tr -d "-" | awk '{print $1+1}')
		START5=$(echo $human_POS | awk '{print $1}')
		END5=$(echo $human_POS | awk '{print $1+5}')
		END_5=$(echo $human_POS | awk '{print $1}')
		START_5=$(echo $human_POS | awk '{print $1-5}')
		human_Direction="-"
	fi
	source ~/conda-pack/survirus/bin/activate
	human_BP_BASE=$(samtools faidx $HOST_VIRUS_FASTA ${human_CHR}:${human_POS}-${human_POS} | tail -n +2 | tr '[:lower:]' '[:upper:]')

	hpv=$(echo $line | awk '{print $3}')
	hpv_CHR=$(echo $hpv | awk -F ":" '{print $1}')
	hpv_len=$(seqkit grep -nrp $hpv_CHR $survirus_virus_fasta | seqkit stats | awk '{print $8}' | tail -n 1 | tr -d ,)
	
	if [[ "$hpv" =~ "+" ]]; then
		hpv_POS=$(echo $hpv | awk -F ":" '{print $2}' | tr -d "+")
	else
		hpv_POS=$(echo $hpv | awk -F ":" '{print $2}' | tr -d "-" | awk '{print $1+1}')
	fi
	hpv_POS2=$(echo $hpv_POS | awk '{print $1+4000}')
	if [[ $hpv_POS2 -gt $hpv_len ]]; then
		hpv_POS3=$(echo $hpv_POS2 $hpv_len | awk '{print $1-$2}')
	else
		hpv_POS3=$(echo $hpv_POS2)
	fi
	hpv_BP_BASE=$(samtools faidx $survirus_virus_fasta ${hpv_CHR}:${hpv_POS3}-${hpv_POS3} | tail -n +2 | tr '[:lower:]' '[:upper:]')

	if [[ "$hpv" =~ "+" ]]; then
		hpv2=$(echo -e $hpv_CHR":+"$hpv_POS3)
	else
		hpv_POS4=$(echo $hpv_POS3 | awk '{print $1-1}')
		hpv2=$(echo -e $hpv_CHR":-"$hpv_POS4)
	fi

	ID=$(echo -e $human"~"$hpv2)
	if [[ "$human" =~ "+" ]]; then
		if [[ "$hpv" =~ "+" ]]; then
			human_ALT=$(echo -e $human_BP_BASE"]"$hpv_CHR":"$hpv_POS3"]")
			hpv_ALT=$(echo -e $hpv_BP_BASE"]"$human_CHR":"$human_POS"]")
		else
			human_ALT=$(echo -e $human_BP_BASE"["$hpv_CHR":"$hpv_POS3"[")
			hpv_ALT=$(echo -e "]"$human_CHR":"$human_POS"]"$hpv_BP_BASE)
		fi
	else
		if [[ "$hpv" =~ "+" ]]; then
			human_ALT=$(echo -e "]"$hpv_CHR":"$hpv_POS3"]"$human_BP_BASE)
			hpv_ALT=$(echo -e $hpv_BP_BASE"["$human_CHR":"$human_POS"[")
		else
			human_ALT=$(echo -e "["$hpv_CHR":"$hpv_POS3"["$human_BP_BASE)
			hpv_ALT=$(echo -e "["$human_CHR":"$human_POS"["$hpv_BP_BASE)
		fi
	fi
	source ~/conda-pack/pysamstats/bin/activate
	Coverage5bp_Median=$(pysamstats --type coverage fastq/bam_0/retained-pairs.remapped.cs.bam --truncate --no-del --pad --chromosome $human_CHR --start $START5 --end $END5\
             | tail -n +2 | awk '{print $3}' | sort -V | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}' )
	Split5bp_Max=$(pysamstats --type coverage_ext fastq/bam_0/retained-pairs.remapped.cs.bam --truncate --no-del --pad --chromosome $human_CHR --start $START5 --end $END5\
             | tail -n +2 | awk '{print $9}' | awk 'NR==1 {max=$1} {if($1>max) max=$1} END{print max}')
	Split_5bp_Min=$(pysamstats --type coverage_ext fastq/bam_0/retained-pairs.remapped.cs.bam --truncate --no-del --pad --chromosome $human_CHR --start $START_5 --end $END_5\
             | tail -n +2 | awk '{print $9}' | awk 'BEGIN{min=1000000} {if($1 !="" && min>$1) min=$1} END{print min}')
	SUPPORTING_PAIRS=$(echo $line | awk '{print $4}' | awk -F "=" '{print $2}')
	SPLIT_READS=$(echo $line | awk '{print $5}' | awk -F "=" '{print $2}')
	DISCORDANT_PAIRS=$(echo $SUPPORTING_PAIRS $SPLIT_READS | awk '{print $1-$2}')
	SPLIT_READS_RATIO=$(echo $SPLIT_READS $Coverage5bp_Median $Split5bp_Max $Split_5bp_Min | awk '{print $1/($2-($3-$4)+$1)}')
	DISCORDANT_PAIRS_RATIO=$(echo $DISCORDANT_PAIRS $InsertSize $ReadLength $Coverage5bp_Median $Split5bp_Max $Split_5bp_Min | awk '{print $1*(2*$3/($2-2*$3))/($4-($5-$6)+$1*(2*$3/($2-2*$3)))}')
	SUPPORTING_PAIRS_RATIO=$(echo $SUPPORTING_PAIRS $InsertSize $ReadLength $Coverage5bp_Median $Split5bp_Max $Split_5bp_Min | awk '{print $1*(2*$3/$2)/($4-($5-$6)+$1*(2*$3/$2))}')
	Split5bp_Max_RATIO=$(echo $Split5bp_Max $Split_5bp_Min $Coverage5bp_Median | awk '{print ($1-$2)/$3}')
	VALUE=$(echo $SPLIT_READS_RATIO $Split5bp_Max_RATIO | awk '{if($2 >= $1*1.5){print 0}else{print 1}}')
	MEDIAN=$(echo -e $DISCORDANT_PAIRS_RATIO"\n"$SUPPORTING_PAIRS_RATIO"\n"$Split5bp_Max_RATIO | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}')
	VALUE2=$(echo $SPLIT_READS_RATIO $MEDIAN | awk '{if($2 >= $1*1.5){print 0}else{print 1}}')
	if [[ $SPLIT_READS == 0 ]]; then
		SCORE=0
		if [[ $SUPPORTING_PAIRS == 0 ]]; then
			Predict_VAF=$(echo $Split5bp_Max_RATIO)
		else
			# median
			Predict_VAF=$(echo -e $DISCORDANT_PAIRS_RATIO"\n"$SUPPORTING_PAIRS_RATIO"\n"$Split5bp_Max_RATIO | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}')
		fi
	elif [[ $VALUE == 0 && $VALUE2 == 0 ]]; then
		SCORE=10
		Predict_VAF=$(echo -e $DISCORDANT_PAIRS_RATIO"\n"$SUPPORTING_PAIRS_RATIO"\n"$Split5bp_Max_RATIO | sort -n | awk '{v[i++]=$1;}END {x=int((i+1)/2); if(x<(i+1)/2) print (v[x-1]+v[x])/2; else print v[x-1];}')
	else
		SCORE=100
		Predict_VAF=$(echo $SPLIT_READS_RATIO)
	fi

	echo -e $human_CHR"\t"$human_POS"\t"$ID"-h\t"$human_BP_BASE"\t"$human_ALT"\t"$SCORE"\tPASS\tSVTYPE=BND\tSPLIT_READS:SUPPORTING_PAIRS:Coverage5bp_Median:Split5bp_Max:Split_5bp_Min:InsertSize:SPLIT_READS_RATIO:DISCORDANT_PAIRS_RATIO:SUPPORTING_PAIRS_RATIO:Split5bp_Max_RATIO:Predict_VAF\t"$SPLIT_READS":"$SUPPORTING_PAIRS":"$Coverage5bp_Median":"$Split5bp_Max":"$Split_5bp_Min":"$InsertSize":"$SPLIT_READS_RATIO":"$DISCORDANT_PAIRS_RATIO":"$SUPPORTING_PAIRS_RATIO":"$Split5bp_Max_RATIO":"$Predict_VAF\
	 >> survirus_result/fastq-4000.vcf
	echo -e $hpv_CHR"\t"$hpv_POS3"\t"$ID"-v\t"$hpv_BP_BASE"\t"$hpv_ALT"\t"$SCORE"\tPASS\tSVTYPE=BND\tSPLIT_READS:SUPPORTING_PAIRS:Coverage5bp_Median:Split5bp_Max:Split_5bp_Min:InsertSize:SPLIT_READS_RATIO:DISCORDANT_PAIRS_RATIO:SUPPORTING_PAIRS_RATIO:Split5bp_Max_RATIO:Predict_VAF\t"$SPLIT_READS":"$SUPPORTING_PAIRS":"$Coverage5bp_Median":"$Split5bp_Max":"$Split_5bp_Min":"$InsertSize":"$SPLIT_READS_RATIO":"$DISCORDANT_PAIRS_RATIO":"$SUPPORTING_PAIRS_RATIO":"$Split5bp_Max_RATIO":"$Predict_VAF\
	 >> survirus_result/fastq-4000.vcf
done

# UNION (SPLITがある方を上にsortしたうえでIDの重複を削除)
cat survirus_result/fastq.vcf survirus_result/fastq-4000.vcf\
 | sort -V -r | awk -F "\t" '!colname[$3]++{print}' | sort -V > survirus_result/union.vcf
