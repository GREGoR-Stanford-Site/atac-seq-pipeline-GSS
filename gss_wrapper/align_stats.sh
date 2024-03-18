#!/bin/bash

# To run: 'bash ./alignstats.sh <Directory> <Cores>'

module load samtools
workdirs=$1 #./gss_atac_persample_workdirs
cores=$2

set -e

qc_dir='./QC'
idxstats="$qc_dir/idxstats"

if [ ! -d "$qc_dir" ]; then
    mkdir -p "$qc_dir"
    echo "Directory created: $qc_dir"
fi

if [ ! -d "$idxstats" ]; then
    mkdir -p "$idxstats"
    echo "Directory created: $idxstats"
fi

align_stats() {
	local bam=$1 # Unfiltered bam
	local idxstats="$2"
	local gss_id=$(basename $bam | sed "s/_S*.*//")
	local primary=${idxstats}/${gss_id}_primary.bam

	# already sorted
	# keep only primary alignments
	samtools view -b -F 0x900 ${bam} -o ${primary}
	# index
	samtools index ${primary}
	samtools idxstats ${primary} > ${idxstats}/${gss_id}_chrinfo.txt
	rm ${primary} ${primary}.bai

	# get counts
	local total=$(awk '{sum+=$3;}END{print sum;}' ${idxstats}/${gss_id}_chrinfo.txt)
	local y=$(grep -E "^chrY[^_]" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3)
	local x=$(grep -E "^chrX" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3)
	local mt=$(grep -E "^chrM" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3)
	local auto=$(grep -E "^chr[0-9]" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3 | awk '{sum+=$1;}END{print sum;}')
	local contig=$(grep -E "_" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3 | awk '{sum+=$1;}END{print sum;}')
	local unmapped=$(grep -E "\*" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 4 | awk '{sum+=$1;}END{print sum;}')

	# get fractions
	local pct_y=$(echo "scale=5; ${y}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_x=$(echo "scale=5; ${x}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_mt=$(echo "scale=5; ${mt}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_auto=$(echo "scale=5; ${auto}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_contig=$(echo "scale=5; ${contig}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_unmapped=$(echo "scale=5; ${unmapped}/(${total}+${unmapped})*100" | bc -l | sed 's/^\./0./')

	# output to file
	echo "gss_id,total_primary_alignments,pct_chrX,pct_chrY,pct_chrM,pct_auto,pct_contig,pct_unmapped" > ${idxstats}/${gss_id}_chrinfo.csv
	echo ${gss_id},${total},${pct_x},${pct_y},${pct_mt},${pct_auto},${pct_contig},${pct_unmapped} >> ${idxstats}/${gss_id}_chrinfo.csv


}

export -f align_stats
#bams=$(ls "${workdirs}"/*/*/*/call-align/*/execution/*.srt.bam) 
bams=$(ls "${workdirs}"/*/*/*/call-align/*/execution/*.srt.bam | head -n 4) #DEBUGGING
parallel --verbose --bar --jobs ${cores} align_stats {} "$idxstats" ::: $bams

# Collapse
rm -f "$qc_dir/merged_chr_info.csv"
cat ${idxstats}/*_chrinfo.csv | grep -v "^gss_id" | sed "1igss_id,total_primary_alignments,pct_chrX,pct_chrY,pct_chrM,pct_auto,pct_contig,pct_unmapped" > "${qc_dir}/merged_chr_info.csv"



