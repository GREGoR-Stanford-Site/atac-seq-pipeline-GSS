#!/bin/bash

# To run: 'bash ./alignstats.sh <Directory> <Cores>'

# NOTE: Make sure you run this in an environment where mosdepth is installed

module load samtools
workdirs=$1 #./gss_atac_persample_workdirs
cores=$2
cores_sub=$3

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
	# echo "1: ${1}" 
	# echo "2: ${2}"
	# echo "3: ${3}"
	# echo "4: ${4}"
	local bam_unfiltered=$1 # Unfiltered bam
	local bam_filtered=$2
	local idxstats="$3"
	local cores_sub="$4"
	local gss_id=$(basename $bam_unfiltered | sed "s/_S*.*//")
	local primary=${idxstats}/${gss_id}_primary.bam

	# already sorted
	# keep only primary alignments
	samtools view -@ ${cores_sub} -b -F 0x900 ${bam_unfiltered} -o ${primary}
	# index
	samtools index -@ ${cores_sub} ${primary}
	samtools idxstats -@ ${cores_sub} ${primary} > ${idxstats}/${gss_id}_chrinfo.txt
	rm ${primary} ${primary}.bai

	# get counts (primary alignments)
	local total=$(awk '{sum+=$3;}END{print sum;}' ${idxstats}/${gss_id}_chrinfo.txt)
	local y=$(grep -E "^chrY[^_]" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3)
	local x=$(grep -E "^chrX" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3)
	local mt=$(grep -E "^chrM" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3)
	local auto=$(grep -E "^chr[0-9]" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3 | awk '{sum+=$1;}END{print sum;}')
	local contig=$(grep -E "_" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 3 | awk '{sum+=$1;}END{print sum;}')
	local unmapped=$(grep -E "\*" ${idxstats}/${gss_id}_chrinfo.txt | cut -f 4 | awk '{sum+=$1;}END{print sum;}')

	# get fractions (primary alignments)
	local pct_y=$(echo "scale=5; ${y}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_x=$(echo "scale=5; ${x}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_mt=$(echo "scale=5; ${mt}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_auto=$(echo "scale=5; ${auto}/${total}*100" | bc -l | sed 's/^\./0./')
	local pct_contig=$(echo "scale=5; ${contig}/${total}*100" | bc -l | sed 's/^\./0./')
	#local pct_unmapped=$(echo "scale=5; ${unmapped}/(${total}+${unmapped})*100" | bc -l | sed 's/^\./0./')

	# Percent unique, secondary and unmapped alignments (unfiltered vs primary)
	samtools idxstats -@ ${cores_sub} ${bam_unfiltered} > ${idxstats}/${gss_id}_all_chrinfo.txt
	local total_all=$(awk '{sum+=$3;}END{print sum;}' ${idxstats}/${gss_id}_all_chrinfo.txt)

	local pct_primary_alignments_unfiltered=$(echo "scale=5; $total / $total_all * 100" | bc)
	local pct_multimapped_unfiltered=$(echo "scale=5; ($total_all - $total) / $total_all * 100" | bc)
	local pct_unmapped_unfiltered=$(echo "scale=5; ${unmapped}/${total_all}*100" | bc -l | sed 's/^\./0./')

	# Average coverage unfiltered bam
	mosdepth --threads ${cores_sub} --no-per-base ${gss_id} ${bam_unfiltered}
	local avg_coverage_unfiltered=$(grep 'total' ${gss_id}.mosdepth.summary.txt | awk '{print $4}')
	mv ${gss_id}.mosdepth.summary.txt ${idxstats}/${gss_id}.mosdepth.summary.txt
	mv ${gss_id}.mosdepth.global.dist.txt ${idxstats}/${gss_id}.mosdepth.global.dist.txt

	# Average coverage filtered bam
	mosdepth --threads ${cores_sub} --no-per-base ${gss_id} ${bam_filtered}
	local avg_coverage_filtered=$(grep 'total' ${gss_id}.mosdepth.summary.txt | awk '{print $4}')
	mv ${gss_id}.mosdepth.summary.txt ${idxstats}/${gss_id}_filtered.mosdepth.summary.txt
	mv ${gss_id}.mosdepth.global.dist.txt ${idxstats}/${gss_id}_filtered.mosdepth.global.dist.txt

	# output to file
	echo "gss_id,total_primary_alignments,pct_chrX_primary,pct_chrY_primary,pct_chrM_primary,pct_auto_primary,pct_contig_primary,pct_primary_alignments_unfiltered,pct_multimapped_unfiltered,pct_unmapped_unfiltered,avg_coverage_unfiltered,avg_coverage_filtered" > ${idxstats}/${gss_id}_chrinfo.csv
	echo ${gss_id},${total},${pct_x},${pct_y},${pct_mt},${pct_auto},${pct_contig},${pct_primary_alignments_unfiltered},${pct_multimapped_unfiltered},${pct_unmapped_unfiltered},${avg_coverage_unfiltered},${avg_coverage_filtered} >> ${idxstats}/${gss_id}_chrinfo.csv

}

export -f align_stats

#bams=$(ls "${workdirs}"/*/*/*/call-align/*/execution/*.srt.bam) 
echo -n "Finding call-align bams... "
bams_unfiltered=$(ls "${workdirs}"/*/*/*/call-align/*/execution/*.srt.bam)
echo " ✓"
echo -n "Finding call-filter bams..."
bams_filtered=$(ls "${workdirs}"/*/*/*/call-filter/*/execution/*.bam)
echo " ✓"
echo "Running align_stats()"

echo "jobs: ${cores}"
echo "cores sub: ${cores_sub}"

#parallel --verbose --bar --jobs ${cores} align_stats {} {} "$idxstats" ${cores_sub} ::: $bams_unfiltered $bams_filtered #:::+ $bams_filtered
#parallel --verbose --bar --jobs ${cores} align_stats {1} {2} "$idxstats" ${cores_sub} ::: $bams_unfiltered ::: $bams_filtered
parallel --verbose --bar --jobs ${cores} align_stats {1} {2} "$idxstats" ${cores_sub} :::+ $bams_unfiltered :::+ $bams_filtered

# Collapse
rm -f "$qc_dir/merged_chr_info.csv"
cat ${idxstats}/*_chrinfo.csv | grep -v "^gss_id" | sed "1igss_id,total_primary_alignments,pct_chrX_primary,pct_chrY_primary,pct_chrM_primary,pct_auto_primary,pct_contig_primary,pct_primary_alignments_unfiltered,pct_multimapped_unfiltered,pct_unmapped_unfiltered,avg_coverage_unfiltered,avg_coverage_filtered" > "${qc_dir}/merged_chr_info.csv"



