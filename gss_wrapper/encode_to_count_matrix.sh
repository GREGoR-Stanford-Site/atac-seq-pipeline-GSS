#!/bin/bash

# Jeren Olsen
# Adapted from MoTrPAC - Nicole Gay, Anna Scherbina

#Tip: Run with atleast 128GB / cpu. Would get seg faults and other strange errors when running with less. Sorting the tagAlign file before bedtools coverage can reduce this to 64GB.
   # zcat ...tn5.tagAlign.gz | sort -dsk1,1 -k2n,2 -k3nr,3  > ...tn5.tagAlign.sorted
#Example:  'bash ./encode_to_count_matrix.sh 4' (512GB available)

module load bedtools

cores=$1

base_dir=$(pwd)
in_dir="${base_dir}/croo"
out_dir="${base_dir}/postprocessing/merged_peaks"


truncate_script="${base_dir}/truncate_narrowpeak_200bp_summit.py"

if [ ! -d "$out_dir" ]; then
	mkdir -p "$out_dir"
	echo "Directory created: $out_dir"
fi

rm -f "${out_dir}/overlap.optimal_peak.narrowPeak.bed.gz"
touch "${out_dir}/overlap.optimal_peak.narrowPeak.bed.gz"

cd "$in_dir"

#concatenate peaks (narrowpeak.gz)
find . -name "*overlap.optimal_peak.narrowPeak.gz" -exec cat {} >> "${out_dir}/overlap.optimal_peak.narrowPeak.bed.gz" \;
#cat "${in_dir}"/croo/final/peak/*narrowPeak.gz >>"${out_dir}/overlap.optimal_peak.narrowPeak.bed.gz"
echo "Success! concatenated peak files from each sample"

#truncate peaks to 200 bp around summit
python $truncate_script --infile "${out_dir}/overlap.optimal_peak.narrowPeak.bed.gz" --outfile "${out_dir}/overlap.optimal_peak.narrowPeak.200.bed.gz"
echo "Success! finished truncating peaks"

# sort and merge peaks --> master peak file
zcat "${out_dir}/overlap.optimal_peak.narrowPeak.200.bed.gz" | bedtools sort | bedtools merge >"${out_dir}/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed"
echo "Success! Finished sorting and merging"

# intersect with tagalign files
counts_matrix_out="${out_dir}/counts_matrix"
mkdir -p "$counts_matrix_out"

intersect_tag() {
  local TAG=$1
  local GSS_ID

  GSS_ID=$(basename "$TAG" | sed 's/^\([^_]*\).*/\1/')
  echo "$GSS_ID" >"${counts_matrix_out}/counts.${GSS_ID}.txt"
  bedtools coverage -nonamecheck -counts -a "${out_dir}/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed" -b "$TAG" | cut -f4 >>"${counts_matrix_out}/counts.${GSS_ID}.txt"
}

export in_dir
export out_dir
export counts_matrix_out
export -f intersect_tag

# shellcheck disable=SC2046

# parallel --verbose --progress --bar --jobs "$cores" intersect_tag ::: $(ls "${in_dir}"/tagalign/*tagAlign.gz)
echo "Running Parallel intersect tag"
tag_aligns=$(find -name "*.tn5.tagAlign.gz" -exec readlink -f {} \;)
parallel --verbose --progress --bar --jobs "$cores" intersect_tag ::: $tag_aligns
echo "Done "

echo -e $'chrom\tstart\tend' >"${out_dir}/index"
cat "${out_dir}/overlap.optimal_peak.narrowPeak.200.sorted.merged.bed" >>"${out_dir}/index"

#split the results counts matrix by tissue
#to do : reimplement in python 

cd "$counts_matrix_out"
ls ./* | awk -F "." '{print $2}' | awk '{print substr($1,8,2)}' | cut -f1 | sort | uniq >>"${out_dir}/tmp_tids.txt"

while IFS= read -r line; do
  rm -f "${out_dir}/T${line}.atac.counts_matrix.txt"
  rm -f "${out_dir}/T${line}.atac.counts_matrix.txt.gz"
  paste "${out_dir}"/index counts.*"${line}"??.txt >"${out_dir}/T${line}.atac.counts_matrix.txt"
  gzip "${out_dir}/T${line}.atac.counts_matrix.txt"
done <"${out_dir}/tmp_tids.txt"

rm "${out_dir}/tmp_tids.txt"
rm "${out_dir}/index"

echo "Success generating counts matrix"

#Error: line number 10059195 of file /oak/stanford/groups/smontgom/jolsen98/gss_encode_atacseq_dev/gss_atac_persample_workdirs_2_16_24/GSS189457/atac/04e45407-baee-4fc8-893e-c5a5e4bf419a/call-bam2ta/shard-0/execution/glob-199637d3015dccbe277f621a18be9eb4/GSS189457_S34_R1_001.trim.srt.nodup.no_chrM_MT.tn5.tagAlign.gz has 2 fields, but 6 were expected.



