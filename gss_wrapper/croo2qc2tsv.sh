#!/bin/bash

work_dirs=$1

qc_dir='./QC'
croo_dir='./croo'
echo "Searching for metadata.json files ..."
files=$(find "${work_dirs}" -mindepth 3 -maxdepth 5 -type f -wholename "*/GSS*/atac/*/metadata.json" -not -path '*/local_loc_dir/*')

files_list='./files_list.txt'
echo "$files" > "$files_list"

qc_json_list='./qc_json_list.txt'
touch $qc_json_list

croo_out=$croo_dir

if [ ! -d "$qc_dir" ]; then
    # Create directory if it does not exist
    mkdir -p "$qc_dir"
    echo "Directory created: $qc_dir"
fi

if [ ! -d "$croo_out" ]; then
    # Create directory if it does not exist
    mkdir -p "$croo_out"
    echo "Directory created: $croo_out"
fi

# TODO: check that unique IDs do not appear more than once

while IFS= read -r file; do
    echo "Processing file: $file"
    gss_id=$(echo "$file" | cut -d'/' -f2)
    echo "gss_id: $gss_id"
    croo --method link "$file" --out-dir "$croo_out/croo_$gss_id"
    echo "$croo_out/croo_${gss_id}/qc/qc.json" >> "$qc_json_list"
done < "$files_list"

rm -f "$qc_dir/qc2tsv_out.tsv"
qc2tsv --file "$qc_json_list" --collapse-header > "$qc_dir/qc2tsv_out.tsv"

rm -f $qc_json_list
rm -f $files_list