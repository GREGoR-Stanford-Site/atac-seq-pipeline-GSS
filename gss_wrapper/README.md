# !Rough draft / WIP!

# GSS ATAC-Seq Pipeline README

This document outlines the steps involved in the GSS ATAC-Seq pipeline from raw BCLs to QC'd bams and peak files.

## Steps:

### 0. Data Acquisition

- BCL data for a new run is deposited into `/oak/stanford/groups/smontgom/gss_atacseq`.
  - Example: `/oak/stanford/groups/smontgom/gss_atacseq/240119_A00509_0854_AHNM2CDMXY`.

### 1. Prepare Sample Sheet

- Acquire clean samplesheet, tracking sheet, and possibly dual index primer list if tracking sheet is incomplete.
- Create samplesheet with `ATAC-Seq SampleSheet.ipynb`.
  - Tracking Sheet: [Google Sheets Link](https://docs.google.com/spreadsheets/d/113Iaqh3nERamkSqu-e5iAnM3s8UKfzIlcdwLMGf2jFI/edit?pli=1#gid=2065782038) (or ask Kevin).

### 2. Demultiplexing

- Create a directory for this run (Example: `gss_atacseq/GSS_RUN1_DEMUX`).
- Run the (modified) nextflow demultiplex pipeline:

```bash
nextflow run https://github.com/JamWithBread/nf_demultiplex \
-profile singularity \
--demultiplexer bcl2fastq \
--input nf_full_SampleSheet.csv \
--outdir ./run1_full_output_NoLaneSplitting \
-c demux_process_high_override.config
```

> **i)** Modified to force bcl2fastq module to not use lane splitting, passing no lane splitting in samplesheet doesn't work
> 
> **ii)** nf_full_SampleSheet.csv -> See nf-core for syntax: https://nf-co.re/demultiplex/1.4.1/docs/usage
> 
> Example:
> 
> id,samplesheet,lane,flowcell
> 
> 240119_A00509_0854_AHNM2CDMXY,/oak/stanford/groups/smontgom/gss_atacseq/testing/2024-01-31_full_demux_sample_sheet.csv,,/oak/stanford/groups/smontgom/gss_atacseq/240119_A00509_0854_AHNM2CDMXY
> 
>> 2024-01-31_full_demux_sample_sheet.csv was created in step 1
>
> **iii)** -c demux_process_high_override.config. Increase resources for bcl2fastq step. Also limits nextflow executor job submission rate (Pipeline crashes when all jobs are submitted at once)

### 3. Ask someone in wetlab for wetlab QCs.
- Examples:
	- GREGOR_ATAC_Libraries_BA.pdf (gregor_atac_pbmc_2100 expert_High Sensitivity DNA Assay_DE04103588_2023-12-04_15-13-23.xad)
	- 11-30-23_Library_Summary_PBMC_ATAC_GREGoR.xlsx  

### 4. Run ENCODE ATAC-Seq Pipeline & QC.
1. **Create New Directory and Setup Paths**
   
   - Clone this repo; [GSS_ATAC-Seq](https://github.com/GREGoR-Stanford-Site/GSS_ATAC-Seq)
   - Rename folder of cloned repo to something indicating current gss atac run, ie GSS_RUN1_encode_atacseq
   - cd GSS_RUN1_encode_atacseq
   - Clone the [Encode ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).
   - Setup paths in `encode_atacseq_auto.py`.

   **Clean Directory Structure:**
```bash
GSS_RUN1_encode_atacseq
	├── align_stats2.sh
	├── atac_no_docker.wdl 
	├── atac_samplesheet_template.json
	├── atac-seq-pipeline -> from ENCODE
	├── croo2qc2tsv.sh
	├── encode_atacseq_auto.py
	├── encode_to_count_matrix.sh
	├── get_html_qcs.py
	├── send_run_summary.py
	└── truncate_narrowpeak_200bp_summit.py
```

2. **Execution Steps**

- Run `encode_atacseq_auto.py` (takes about 24 hours)
- Run `croo2qc2tsv.sh` to get aggregated QC TSV.
  ```
  ./croo2qc2tsv.sh gss_atac_persample_workdirs_2_16_24
  ```
- Run alignment stats.
  ```
  bash ./align_stats2.sh gss_atac_persample_workdirs_2_16_24 64
  ```
- Run `get_html_qcs.py`. Puts all QC HTMLs in one place
  ```
  python3 get_html_qcs.py gss_atac_persample_workdirs_2_16_24
  ```
- (OPTIONAL) Run `encode_to_count_matrix.sh` for peak merging. Read the top of the file for memory info
  > Used for downstream analysis (ie corroboration with RNASeq). Could get creative with QC.
  ```
  bash ./encode_to_count_matrix.sh 4
  ```
- Get all frag length distributions into one place.
  > flip through them in file viewer and make note of general pattern of NFR/mono-nuc/di-nuc densitites for given run (and major outliers)
  ```
  find ./<workdir> -name "*no_chrM_MT.fraglen_dist.png" | grep "call-qc_report" | xargs -I {} cp {} ./QC/frag_len_dists
  ```
  
  **After All is Said and Done:**
```bash
GSS_RUN1_encode_atacseq
	├── align_stats2.sh
	├── atac_no_docker.wdl
	├── atac_run_summaries 
	│ 	└── (64 run summaries)
	├── atac_samplesheet_template.json
	├── atac-seq-pipeline
	├── croo 
	│       └── (64 croo dirs)
	├── croo2qc2tsv.sh
	├── encode_atacseq_auto.py
	├── encode_to_count_matrix.sh
	├── get_html_qcs.py
	├── gss_atac_persample_workdirs
	│ 	└── (64 sample directories -> **all ATAC outputs**)
	├── postprocessing
	│ 	└── (Merged peak sets, peak counts matrix)
	├── QC
	│     ├── (64 qc_htmls)
	│     ├── (64 frag length dists)
	│     ├── qc2tsv_out.tsv (aggregated QC)
	│     └── merged_chr_info.csv (align stats)
	├── sbatch_submission_log.txt 
	├── send_run_summary.py
	└── truncate_narrowpeak_200bp_summit.py
```

3. **Post-Processing and Summary**

- [ ] TODO: Do something QC related with the peaks count matrix.
- [ ] TODO: Get set of markers (ie. immune genes) to verify expected TSS sites / presence of peak in bigwig:peak.bed file pairs for each sample (currently doing manually in IGV)
- Generate summary QC report. For now, load `QC/merged_chr_info.csv` and `QC/qc2tsv_out.tsv` into notebook `ATAC_QC.ipynb`.
  - [ ] TODO: Run notebook as Python script with HTML output.
.
