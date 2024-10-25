import sys 
import os
import subprocess
import shutil
import time
from jinja2 import Template
import argparse

# Run normally:
    # python3 encode_atacseq_auto.py \
            # --run_title 'GSS ATACSEQ BATCH 1' \
            # --workdirs_path './gss_atac_persample_workdirs_testing' \
            # --demux_samplesheet_path '/oak/stanford/groups/smontgom/gss_atacseq/testing/2024-01-31_full_demux_sample_sheet.csv' \
            # --data_dir '/oak/stanford/groups/smontgom/gss_atacseq/testing/run1_full_output_NoLaneSplitting/240119_A00509_0854_AHNM2CDMXY/GSS_ATACSEQ_RUN1' 

# Run and overwrite specific IDs:
    # python3 encode_atacseq_atuo.py \
            # --run_title 'GSS ATACSEQ BATCH 1' \
            # --workdirs_path './gss_atac_persample_workdirs_testing' \
            # --demux_samplesheet_path '/oak/stanford/groups/smontgom/gss_atacseq/testing/2024-01-31_full_demux_sample_sheet.csv' \
            # --data_dir '/oak/stanford/groups/smontgom/gss_atacseq/testing/run1_full_output_NoLaneSplitting/240119_A00509_0854_AHNM2CDMXY/GSS_ATACSEQ_RUN1' \
            # --overwrite_ids 'GSS1234567,GSS1234568'

# Run with full overwrite:
    # python3 encode_atacseq_atuo.py \
            # --run_title 'GSS ATACSEQ BATCH 1' \
            # --workdirs_path './gss_atac_persample_workdirs_testing' \
            # --demux_samplesheet_path '/oak/stanford/groups/smontgom/gss_atacseq/testing/2024-01-31_full_demux_sample_sheet.csv' \
            # --data_dir '/oak/stanford/groups/smontgom/gss_atacseq/testing/run1_full_output_NoLaneSplitting/240119_A00509_0854_AHNM2CDMXY/GSS_ATACSEQ_RUN1' \
            # --full_overwrite True


# TODO:
    # Shouldn't be required to provide fastqs dir if running with merged bam  (merged_bam_id flag)

def parse_args():
    parser = argparse.ArgumentParser(prog = 'encode_atacseq_auto.py',
                                     formatter_class = argparse.RawTextHelpFormatter, description = 
                                     '\n'
                                     '-Description:\n'
                                     )

    parser.add_argument('--run_title', help = "Title/name for current run. Ie. 'GSS ATACSEQ BATCH 1'. Alphanumeric characters only" )
    parser.add_argument('--workdirs_path', help =  "Path to directory for work dirs for each sample. This is where all work and outputs will be stored.\nExample: './gss_atacseq_workdirs_MM_DD_YY'")
    parser.add_argument('--demux_samplesheet_path', help='Path to samplesheet used in the demultiplex process. If starting from BAM format, pass "None"')
    parser.add_argument('--data_dir', help='Path to directory containing either 1) all fastqs outputted from demultiplex or 2)  a single raw bam file for each sample')
    parser.add_argument('--genomic_start_format', help = 'data type of files in data_dir, ie, fastq or bam.' ,required = True, default = False, choices = ['fastq','bam'])
    parser.add_argument('--full_overwrite', help='Pass argument with True to overwrite all samples in specified workdir', required=False, default=False, choices = ['True', 'False'])
    parser.add_argument('--overwrite_ids', help='Specific GSS IDs you would like to overwrite/rerun. Only these IDs will be run. Pass as GSS123456,GSS123457... [comma separated, no spaces]', required=False, default=[None])
    parser.add_argument('--merged_bam_id', help = 'Add this argument if running pipeline with single merged BAM. Provide unique identifer as argument that exists in the filename of the merged bam', required=False,default=None)
    parser.add_argument('--batch', help = 'Batch number to add to slurm job name to help distinguish jobs from separate batches', required=True)

    args = parser.parse_args()

    if args.full_overwrite == 'True':
        args.full_overwrite = True

    if any(not c.isalnum() for c in args.run_title if c != " "):
        print("Error, --run_title string must not contain special characters")
        raise Exception

    if eval(args.demux_samplesheet_path) != None and  os.path.exists(args.demux_samplesheet_path) == False:
        print(f"Error, demux_samplesheet_path: '{args.demux_samplesheet_path}' Does Not Exist")
        raise Exception

    if os.path.isdir(args.data_dir) == False:
        print(f"Error, data_dir '{args.data_dir}' Does Not Exist")
        raise Exception

    if args.overwrite_ids != [None] and len(args.overwrite_ids.split(" ")) > 1:
        print("Error, --overwrite_ids argument must not contain spaces")
        raise Exception

    if args.genomic_start_format == 'fastq' and args.merged_bam_id != None:
        print("Cannot pass genomic start type as fastq and provide merged_bam_id")
        raise Exception

    return args

# def check_atac_wdl_no_docker(encode_repo, atac_wdl):
# 	#TODO: Consider always copying (what if changes are made to atac_no_docker.wdl)
# 	if not os.path.exists(atac_wdl):
# 		wdl_script = atac_wdl.split('/')[-1]
# 		shutil.copy2(wdl_script, os.path.join(encode_repo, wdl_script))
# 		print(f" ~ Copied {wdl_script} to {encode_repo}")

# 	return 

def prompt_yes_no(message):
    while True:
        response = input(message + " (y/n): ").strip().lower()
        if response in ['y', 'yes']:
            return True
        elif response in ['n', 'no']:
            return False
        else:
            print("Invalid input. Please enter 'y' or 'n'.")


def get_fastq_GSS_IDs(demux_samplesheet_path):
    ids = []
    with open(demux_samplesheet_path) as f:
        for line in f:
            if line.startswith('GSS') and line[3:9].isdigit():
                ids.append(line[0:9])
    return ids

def generate_samplesheet(template_base_path, gss_id_work_dir, genomic_file_type, file_dict, samplesheet_title, samplesheet_description):
    template_path = template_base_path.format(genomic_file_type.lower()) 
    samplesheet_outpath = gss_id_work_dir+'/atac_samplesheet.json'

    bam_data = None
    fastq_data = None    

    template_content = None
    with open(template_path, "r") as f:
        template_content = f.read()

    # Jinja2 template object
    template = Template(template_content)
    data = {
            "pipeline_type": "atac",
            "genome_tsv": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv",
            "paired_end": "true",
            "auto_detect_adapter": "true",
            "enable_xcor": "true",
            "title": f"{samplesheet_title}",
            "description": f"{samplesheet_description}",
            "align_cpu": 6,
            "align_mem_factor": 1.0,
            "filter_cpu": 6,
            "filter_mem_factor": 1.2
            }

    if genomic_file_type.lower() == 'bam':
        bam_data = {"bam": f"{file_dict['bam_path']}"}
        data.update(bam_data)

    else:
        fastq_data = {
                "fastq_rep1_R1": f"{file_dict['fastq_r1_path']}",
                "fastq_rep1_R2": f"{file_dict['fastq_r2_path']}"
                }
        data.update(fastq_data)

    rendered_template = template.render(data)
    with open(samplesheet_outpath, "w") as f:
        f.write(rendered_template)

    return

#def get_gss_atac_fastq_paths(fastqs_dir, gss_id):
#    fastqs = [item for item in os.listdir(fastqs_dir) if gss_id in item]
#    fastqs = [item for item in fastqs if 'L00' not in item]
#    fastq_r1 = [item for item in fastqs if '_1' in item and item[15:] =='fastp.fastq.gz']
#    fastq_r2 = [item for item in fastqs if '_2' in item and item[15:] =='fastp.fastq.gz']
#    fastq_r1_path = os.path.join(fastqs_dir, fastq_r1)
#    fastq_r2_path = os.path.join(fastqs_dir, fastq_r2)
#
#    return fastq_r1_path, fastq_r2_path

def get_gss_atac_paths(data_dir, gss_id, genomic_file_type):
    file_dict = {'bam_path':None,'fastq_r1_path':None, 'fastq_r2_path':None}
    if genomic_file_type.lower() == 'fastq':
        fastqs = [item for item in os.listdir(data_dir) if gss_id in item and 'fastq' in item]
        fastqs = [item for item in fastqs if 'L00' not in item]
        fastq_r1 = [item for item in fastqs if '_1' in item and item[-15:] =='.fastp.fastq.gz'][0]
        fastq_r2 = [item for item in fastqs if '_2' in item and item[-15:] =='.fastp.fastq.gz'][0]
        fastq_r1_path = os.path.join(data_dir, fastq_r1)
        fastq_r2_path = os.path.join(data_dir, fastq_r2)
        file_dict['fastq_r2_path'] = fastq_r2_path 
        file_dict['fastq_r1_path'] = fastq_r1_path
    elif genomic_file_type.lower() == 'bam':
        bam = [item for item in os.listdir(data_dir) if gss_id in item and 'bam' in item]
        if bam == []:
            return None
        bam = bam[0]
        bam_path = os.path.join(data_dir, bam)
        file_dict['bam_path'] = bam_path
    else:
        raise Exception('Invalid genomic_data_type arg')

    return file_dict


def write_sbatch_script(gss_id, gss_id_work_dir, job_name, partition, atac_run_summary_dir, atac_run_summary_script):

    atac_run_summary_script = "./"+atac_run_summary_script.split("/")[-1]

    run_summary_path = atac_run_summary_dir+f'/{gss_id}_run_summary.csv'

    samplesheet_path = './atac_samplesheet.json' 
    sbatch_script_out_path = gss_id_work_dir+"/launch_encode_atac_sbatch.sh"

    sbatch_script = """#!/bin/bash
#SBATCH --output={gss_id_work_dir}/slurm_%x.%j.out
#SBATCH --job-name={job_name} 
#SBATCH --cpus-per-task=2
#SBATCH --partition={partition}
#SBATCH --account=smontgom
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=6G


#. /home/jolsen98/micromamba/etc/profile.d/conda.sh 
micromamba activate caper 

cd {gss_id_work_dir}

export SINGULARITY_BIND="/oak/stanford/groups/smontgom/jolsen98/atac-seq-pipeline-GSS/src:/mnt/src"

caper run atac.wdl -i {samplesheet_path} \
        --singularity https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/atac-seq-pipeline_v2.2.2.sif \
        --local-loc-dir ./local_loc_dir 
#        --leader-job-name {job_name}

#miniwdl run atac.wdl -i {samplesheet_path} 

python3 {atac_run_summary_script} --gss_id {gss_id} --run_summary_path {run_summary_path}
    """.format(gss_id_work_dir=gss_id_work_dir, samplesheet_path=samplesheet_path, job_name=job_name, partition=partition, gss_id=gss_id, run_summary_path=run_summary_path, atac_run_summary_script=atac_run_summary_script)

    with open(sbatch_script_out_path, "w") as f:
        f.write(sbatch_script)
    f.close()

    os.chmod(sbatch_script_out_path, 0o755)

    return


def run_sample(gss_id_work_dir, job_name, partition):
    sbatch_script = gss_id_work_dir+"/launch_encode_atac_sbatch.sh"
    #subprocess.call(["sbatch", sbatch_script])

    output = subprocess.check_output(["sbatch", sbatch_script])
    batch_job_id = output.decode("utf-8")
    batch_job_id = batch_job_id.split(" ")[-1]

    print(f"Launched:\n\t{sbatch_script}\n\tPartition: {partition}\n\tJobname: {job_name}")

    return batch_job_id

def check_run_completed(gss_id, atac_run_summary_dir):
    run_summary_file = [path for path in os.listdir(atac_run_summary_dir) if gss_id in path]
    if len(run_summary_file) == 0:
        return False
    run_summary_file = run_summary_file[0]
    run_summary_path = os.path.join(atac_run_summary_dir,run_summary_file)

    with open(run_summary_path) as f:
        for line in f:
            if '#Status' in line and 'Succeeded' in line:
                return True
    return False


def write_submission_log_file(outpath, run_metadata):
    with open(outpath, "w") as log_file:
        for key, values in run_metadata.items():
            log_file.write(f"Sample: {key}\n")
            for value in values:
                log_file.write(f"{value}\n")
            log_file.write("\n") 

    print(f"Submission log file written at {outpath}")


def main():

    # Parse args
    args = parse_args()
    run_title = args.run_title
    workdirs = args.workdirs_path
    batch_number = args.batch
    demux_samplesheet_path = args.demux_samplesheet_path
    data_dir = args.data_dir
    genomic_file_type = args.genomic_start_format
    merged_bam_id = args.merged_bam_id
    overwrite = args.full_overwrite
    overwrite_ids = args.overwrite_ids
    if type(overwrite_ids) == str:
        overwrite_ids = overwrite_ids.split(',')

    # Overwrite control
    if overwrite:
        if prompt_yes_no(f"Overwrite = True\n\tAre you sure you want to overwrite everything in {workdirs.split('/')[-1]} ?"):
            print("Continuing...")
        else:
            print("Exiting...")
            return

    # Static Paths
    pwd = os.getcwd()
    atac_samplesheet_template = 'atac_samplesheet_template_{}.json'.format(genomic_file_type)
    encode_repo = '..'
    atac_wdl = encode_repo+"/atac_no_docker.wdl" #Needs to be inside of cloned ENCODE repo: atac-seq-pipeline
    atac_run_summary_script = os.path.join(pwd, "send_run_summary.py") 
    atac_run_summary_dir = os.path.join(pwd, "atac_run_summaries")

    ######## Setup #######

    print(f" ~ Running in {workdirs.split('/')[-1]}")
    if merged_bam_id:
        gss_ids = [merged_bam_id]
    elif overwrite_ids != [None]:
        gss_ids = overwrite_ids
    else:
        gss_ids = get_fastq_GSS_IDs(demux_samplesheet_path)

#	check_atac_wdl_no_docker(encode_repo,atac_wdl)

    if not os.path.isdir(workdirs):
        os.makedirs(workdirs)
    if not os.path.isdir(atac_run_summary_dir):
        os.makedirs(atac_run_summary_dir)

    skip_ids = []
    for gss_id in gss_ids: 
        gss_id_work_dir = os.path.join(workdirs, gss_id)
        run_summary_path = atac_run_summary_dir+f"/{gss_id}_run_summary.csv"
        if os.path.exists(gss_id_work_dir) and (overwrite or gss_id in overwrite_ids):
            print(f'Overwritting {gss_id}')
            shutil.rmtree(gss_id_work_dir)
            if os.path.exists(run_summary_path):
                os.remove(run_summary_path)

        elif os.path.exists(gss_id_work_dir) and check_run_completed(gss_id, atac_run_summary_dir): 
            print(f"Completed run found for {gss_id}")
            skip_ids.append(gss_id)
            continue

        #Given above conditions, a directory exists with a partially completed run or a run summary csv was manually deleted -> clear gssid for new run 
            # Or this is the first time running in specified work directory
        else:
            if not os.path.exists(gss_id_work_dir):
                print(f"Preparing first run for {gss_id}")
            else:
                print(f"Overwritting incomplete run for {gss_id}")
                if os.path.exists(run_summary_path):
                    os.remove(run_summary_path)
                if os.path.exists(gss_id_work_dir):
                    shutil.rmtree(gss_id_work_dir)


        os.makedirs(gss_id_work_dir, exist_ok=True)
        shutil.copy2(atac_wdl, os.path.join(gss_id_work_dir, "atac.wdl"))
        shutil.copy2(atac_run_summary_script, os.path.join(gss_id_work_dir, "send_run_summary.py"))

        # Make samplesheet
        file_dict = get_gss_atac_paths(data_dir, gss_id, genomic_file_type)
        if not file_dict:
            print(f"No matching data file found for id {gss_id}, skipping")
            skip_ids.append(gss_id)
            continue
        generate_samplesheet(atac_samplesheet_template, gss_id_work_dir, genomic_file_type = genomic_file_type, file_dict = file_dict,  samplesheet_title = gss_id, samplesheet_description=run_title)

    gss_ids = [gss_id for gss_id in gss_ids if gss_id not in skip_ids]

    print(f"Samples to be submitted: {len(gss_ids)}")


    ######## Execution ########

    i = 0
    submissions_metadata = {gss_id:[] for gss_id in gss_ids}
    for gss_id in gss_ids:
        i+=1
        #Debugging:
        #if i == 3:
        #    break

        gss_id_work_dir = os.path.join(workdirs, gss_id)

        job_name = f"b{batch_number}_{gss_ids.index(gss_id)+1}_{len(gss_ids)}_atac_leader"
        partition = 'batch'

        write_sbatch_script(gss_id, gss_id_work_dir, job_name, partition, atac_run_summary_dir, atac_run_summary_script)
        batch_job_id = run_sample(gss_id_work_dir, job_name, partition)

        time.sleep(30) if i%10==0 else time.sleep(2)

        submissions_metadata[gss_id] = [f"\tjob_name: {job_name}", f"\tjob_id: {batch_job_id}", f"\tpartition: {partition}", f"\tdirectory: {gss_id_work_dir}"]

    log_outpath = "./sbatch_submission_log.txt"
    write_submission_log_file(log_outpath, submissions_metadata)


    return

if __name__ == '__main__':
    main()










