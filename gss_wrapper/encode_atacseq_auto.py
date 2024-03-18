import sys 
import os
import subprocess
import shutil
import time
from jinja2 import Template

# Run normally:
	# python3 encode_atacseq_atuo.py 

# Run and overwrite specific IDs:
	# python3 encode_atacseq_atuo.py GSS123456,GSS135790

# ! Edit Settings section of main function !

def check_atac_wdl_no_docker(encode_repo, atac_wdl):
	#TODO: Consider always copying (what if changes are made to atac_no_docker.wdl)
	if not os.path.exists(atac_wdl):
		wdl_script = atac_wdl.split('/')[-1]
		shutil.copy2(wdl_script, os.path.join(encode_repo, wdl_script))
		print(f" ~ Copied {wdl_script} to {encode_repo}")

	return 

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

def generate_samplesheet(template_path, gss_id_work_dir, fastq_r1_path, fastq_r2_path, samplesheet_title, samplesheet_description):

	samplesheet_outpath = gss_id_work_dir+'/atac_samplesheet.json'

	template_content = None
	with open(template_path, "r") as f:
		template_content = f.read()

	# Jinja2 template object
	template = Template(template_content)
	data = {
		"pipeline_type": "atac",
		"genome_tsv": "https://storage.googleapis.com/encode-pipeline-genome-data/genome_tsv/v4/hg38.tsv",
		"fastq_rep1_R1": f"{fastq_r1_path}",
		"fastq_rep1_R2": f"{fastq_r2_path}",
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

	rendered_template = template.render(data)
	with open(samplesheet_outpath, "w") as f:
		f.write(rendered_template)

	return

def get_gss_atac_fastq_paths(fastqs_dir, gss_id):
	#TODO: improve filtering to remove possible edge cases (regex)
	fastqs = [item for item in os.listdir(fastqs_dir) if gss_id in item]
	fastqs = [item for item in fastqs if 'L00' not in item]
	fastq_r1 = [item for item in fastqs if 'R1' in item][0]
	fastq_r2 = [item for item in fastqs if 'R2' in item][0]
	fastq_r1_path = os.path.join(fastqs_dir, fastq_r1)
	fastq_r2_path = os.path.join(fastqs_dir, fastq_r2)

	return fastq_r1_path, fastq_r2_path


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
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=6G


. /home/jolsen98/micromamba/etc/profile.d/conda.sh 
conda activate encode_env

cd {gss_id_work_dir}

caper run atac.wdl -i {samplesheet_path} \
--singularity https://encode-pipeline-singularity-image.s3.us-west-2.amazonaws.com/atac-seq-pipeline_v2.2.2.sif \
--local-loc-dir ./local_loc_dir \
--leader-job-name {job_name}

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

	# args

	# Can write as only cmdln argument GSSIDs you want overwritten. Example: 'python3 encode_atacseq_atuo.py GSS123456,GSS135790'
	overwrite_ids = [None]
	if len(sys.argv) > 1:
		overwrite_ids = sys.argv[1]
		overwrite_ids = overwrite_ids.split(',') 

	# overwrite = False: resumes processing by default, does not overwrite samples that have finished
	# overwrite = True: Everything is deleted and rerun
	overwrite = False

	pwd = os.getcwd()

	# SET THESE EACH RUN (of 64 samples)
	run_description = "GSS ATACSEQ BATCH 1"
	demux_samplesheet_path = '/oak/stanford/groups/smontgom/gss_atacseq/testing/2024-01-31_full_demux_sample_sheet.csv'
	atac_samplesheet_template = 'atac_samplesheet_template.json'
	fastqs_dir = '/oak/stanford/groups/smontgom/gss_atacseq/testing/run1_full_output_NoLaneSplitting/240119_A00509_0854_AHNM2CDMXY/GSS_ATACSEQ_RUN1'
	encode_repo = '..'
	atac_wdl = encode_repo+"/atac_no_docker.wdl" #Needs to be inside of cloned ENCODE repo: atac-seq-pipeline 
	workdirs = os.path.join(pwd, "gss_atac_persample_workdirs_testing")
	# END OF SETTINGS

	print(f" ~ Running in {workdirs.split('/')[-1]}")

	atac_run_summary_script = os.path.join(pwd, "send_run_summary.py")
	atac_run_summary_dir = os.path.join(pwd, "atac_run_summaries")


	# Overwrite control

	if overwrite:
		if prompt_yes_no(f"Overwrite = True\n\tAre you sure you want to overwrite everything in {workdirs.split('/')[-1]} ?"):
			print("Continuing...")
		else:
			print("Exiting...")
			return


	######## Setup #######

	gss_ids = get_fastq_GSS_IDs(demux_samplesheet_path)

	check_atac_wdl_no_docker(encode_repo,atac_wdl)

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
		fastq_r1_path, fastq_r2_path = get_gss_atac_fastq_paths(fastqs_dir, gss_id)
		generate_samplesheet(atac_samplesheet_template, gss_id_work_dir, fastq_r1_path, fastq_r2_path, samplesheet_title = gss_id, samplesheet_description=run_description)

	gss_ids = [gss_id for gss_id in gss_ids if gss_id not in skip_ids]

	print(f"Samples to be submitted: {len(gss_ids)}")

	# Execution
	i = 0
	submissions_metadata = {gss_id:[] for gss_id in gss_ids}
	for gss_id in gss_ids:
		i+=1

		gss_id_work_dir = os.path.join(workdirs, gss_id)

		job_name = f"{gss_ids.index(gss_id)+1}_{len(gss_ids)}_atac_leader"
		partition = 'batch'

		write_sbatch_script(gss_id, gss_id_work_dir, job_name, partition, atac_run_summary_dir, atac_run_summary_script)
		batch_job_id = run_sample(gss_id_work_dir, job_name, partition)
		if i%10 == 0:
			time.sleep(2)

		time.sleep(30) if i%10==0 else time.sleep(2)

		submissions_metadata[gss_id] = [f"\tjob_name: {job_name}", f"\tjob_id: {batch_job_id}", f"\tpartition: {partition}", f"\tdirectory: {gss_id_work_dir}"]

	log_outpath = "./sbatch_submission_log.txt"
	write_submission_log_file(log_outpath, submissions_metadata)


	return

if __name__ == '__main__':
	main()










