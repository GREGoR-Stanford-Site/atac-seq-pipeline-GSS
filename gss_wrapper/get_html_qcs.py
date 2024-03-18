import os 
import shutil
import sys

# ! Pass work directory as only argument !

work_dir_path = sys.argv[1]

basepath = os.getcwd()
out_dir = basepath+"/QC/qc_html"
outpath = out_dir+"/{}_qc.html"

if not os.path.exists(out_dir):
	os.mkdir(out_dir)

sample_dirs = [os.path.join(work_dir_path, path) for path in os.listdir(work_dir_path)]
sample_dirs = [path for path in sample_dirs if "GSS" in path]
copied = 0
for path in sample_dirs:
	process_dir = os.listdir(path+"/atac")[0]
	full_path = path+f"/atac/{process_dir}/call-qc_report/execution/qc.html"
	if os.path.exists(full_path):
		gss_id = [item for item in full_path.split('/') if 'GSS' in item][0]
		shutil.copy2(full_path, outpath.format(gss_id))
		copied+=1
	else:
		print(f"File does not exist:\n\t{full_path}")

print(f"{copied} qc.htmls copied to {out_dir}")

