import pandas as pd
import os
import argparse
from datetime import datetime
from datetime import date
import json

def parse_args():
    parser = argparse.ArgumentParser(prog = 'send_run_summary.py',
        formatter_class = argparse.RawTextHelpFormatter, description = 
        '\n'
        '-Description:\n'
        )

    parser.add_argument('--gss_id', help = 'gss_id')
    parser.add_argument('--run_summary_path', help =  'Path to output directory of run summary files')

    args = parser.parse_args()

    return args


def locate_metadata_path():

    process_dir = os.listdir('./atac') if os.path.exists('./atac') else False
    if process_dir == False:
        return False
    if len(process_dir) > 1: 
        print(f'Correct process dir is ambiguos, multiple exist: {process_dir}. Find better way to locate atac metadata.json file ')
        raise Error
    process_dir = process_dir[0]
    metadata_path = f"./atac/{process_dir}/metadata.json"

    return metadata_path

def get_time_dif(timestamp1_str,timestamp2_str):
    timestamp1_str = timestamp1_str.split('.')[0]
    timestamp2_str = timestamp2_str.split('.')[0]

    timestamp1 = datetime.strptime(timestamp1_str, '%Y-%m-%dT%H:%M:%S')
    timestamp2 = datetime.strptime(timestamp2_str, '%Y-%m-%dT%H:%M:%S')
    
    time_difference = timestamp2 - timestamp1

    return str(time_difference)

def get_pipeline_status(gss_id, metadata_path, outpath):

    # Check if metadata.json exists
    if not metadata_path or not os.path.exists(metadata_path):
        with open(outpath, 'w') as f:
            f.write(f'Pipeline Failed, {metadata_path} does not exist')
        f.close()
        return
    
    with open(metadata_path, 'r') as f:
        data = json.load(f)
    
    # Sample Level Summary
    status = data['status'] if 'status' in data.keys() else None 
    start = data['start'] if 'start' in data.keys() else None
    end = data['end'] if 'end' in data.keys() else None
    elapsed = get_time_dif(start, end) if start and end else None
    jobs = data['calls'].keys() if 'calls' in data.keys() else None
    
    # Job Level Summary
    cols = ['job', 'status', 'start', 'end', 'elapsed']
    job_report = pd.DataFrame(columns = cols)
    job_report['job'] = list(data['calls'].keys())
    if jobs != None:
        for i, row in job_report.iterrows():
            job_info = data['calls'][row['job']][0]
            job_report.at[i, 'status'] = job_info['executionStatus']
            job_report.at[i, 'start'] = job_info['start']
            job_report.at[i, 'end'] = job_info['end']

            diff = get_time_dif(job_info['start'], job_info['end'])
            job_report.at[i, 'elapsed'] = diff

        job_report = job_report.sort_values(by='start')

    # Write summaries to head directory
    with open(outpath, 'w') as f:
        f.write(f'#ENCODE ATAC RUN SUMMARY {gss_id}\n')
        f.write(f'#Status: {status}\n')
        f.write(f'#Start time: {start}\n')
        f.write(f'#End time: {end}\n')
        f.write(f'#Elapsed: {elapsed}\n')
    f.close()
    job_report.to_csv(outpath, mode = 'a', index=False)

    return


def main():

    args = parse_args()
    outpath = args.run_summary_path
    gss_id = args.gss_id


    metadata_path = locate_metadata_path()
    get_pipeline_status(gss_id, metadata_path, outpath)

    return



if __name__ == '__main__':
    main()



















