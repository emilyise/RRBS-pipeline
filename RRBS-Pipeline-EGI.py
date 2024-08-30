import subprocess
import time
import sys
import os

def submit_job(script_path, project_dir):
    """Submit a job to SLURM using sbatch and include the project directory as an argument."""
    result = subprocess.run(['sbatch', script_path, project_dir], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    job_id = None
    if result.returncode == 0:
        output = result.stdout.decode('utf-8')
        print(f"Job submitted: {output.strip()}")
        job_id = output.strip().split()[-1]
    else:
        print(f"Error submitting job: {result.stderr.decode('utf-8')}")
    return job_id

def wait_for_job(job_id):
    """Wait for the job to complete."""
    while True:
        result = subprocess.run(['squeue', '--job', job_id], stdout=subprocess.PIPE)
        output = result.stdout.decode('utf-8').strip()
        
        if job_id not in output:
            print(f"Job {job_id} completed.")
            break
        else:
            print(f"Job {job_id} is still running...")
            time.sleep(600)  # Check every 10 minutes

def main(project_dir):
    # Define the subdirectory where the .sh files are located
    sh_files_dir = os.path.join(project_dir, 'sh-files')

    # List of sbatch scripts to run sequentially
    sbatch_scripts = ['1-adapter-trimming.sh']

    # Verify project directory exists
    if not os.path.isdir(project_dir):
        print(f"Error: Project directory {project_dir} does not exist.")
        sys.exit(1)
        
    for script in sbatch_scripts:
        script_path = os.path.join(sh_files_dir, script)
        if not os.path.isfile(script_path):
            print(f"Error: Script {script_path} not found.")
            continue

        job_id = submit_job(script_path, project_dir)
        if job_id:
            wait_for_job(job_id)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 sequential_sbatch.py <project_directory>")
        sys.exit(1)
    
    project_directory = sys.argv[1]
    main(project_directory)

