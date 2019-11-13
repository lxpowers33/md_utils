import os
import re
import subprocess
import time

REP_FOLDER_RE = re.compile("rep[0-9]+")

#Styles of submissions
#submitting a new job
#to_resubmit(base_directory, version, conditions, time, dryrun=True, reps=[1,2,3,4,5])
#updating a select set
#to_resubmit(base_directory, version, conditions, time, update='-u', dryrun=False, reps=[1,2,3,4,5])

def main():
	base_directory = '/scratch/PI/rondror/MD_simulations/amber/muscarinic/M1active/M1_Karuna'
	version = 1
	conditions = ['M1_xan_pose6'] #MK6_loop
	time = 500
	to_resubmit(base_directory, version, conditions, time, dryrun=True)


def cancel_pending_jobs_for_conditions(conditions):
	cmd = 'squeue -u $USER -t PENDING --format "%.18i %.9P %.j %.8u %.2t %.10M %.6D %R"'
	pipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout
	jobs = pipe.read()[:-1].split('\n')
	for job_line in jobs:
		job_info = parse_slurm_job_line(job_line)
		print(job_info)
		for c in conditions:
			if (c == job_info['name'][6:-5] or c == job_info['name'][6:-6]):
				print('stopping {} {}'.format(job_info['name'],job_info['id']))
				os.system('scancel {}'.format(job_info['id']))
	return jobs

def parse_slurm_job_line(job_line):
	l = job_line.split(' ')
	l = [i for i in l if i != '']
	job_info = {}
	job_info['id'] = l[0]
	job_info['name'] = l[2]
	return job_info

resubmit_folder = '/home/users/lxpowers/general_code/resubmit'

def to_resubmit(base_directory, version, conditions, time_run, update = '', dryrun=True, reps=[], immediate=False):
	for c_name in conditions:
		version_folder = base_directory + '/' + c_name + '/v' + str(version) + '/'
		if reps == []: 
			rep_folders = get_replicate_folders(version_folder)
		else:
			rep_folders = [version_folder+'rep'+str(r) for r in reps]
			for rep_folder in rep_folders:
				print(rep_folder)	
				cmd = 'python {}/resubmit.py  --userdb {}  {}/mdinfo {}/sim_mdstep.sbatch {}'.format(resubmit_folder,
					update, rep_folder, rep_folder, time_run)																		  
				print('sbatch {}/sim_mdstep.sbatch'.format(rep_folder))
				print(cmd)
				if not dryrun:
					os.system(cmd)
					if (time_run != 0 and immediate): 
						os.chdir(rep_folder)
						os.system('sbatch sim_mdstep.sbatch')
						os.system('sbatch {}/sim_mdstep.sbatch'.format(rep_folder))
					else:
						print('not running')
						time.sleep(1)

def get_replicate_folders(folder):
	'''
	Return a list of all replicate folders in the version folder
	:param folder: list corresponding to version directory
	:return: list of strings
	'''
	rep_folders = []
	listing = os.listdir(folder)
	for name in listing:
		if REP_FOLDER_RE.match(name) is not None:
			rep_folders.append(name)
	rep_folders = [os.path.join(folder, name) for name in rep_folders]
	return sorted(rep_folders)

if __name__== "__main__":
	main()