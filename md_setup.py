import os 
from datetime import datetime

Condition_purpose = "Describe conditions and file system"
Version_purpose = "Document versions"
Prep_purpose = "Document steps for preparing system"
version_starter_text = "v1 is my first version of simulations. Any alterations made in later versions will be described here."
prep_starter_text = "<Description of all steps goes here>"

def test():
	print("test worked!! :)")

def new_folder(path):
	#check that folder does not already exist
	if not os.path.exists(path):
		os.makedirs(path)

def new_readme(path,maintext,purpose,author,project): 
	f = open(path+"/README.txt","w+")
	f.write("PURPOSE: %s \n" % (purpose))
	f.write("PROJECT: %s \n" % (project))
	f.write("AUTHOR: %s \n" % (author))
	f.write("PATH: %s \n" % (path))
	f.write("DATE: %s \n \n" % (datetime.now().strftime("%Y-%m-%d")))
	f.write(maintext)
	f.close() 

def build_new_condition_directory(directory, c, config):
	#build condition directory
	directory_condition = directory+'/'+c["name"]
	new_folder(directory_condition)
	#build version directory
	directory_version = directory_condition +'/v1'
	new_folder(directory_version)
	#make the readme file for versioning
	new_readme(directory_condition, version_starter_text, Version_purpose, config['creator'], config['project'])
	#build prep folder for condition version
	directory_prep = directory_version+'/prep'
	new_folder(directory_prep)
	#make the readme file for prep 
	new_readme(directory_prep, prep_starter_text, Prep_purpose, config['creator'], config['project'])
	#build rep directories
	for i in range(1,c["reps"]+1):
		new_folder(directory_version+'/rep'+str(i))

def build_new_simulation_directory(directory, config):
	#check that the initial path exists
	if not os.path.exists(directory):
		raise Exception('The directory you provided does not exist \n %s' % (directory))
	for c in config["conditions"]:
		build_new_condition_directory(directory, c, config)
	#build the readme file for conditions and versioning methods
	new_readme(directory, config["conditions_text"], Condition_purpose, config['creator'], config['project'])

from os import listdir

def setup_simulation_files(directory, condition, version, reps, input_dir):
	'''
	Usage: setup_simulation_files(directory, 'MK6_loop', 'v1', 5, 'input_dir')
	'''
	f = open("setup_condition.sh","w+")
	for c in range(1,reps+1):

		target_dir = directory + '/' + condition + '/' + version +'/rep' + str(c)
		
		#copy all the input files 
		input_files = [f for f in listdir(input_dir) if f.endswith(".in")]
		for filename in input_files: 
			f.write("cp %s %s \n" % (input_dir+'/'+filename, target_dir+'/'+filename))
		
		#change the job name on the 4th line of sbatch
		jobname = 'GRP40_'+condition+'_'+version+'_'+reps
		f.write('sed s/# --job-name=gpr40_binary/--job-name=%s/ %s/sim_mdstep.sbatch > %s/sim_mdstep.sbatch \n' % (jobname, input_dir, target_dir))
		f.write('chmod +x %s/sim_mdstep.sbatch \n' % (target_dir)) #give it executible persmission
		
		#do the symlinks for the system files from the prep folder
		prep_dir = target_dir = directory + '/' + condition + '/' + version +'/prep'
		f.write('ln -s %s/system.prmtop %s/system.prmtop' %(prep_dir, target_dir))
		f.write('ln -s %s/system.inpcrd %s/system.inpcrd' %(prep_dir, target_dir))
		f.write('ln -s %s/system.psf %s/system.psf' %(prep_dir, target_dir))
	f.close() 
