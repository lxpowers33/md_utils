import os 
from datetime import datetime
import sys

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

def setup_simulation_files(directory, condition, version, reps, input_dir, project):
    '''
    Usage: setup_simulation_files(directory, 'MK6_loop', 'v1', 5, 'input_dir')
    '''
    f = open("setup_condition.sh","w+")
    for c in range(1,reps+1):

        target_dir = directory + '/' + condition + '/' + version +'/rep' + str(c)

        #copy all the input files
        input_files = [n for n in listdir(input_dir) if n.endswith(".in")]
        for filename in input_files:
            f.write("cp %s %s \n" % (input_dir+'/'+filename, target_dir+'/'+filename))

        #change the job name on the 4th line of sbatch
        jobname = project+'_'+condition+'_'+version+'_'+str(c)
        f.write("sed 's/#SBATCH --job-name=gpr40_binary/#SBATCH --job-name=%s/' %s/sim_mdstep.sbatch > %s/sim_mdstep.sbatch \n" % (jobname, input_dir, target_dir))
        f.write('chmod +x %s/sim_mdstep.sbatch \n' % (target_dir)) #give it executible persmission

        #do the symlinks for the system files from the prep folder
        prep_dir  = directory + '/' + condition + '/' + version +'/prep'
        for ext in ['prmtop', 'inpcrd', 'psf']:
            f.write('cp {}/{}_dabbled.{} {}/system.{} \n'.format(prep_dir, condition, ext, prep_dir, ext))

        f.write('ln -s %s/system.prmtop %s/system.prmtop \n' %(prep_dir, target_dir))
        f.write('ln -s %s/system.inpcrd %s/system.inpcrd \n' %(prep_dir, target_dir))
        f.write('ln -s %s/system.psf %s/system.psf \n' %(prep_dir, target_dir))
    f.close()

def add_simulation_reps(directory, condition, version, reps, input_dir, project):
    '''
    Usage: setup_simulation_files(directory, 'MK6_loop', 'v1', 5, 'input_dir')
    '''
    f = open("setup_condition.sh","w+")
    for c in reps:
        target_dir = directory + '/' + condition + '/' + version +'/rep' + str(c)
        #make the replicate files
        new_folder(target_dir)

        #copy all the input files
        input_files = [n for n in listdir(input_dir) if n.endswith(".in")]
        for filename in input_files:
            f.write("cp %s %s \n" % (input_dir+'/'+filename, target_dir+'/'+filename))

        #change the job name on the 4th line of sbatch
        jobname = project+'_'+condition+'_'+version+'_'+str(c)
        f.write("sed 's/#SBATCH --job-name=gpr40_binary/#SBATCH --job-name=%s/' %s/sim_mdstep.sbatch > %s/sim_mdstep.sbatch \n" % (jobname, input_dir, target_dir))
        f.write('chmod +x %s/sim_mdstep.sbatch \n' % (target_dir)) #give it executible persmission

        #do the symlinks for the system files from the prep folder
        prep_dir  = directory + '/' + condition + '/' + version +'/prep'
        f.write('ln -s %s/system.prmtop %s/system.prmtop \n' %(prep_dir, target_dir))
        f.write('ln -s %s/system.inpcrd %s/system.inpcrd \n' %(prep_dir, target_dir))
        f.write('ln -s %s/system.psf %s/system.psf \n' %(prep_dir, target_dir))
    f.close()

def write_reimage_sh(commands, name):
    with open(name, 'w') as f:
        f.write('#!/bin/bash\n')
        for command in commands:
            f.write(command + '\n')


def reimage(base_directory, c_name, version, reps, stride):
    jobs = []
    for r in range(1, reps+1):
        version_folder = base_directory + '/' + c_name + '/v' + str(version) + '/'
        rep_folder = version_folder + 'rep{}'.format(r)
        coord_file = rep_folder + '/system.inpcrd'
        parm_file = rep_folder  + '/system.psf'
        target = rep_folder
        cmd = 'reimage.py {} {} {} --stride {}  --include-eq {}'.format(parm_file, coord_file, rep_folder, stride, target)
        print(cmd)
        jobs.append(cmd)

    write_reimage_sh(jobs, 'reimage.sh')
    print('sbatch -t 00:20:00 -p rondror reimage.sh')
    os.system('sbatch -t 00:20:00 -p rondror reimage.sh')

if __name__ == '__main__':
    #first task is to build the new condition folder
    task = sys.argv[1]
    if task == 'help':
        print("         newc to create a new folder \n\
        arguments are: newc PENT_MK6 3 GPR40 (folder name, reps, project) \n\n\
        dabble to run dabble with default arguments (in prep folder) \n\
        arguments are: dabble PENT_MK6 \n\n\
        sim_files to build the simulation files in rep directories \n\
        arguments are: sim_files PENT_MK6 v1 GPR40 membrane (or water)\n\n\
        run_add to add to list of simulations to keep running \n\
        arguments are: run_add PENT_MK6 1 GPR40\n\n\
        reimage to reimage condition \n\
        arguments are: reimage PENT_MK6 1 5")
    elif task == 'newc':
        #arguments are  newc PENT_MK6 3 GPR40
        gen_config = {"creator": "Alex Powers", "project": sys.argv[4]}
        cond_config = {"name": sys.argv[2], "reps": int(sys.argv[3])}
        directory = os.getcwd()
        build_new_condition_directory(directory, cond_config, gen_config)
    elif task == 'dabble':
        #collect name of to_dabble.mae file
        file_prefix = sys.argv[2]
        files = os.listdir(os.getcwd())
        strfiles = [f for f in files if f.endswith('.str')]
        strs = ''.join(["-str {} ".format(f) for f in strfiles])
        cmd = ''
        if sys.argv[3] == '1':
            print('option 1 normal size is 80')
            cmd = "~/miniconda3/envs/conda2.7/bin/dabble -ff charmm36m --absolute-x 90.0 --absolute-y 90.0 -w 10 --hmr \
                    -i {}_to_dabble.mae {} -o {}_dabbled.prmtop -O | tee dabble.log".format(file_prefix, strs,
                                                                                            file_prefix)
        if sys.argv[3] == '2':
            print('option 2')
            cmd = "~/miniconda3/envs/conda2.7/bin/dabble -ff charmm36m --absolute-x 122 --absolute-y 122 -w 10 --hmr \
                                -i {}_to_dabble.mae {} -o {}_dabbled.prmtop -O | tee dabble.log".format(file_prefix,
                                                                                                        strs, file_prefix)
        
	if sys.argv[3] == '3':
            print('option 3')
            cmd = "~/miniconda3/envs/conda2.7/bin/dabble -ff charmm36m -M TIP3 -w 14 --hmr \
                    -i {}_to_dabble.mae {} -o {}_dabbled.prmtop -O | tee dabble.log".format(file_prefix, strs,
                                                                                            file_prefix)
	print(cmd)
        os.system(cmd)
    elif task == 'sim_files':
        directory = os.getcwd()
        condition = sys.argv[2]
        version = sys.argv[3]
        project = sys.argv[4]
        target_dir = directory + '/' + condition + '/' + version
        reps = len([f for f in os.listdir(target_dir) if f.startswith("rep")]) #get number of prep directories
        if sys.argv[5] == 'membrane':
        	input_dir = '/home/users/lxpowers/projects/GPR40/input_files'
        if sys.argv[5] == 'water':
        	input_dir = '/home/users/lxpowers/projects/CREB/input_files'
        setup_simulation_files(directory, condition, version, reps, input_dir, project)
        os.system('chmod +x setup_condition.sh')
    elif task == 'run_add':
        directory = os.getcwd()
        c = sys.argv[2]
        v = str(int(sys.argv[3]))
        p = sys.argv[4]

        target_dir = directory + '/' + c + '/v' + v
        nreps = len([f for f in os.listdir(target_dir) if f.startswith("rep")])
        r_string = ','.join([str(i) for i in range(1,nreps+1)])

        add_file = '/home/users/lxpowers/analyze_running/run_config.py'
        with open(add_file) as f:
            content = f.readlines()
        print(content)
        content.append(']\n')
        newc =  "{{'c':['{}'], 'v':[{}], 'r':[{}], 'p':'{}', 'd': '{}'}},\n".format(c, v, r_string, p, directory)
        content[-2] = newc
        with open(add_file, "w") as f:
            for line in content:
                # write line to output file
                f.write(line)
    elif task == 'reimage':
        directory = os.getcwd()
        c = sys.argv[2]
        v = str(int(sys.argv[3]))
        s = sys.argv[4]
        target_dir = directory + '/' + c + '/v' + v
        nreps = len([f for f in os.listdir(target_dir) if f.startswith("rep")])
        reimage(directory, c, v, nreps, s)
    else:
        print('invalid command')


