import numpy as np
import pandas as pd
from IPython import embed
from os import listdir
from os.path import isfile, join
from vmd import *
import utils
import datetime
import json

####################################################
#####             Load Trajectories            #####
####################################################
'''
Assumes the file structure
prep/
	5tyz_AP8_dabble.psf
rep1/
	traj_combined_stride1_reimaged.nc
rep2/
	traj_combined_stride1_reimaged.nc
...
etc. 

Example Usage
base = '/scratch/PI/rondror/MD_simulations/amber/GPR40/jpaggi/AP8/v1'
psf = '5tyz_AP8_dabble.psf'
nc = 'traj_combined_stride1_reimaged.nc'
load_replicates(base, psf, nc, replicates=range(1, 6), stride = 5, stop = 5500*5)
'''

def load_traj(struct_file_name, traj_file_name,
              start=0, stop=-1, stride=1, smoothing=0):
    filetype='netcdf'
    trajid = molecule.load('psf', struct_file_name)
    molecule.read(trajid, filetype, traj_file_name,
                  beg=start, end=stop, skip=stride, waitfor=-1)

def load_replicates(base, psf, nc,
                    replicates=range(1, 6), stride = 1, stop = -1):
    sname = join(base, 'prep', psf)
    print(sname)
    for rep in map(str, replicates):
        tname = join(base, 'rep'+rep, nc)
        print(tname)
        load_traj(sname, tname, stop = stop, stride = stride)

####################################################
#####             Alignment and RMSD           #####
####################################################
'''
Example usage: 
tm1to4 = 'protein and name CA and (resid 4 to 37 or resid 39 to 69 or resid 75 to 110 or resid 119 to 146)'
align_to_initial(tm1to4, range(5))
mk6_rmsds = [rmsd_from_initial('noh resname MK6', molid) for molid in range(5) if molecule.numframes(molid)]
'''

def align_to_initial(sel, molids):
    for molid in molids:
        ref = atomsel(sel, molid=molid, frame=0)
        for t in range(molecule.numframes(molid)):
            sel_t = atomsel(sel, molid=molid, frame=t)
            M = sel_t.fit(ref)
            atomsel('all',molid=molid, frame=t).move(M)

def rmsd_from_initial(sel, molid):
    sel_0 = atomsel(sel, molid=molid,frame=0)
    rmsds = np.zeros((molecule.numframes(molid),))
    for t in range(molecule.numframes(molid)):
        sel_t = atomsel(sel, molid=molid, frame=t)
        rmsds[t] = sel_0.rmsd(sel_t)
    return rmsds

def rmsd_pd_wrapper(sel, molid):
    data = rmsd_from_initial(sel, molid)
    dataout = pd.DataFrame(data, columns = ['rmsd'])
    return dataout

####################################################
#####            Distance analysis             #####
####################################################

def atom_xyz(atoms, mol_id):
    '''
    [{'name':'O1', 'sel':''}, {'name':'O', 'sel':''}]
    '''
    #each selection is a name and a sel
    labels = []
    n_atoms = len(atoms)
    end_frame = vmd.molecule.numframes(mol_id)
    data = np.zeros((end_frame, n_atoms*3))

    for i in range(0,n_atoms):
        labels = labels + [atoms[i]['name']+':'+d for d in ['x','y','z']]
        for f in range(0,end_frame): 
            sel = vmd.atomsel(atoms[i]['sel'], molid=mol_id, frame=f)
            data[f,3*i:(3*i+3)] = sel.get('x') + sel.get('y') + sel.get('z')

    dataout = pd.DataFrame(data, columns = labels)
    return dataout

def pair_dist_pd_wrapper(selections, molid): 
    """
    Selections [{'sel1':'', 'sel2':'', 'name':''},{'sel1':'', 'sel2':''}]
    """
    data = np.zeros((molecule.numframes(molid),len(selections)))
    for t in range(molecule.numframes(molid)):
        i = 0
        for item in selections: 
            sel1 = atomsel(item['sel1'], molid=molid, frame=t)
            sel2 = atomsel(item['sel2'], molid=molid, frame=t)
            data[t, i] = utils._calcdist(sel1, sel2)
            i = i + 1
    dataout = pd.DataFrame(data, columns = [p['name'] for p in selections])
    return dataout

def projection_metric(selections, mol_id):
    """
    Selections [{'ref1':'', 'ref2':'', 'target':'', 'name':''},{'sel1':'', 'sel2':''}]
    
    Measures the projection from target-ref1 onto ref2-ref1
    ref1------x ref2
      \   ^
       \  |
        \ |
         x
        target

    """
    data = np.zeros((molecule.numframes(molid),len(selections)))
    for t in range(molecule.numframes(molid)):
        i = 0
        for item in selections: 
            ref1 = utils._vectorize_coords(atomsel(item['ref1'], molid=molid, frame=t))
            ref2 = utils._vectorize_coords(atomsel(item['ref2'], molid=molid, frame=t))
            target = utils._vectorize_coords(atomsel(item['target'], molid=molid, frame=t))
            data[t, i] = utils.v_projection(target - ref1, ref2 - ref1) #project 1st onto 2nd
            i = i + 1
    dataout = pd.DataFrame(data, columns = [p['name'] for p in selections])
    return dataout


####################################################
#####            Analysis over Conditions      #####
####################################################

TM1TO4 = 'protein and name CA and (resid 4 to 37 or resid 39 to 69 or resid 75 to 110 or resid 119 to 146)'
def get_data_over_selections(info, rep, align_sel=TM1TO4): 
    '''
    conditions_info should contain name, psf, paths, reps, selections
    '''
    sname = join(info['path'], 'prep', info['psf'])
    tname = join(info['path'], 'rep'+str(rep), info['nc'])
    if not os.path.isfile(sname):
        print('file does not exist '+sname)
        return [], True
    if not os.path.isfile(tname):
        print('file does not exist '+tname)
        return [], True

    #load and align trajectories
    load_traj(sname, tname, stop = -1, stride = 5)

    #get the molid (just lodaded will be on top)
    mol_id = vmd.molecule.get_top()

    #align if desired
    if (align_sel != None):
        align_to_initial(align_sel, [mol_id])

    #perform all the anlyses using the selections in info 
    datasets = []
    for i in range(0, len(info['analyses'])):
        analysis = info['analyses'][i]
        datasets.append(analysis(info['selections'][i], mol_id)) 

    #combine data from the different analyses
    dataout = pd.concat(datasets, axis=1)

    #delete the molecule
    vmd.molecule.delete(mol_id)

    return dataout, False

def run_analyses(working_dir, save_dir, save_name, conditions, align_sel):
    '''
    Output 
        /savedir
            /conditions[0]['name']
                save_name_1.pkl
                save_name_2.pkl
                ...
            /conditions[0]['name']
            etc. 

    will write pkl files containing pandas dataframes
    Each column is a analyses and selection
    Each row is a nanosecond
    '''
    #write log file for analyses
    f = open(working_dir+'/'+save_name+".log", "a")
    f.write('\n \n \n'+'Log for data in '+save_name)
    f.write(datetime.datetime.now().strftime("%c")+'\n')
    f.write(json.dumps([c['name'] for c in conditions])+'\n')
    f.write(json.dumps(analyses)+'\n')
    f.write(json.dumps(conditions[0]['selections'])+'\n')
    f.close()
    #collect data
    for condition in conditions: 
        condition_dir = '{}/{}/{}'.format(working_dir, save_dir, condition['name'])
        for rep in range(1,condition['reps']+1): 
            dataout, err = get_data_over_selections(condition, rep, align_sel)
            if (not os.path.exists(condition_dir)): 
                os.makedirs(condition_dir)
            #pickle and save
            if (not err):
                dataout.to_pickle(condition_dir+'/{}_{}.pkl'.format(save_name, rep))

