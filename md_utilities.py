import numpy as np
import pandas as pd
from IPython import embed
from os import listdir
from os.path import isfile, join
from vmd import *
from utils import *

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
#####            Analysis over Conditions      #####
####################################################

TM1TO4 = 'protein and name CA and (resid 4 to 37 or resid 39 to 69 or resid 75 to 110 or resid 119 to 146)'
def get_data_over_selections(info, rep, analysis, dim, align_sel=TM1TO4): 
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
    #check that the mol id is 0
    mol_id = vmd.molecule.get_top()
    print('mold ID is '+str(mol_id))

    align_to_initial(align_sel, [mol_id])
    
    dataout = analysis(info['selections'], mol_id)
    vmd.molecule.delete(mol_id)

    return dataout, False

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
