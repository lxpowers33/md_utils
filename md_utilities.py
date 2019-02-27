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

def rmsd_pd_wrapper(selections, molid):
    nframes = molecule.numframes(molid)
    data = np.zeros((nframes,len(selections)))
    i = 0
    for item in selections: 
        data[:,i] = rmsd_from_initial(item['sel1'], molid)
        i = i + 1
    dataout = pd.DataFrame(data, columns = [p['name']+'_rmsd' for p in selections])
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

def projection_metric(selections, molid):
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
            data[t, i] = utils.v_projection(target.T - ref1.T, ref2.T - ref1.T) #project 1st onto 2nd
            i = i + 1
    dataout = pd.DataFrame(data, columns = [p['name'] for p in selections])
    return dataout

####################################################
#####            Analysis over Conditions      #####
####################################################

TM1TO4 = 'protein and name CA and (resid 4 to 37 or resid 39 to 69 or resid 75 to 110 or resid 119 to 146)'
def get_data_over_selections(info, align_sel=TM1TO4): 
    '''
    conditions_info should contain name, psf, paths, reps, selections
    '''
    #get the molid (just lodaded will be on top)
    mol_id = vmd.molecule.get_top()

    #align if desired
    if (align_sel != None):
        align_to_initial(align_sel, [mol_id])

    #perform all the anlyses using the selections in info 
    datasets = []
    for i in range(0, len(info['analyses'])):
	print('running analysis')
        analysis = info['analyses'][i]
	print(analysis)
        #perform analysis level alignment
        if 'alignments' in info:
            alignment_sel = info['alignments'][i]
            if alignment_sel != '':
                align_to_initial(alignment_sel, [mol_id])
        datasets.append(analysis(info['selections'][i], mol_id)) 

    #combine data from the different analyses
    dataout = pd.concat(datasets, axis=1)

    return dataout, False

def run_analysis_crystal(working_dir, save_dir, save_name, info, align_sel):
    out_dir = '{}/{}'.format(working_dir, save_dir)
    mol_id = molecule.load(info['file_type'], info['struct_file_name'])
    dataout, err = get_data_over_selections(info, align_sel)
    vmd.molecule.delete(mol_id)
    #Pickle and save
    if (not err):
        dataout.to_pickle(out_dir+'/{}.pkl'.format(save_name))

def run_analysis_traj(working_dir, save_dir, save_name, conditions, align_sel):
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
    f.write(json.dumps([str(c['analyses']) for c in conditions])+'\n')
    f.write(json.dumps(conditions[0]['selections'])+'\n')
    f.close()
    #collect data
    for condition in conditions: 
        condition_dir = '{}/{}/{}'.format(working_dir, save_dir, condition['name'])
        for rep in range(1,condition['reps']+1): 
            #Load the trajectory
            sname = join(condition['path'], 'prep', condition['psf'])
            tname = join(condition['path'], 'rep'+str(rep), condition['nc'])
            if not os.path.isfile(sname):
                print('file does not exist '+sname)
                return [], True
            if not os.path.isfile(tname):
                print('file does not exist '+tname)
                return [], True
            #load and align trajectories
            load_stride = int(5/condition['stride']) #for analysis every nanosecond
            #note 'stride' should be what you reimaged at 
            load_traj(sname, tname, stop = -1, stride = load_stride) 
            print('loaded trajectory successfully')
	    #get the data
            dataout, err = get_data_over_selections(condition, align_sel)
            #Delete the molecule
            mol_id = vmd.molecule.get_top()
	    vmd.molecule.delete(mol_id)
            #Make the save directory if it doesn't exist
            if (not os.path.exists(condition_dir)): 
                os.makedirs(condition_dir)
            #Pickle and save
            if (not err):
                dataout.to_pickle(condition_dir+'/{}_{}.pkl'.format(save_name, rep))

####################################################
#####            Development                   #####
####################################################


def calc_average_structure(molids, psf, minframe=0):
    """
        # From Robin - suggestion: should calculate per replicate average structure.
        # Then, compare each average structure to the per condition average to determine variation across replicates.
    Calculates the average structure for a given trajectory and psf

    Args:
        molids (list of int): Trajectory IDs to average
        psf (str): Path to psf file describing this topology
        minframe (int): Frame to start computation from
    """
    data = []
    start_frame = minframe
    total_frames = 0
    for m in molids:
        if start_frame >= molecule.numframes(m): continue
        summed = vmdnumpy.timestep(m, start_frame)
	total_frames = total_frames + 1
	for f in range(start_frame+1, molecule.numframes(m)):
            summed  = vmdnumpy.timestep(m, f) + summed
	    total_frames = total_frames + 1		
    avg = summed/total_frames

    # Now we have average coords, so set them in a new molecule
    if "_trans" in psf:
        pdb = psf.replace("_trans.psf", ".pdb")
    else:
        pdb = psf.replace(".psf", ".pdb")

    outid = molecule.load('psf', psf, 'pdb', pdb)
    atomsel("all", outid).set('x', avg[:,0])
    atomsel("all", outid).set('y', avg[:,1])
    atomsel("all", outid).set('z', avg[:,2])
    return outid

def calc_rmsf_to_average(molid, avg, selstr, minframe=0, calc_rmsd=False):
    """
        # From Robin
    Calculates the RMSF of an atom selection 

    Args:
        molid (int): VMD molecule ID to calculate
        avg (int): VMD molecule ID of average structure
        selstr (str): Selection to compute RMSF over
        minframe (int): Frame to start computation from
    """
    mask = vmdnumpy.atomselect(avg, 0, selstr)
    ref = np.compress(mask, vmdnumpy.timestep(avg,0), axis=0)

    if molecule.numframes(molid) <= minframe:
        print("Only %d frames in %d" % (molecule.numframes(molid), molid))
        return None

    N = len(ref)
    rmsf = np.zeros(N)
    rmsds = []
    for f in range(minframe, molecule.numframes(molid)):
        frame = np.compress(mask, vmdnumpy.timestep(molid, f), axis=0)
        rmsf += np.sum((frame-ref)**2, axis=1) #removed the extra squareroot
        if calc_rmsd:
            MSD = np.sum(np.sum((frame-ref)**2, axis=1))/N
            rmsds.append(np.sqrt(MSD))

    rmsf /= (molecule.numframes(molid)-minframe)
    rmsf = np.sqrt(rmsf)

    return rmsf, rmsds

def rmsd_average_wrapper(selections, molid):
    """
    Selections [{'sel1':'', 'name':''},{'sel1':'', 'name':''}]
    """
    avg_id = calc_average_structure([molid], selections[0]['psf'], minframe=0) 
    data = np.zeros((molecule.numframes(molid),len(selections)*2))

    i = 0
    for sel in selections:
        rmsf, rmsds = calc_rmsf_to_average(molid, avg_id, sel['sel1'], minframe=0, calc_rmsd=True)
        data[:,i] = rmsds
        i = i + 1
        data[0,i] = np.mean(rmsf)
        i = i + 1

    names = []
    for p in selections:
        names.append(p['name']+'_avg_rmsd')
        names.append(p['name']+'_avg_rmsf')
    dataout = pd.DataFrame(data, columns = names)
    vmd.molecule.delete(avg_id)
    return dataout

def fast_pair_distance(selections, molid):
    #make a dictionairy of base atoms 
    #generate a mask of the relevent atoms 
    #extract those positions for a particular frame 
        #frame = numpy.compress(mask, vmdnumpy.timestep(molid, f), axis=0)
        #add positions to array with time dimension 
    #for each pair use the relevent row to perform the operations
        #numpy.sqrt(numpy.sum((data[i,:,:] - data[j,:,:])**2, axis=1))
    """
    Selections [{'sel1':'', 'sel2':'', 'name':''},{'sel1':'', 'sel2':''}]
    """
    unique_atoms = set()
    #dump all the selections into a set
    for pair in selections:
        unique_atoms.add('( '+pair['sel1']+' )')
        unique_atoms.add('( '+pair['sel2']+' )')

    sel_to_index = {}
    for atom in unique_atoms:
        sel_to_index[atom] = np.nonzero(vmdnumpy.atomselect(molid, 0, atom))[0][0] 

    sel_all = ' or '.join(list(unique_atoms))
    mask = vmdnumpy.atomselect(molid, 0, sel_all)
    atom_index = np.nonzero(mask)[0]
    index_to_i = {x:i for i,x in enumerate(atom_index)}

    n_frames = molecule.numframes(molid)
    n_atoms = len(unique_atoms)
    positions = np.zeros((n_atoms, 3, n_frames))
    for f in range(0, n_frames):
        positions[:,:,f] = np.compress(mask, vmdnumpy.timestep(molid, f), axis=0)

    data = np.zeros((n_frames,len(selections)))
    s = 0
    for item in selections: 
        i = index_to_i[sel_to_index['( '+item['sel1']+' )']]
        j = index_to_i[sel_to_index['( '+item['sel2']+' )']]
        data[:, s] = np.sqrt(np.sum((positions[i,:,:] - positions[j,:,:])**2, axis=0))
        s = s + 1
    dataout = pd.DataFrame(data, columns = [p['name'] for p in selections])
    return dataout

