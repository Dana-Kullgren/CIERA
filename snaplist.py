import numpy as np, pandas as pd
from h5py import File as h5pyFile
default_sim_path = '/projects/b1095/newlin/cmc/IMF_fbh_grid/rundir/rv1/rg8/z0.002/n8e5/w5fb0.05fbh0.05alpha3-2.3/0_0'

##### Utility/Support Functions #####
def conv(unit,filepath):
    ''' Returns the unit conversion multiplier given a simulation's *.conv.sh file and a unit from the below dictionary. '''
    dict = {'m_cgs':5,'m':7,'mstar_cgs':9,'mstar':11,'l_cgs':13,'l':15,'t_cgs':17,'t':19,'tnbody_cgs':21,'tnbody':23}
    from re import findall
    with open(filepath,'r') as f: head = [next(f) for x in range(24)]
    findconv = findall('\d+[\.]?\d*e\+\d*',head[dict[unit]])
    if   len(findall('\d+[\.]?\d*e\+\d*',head[dict[unit]])) == 1: return(expnot(findconv[0])) # If conversion needs to be transformed from exponential to standard notation
    elif len(findall('\d+[\.]?\d*e\+\d*',head[dict[unit]])) == 0: return(float(findall('\d+[\.]?\d*',head[dict[unit]])[0])) # Conversion already in standard notation
def conv_all(filepath): return((conv('m',filepath), conv('mstar',filepath), conv('t',filepath), conv('tnbody',filepath), conv('l',filepath)))
def conv_all_cgs(filepath): return((conv('m_cgs',filepath), conv('mstar_cgs',filepath), conv('t_cgs',filepath), conv('tnbody_cgs',filepath), conv('l_cgs',filepath)))

'''
Usage of conv:
Say you want to plot the positions of stars (array-like 'r') from a snapshot of a CMC simulation (specified by its filepath/directory 'path')
To convert r from CMC code units to parsecs, multiply it by l_conv defined below:
  (m_conv, mstar_conv, t_conv, tnbody_conv, l_conv) = conv_all(path)
To convert r from CMC cod units to cm, multiply it by l_conv defined below:
  (m_conv, mstar_conv, t_conv, tnbody_conv, l_conv) = conv_all_cgs(path)
Note that many columns in CMC output already have a unit specified in their column header. For those that don't, the general rule is as follows:
  - Masses: multiply lengths by m_conv, except for very rare cases (e.g., collision files), when mstar_conv may be neaded
  - Lengths: multiply by l_conv
  - Times: multiply by t_conv
  - Velocities and Specific Angular Momentum J: multiply by (l_conv/tnbody_conv)
  - Specific Energies: multiply by (l_conv/tnbody_conv)**2
'''


'''
Fastest way to load CMC .dat files (e.g. the initial.dyn.dat file):
  data = pd.read_csv(path+'/initial.dyn.dat',skiprows=1,header=None,delim_whitespace=True,dtype='str',usecols=(0,1,2)).values.T
  data[np.where(esc_data == 'na')] = 'nan'
  data = data.astype(float)
Make sure to set skiprows correctly to the number of rows at the start of the file that are commented out
and set usecols to load only the columns that you need (for fastest loading). To load just one column, set usecols='(<colnum>,)'
'''

def list_snapshots(path=default_sim_path):
    ''' Return the snapshot keys and ages for the specified CMC simulation '''
    f = h5pyFile(path+'/initial.window.snapshots.h5', 'r')
    keys = f.keys()
    snapno_of_window = [key.split('(')[0] for key in keys]
    age = [key.split('=')[1].split('G')[0] for key in keys]
    snapshot_keys_ordered_by_age = np.array(['%s(t=%sGyr)'%(s,a) for (a,s) in sorted(zip(age,snapno_of_window))],dtype='str')
    snapshot_ages_Myr = np.array(list(map(float,sorted(age)))) * 1e3
    f.close()
    return (snapshot_keys_ordered_by_age, snapshot_ages_Myr)

def load_snapshot(path=default_sim_path,key=None,return_dataframe=0,columns=np.arange(0,62)):
    ''' Load a CMC snapshot
    Input:
      - path (str, optional): filepath to the directory containing the CMC simulation, defaults to default_sim_path
      - key (str, optional): the key corresponding to the desired simulation snapshot, defaults to last snapshot
      - return_dataframe (bool, optional): if 1, returns snapshot as pandas dataframe; if 0 (default), returns snapshot as a numpy array
      - columns (array-like, optional): if return_dataframe=0, returns only the desired columns from the simulation snapshot, defaults to all columns (which takes more time to load)
                                        the columns can be specified as ints, e.g., for the first three columns [0,1,2]
    '''
    if key == None:
        key = list_snapshots(path)[0][-1]
        #snapshot_age = 1e3 * float(key.split('=')[1].split('G')[0]) # Units: [Myr]
    if return_dataframe == 1: # return raw dataframe (can be more annoying to deal with than it's worth)
        return pd.read_hdf(path+'/initial.window.snapshots.h5',key=key)
    else: # return a simple numpy array of the desired snapshot columns
        return pd.read_hdf(path+'/initial.window.snapshots.h5',key=key).values.T[np.array(columns).astype(int)]
