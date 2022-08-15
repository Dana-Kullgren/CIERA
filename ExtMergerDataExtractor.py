
#-Looks at properties of BH Mergers 
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from time_conversion import read_units
from inspiral_time import inspiral_time_peters
import os,sys

# define paths
path0 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-1.6/1_1/'
path1 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.0/0_0/'
path2 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.3/0_0/'
path3 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.6/0_0/'
path4 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-3.0/0_0/'

paths = [path0, path1, path2, path3, path4]
alpha_values = ['1.6', '2.0', '2.3', '2.6', '3.0']
DataDir = './data/'

for idx in range(len(paths)):
    #--- Creates data files of all mergers involving BHs that merge outside the cluster 
    columns = ['t_merge', 'M1', 'M2','V_kick [km/s]','a[AU]', 'e', 'Type','Outcome', 'Spin1', 'Spin2']

    df = pd.DataFrame(columns=columns)

    outcomes = []
    times = []
    mass1 = []
    mass2 = []
    types = []
    v_kicks = []
    a_s = []
    e_s = []
    bin_flags = []
    spins1 = []
    spins2 = []

    path = paths[idx]
    print(f"Now collecting data for alpha value {alpha_values[idx]}")
    units=read_units(path+'initial')
    time_units = units[0]['t_myr']
    mass_units = units[0]['mstar_msun']

    # creating arrays with ejected info to compare:
    ejected_data = np.genfromtxt(path + 'initial.esc.dat')
    ejected_t = ejected_data[:,1]
    ejected_binflag = ejected_data[:,14]

    for l in range(len(ejected_binflag)):
        if (ejected_binflag[l] == 1) : #binary
            ejected_a = (ejected_data[l,19])
            ejected_e = (ejected_data[l,20])
            mass1s     = (ejected_data[l,15])
            mass2s     = (ejected_data[l,16])
            type1     = (ejected_data[l,22])
            type2     = (ejected_data[l,23])
            id1 = (ejected_data[l,17])

            if (mass1s >= 0.0 and type1 == 14) and (mass2s >= 0.0 and type2== 14): #BBH being ejected
                print(id1, mass1s,mass2s)
                inspiral_time = inspiral_time_peters(ejected_a,ejected_e,mass1s,mass2s) #calculating inspiral time in Gyr
                tot_time = ejected_t[l]*time_units + inspiral_time*10**3 #in myr
                if tot_time <= 13700.0:
                    # mass1 = np.append(mass1, ejected_data[l,15])
                    # mass2 = np.append(mass2, ejected_data[l,16])
                    # times = np.append(times, tot_time)
                    # types = np.append(types, 'Ejected Binary')
                    # outcomes = np.append(outcomes, 'N/A')
                    # a_s = np.append(a_s, ejected_a)
                    # e_s = np.append(e_s, ejected_e)
                    # spins1 = np.append(spins1, ejected_data[l,58])
                    # spins2 = np.append(spins2, ejected_data[l,59])
                    # v_kicks = np.append(v_kicks, 'N/A')

                    mass1.append(ejected_data[l,15])
                    mass2.append(ejected_data[l,16])
                    times.append(tot_time)
                    types.append('Ejected Binary')
                    outcomes.append('N/A')
                    a_s.append(ejected_a)
                    e_s.append(ejected_e)
                    spins1.append(ejected_data[l,58])
                    spins2.append(ejected_data[l,59])
                    v_kicks.append('N/A')
                    
    # print(len([model]*len(ejected_times)), len(ejected_times), len(ejected_m1),len(ejected_m2))
    df2 = pd.DataFrame({'t_merge': times, 'M1': mass1, 'M2': mass2,'V_kick [km/s]': v_kicks,'a[AU]': a_s, 'e': e_s, 'Type': types, 'Outcome':outcomes, 'Spin1':spins1, 'Spin2': spins2})
    df = df.append(df2,ignore_index = True)[df2.columns.tolist()]

    df[columns].to_csv(f"{DataDir}bbh_ejected_{alpha_values[idx]}.csv", index=False)
