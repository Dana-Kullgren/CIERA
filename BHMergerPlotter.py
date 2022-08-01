# import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from re import *

##### Utility/Support Functions #####
def expnot(string): # Converts from shorthand exponential notation (‘4e3’, ‘1e+2’ or ‘3.21e-1’) to a normal float (4000., 100., 0.321)
    if   len(findall('e\+',string)) == 1: n = float(string.split('e+')[0])*10**float(string.split('e+')[1])
    elif len(findall('e\-',string)) == 1: n = float(string.split('e-')[0])*10**(-1.*float(string.split('e-')[1]))
    elif len(findall('e',string))   == 1: n = float(string.split('e')[0])*10**float(string.split('e')[1])
    elif len(findall('e',string))   == 0: n = float(string)
    else: n = "improper input"
    return n

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

def change_to_floats(col):
    '''Remove 'na' values and change types in array to floats
    Input:
      - col (array-like): the array that will have its entries changed to floats
    '''
    col[np.where(col == 'na')] = 'nan'
    col = col.astype(np.float)
    col = np.array(col)
    return col

# define paths and directories
path0 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-1.6/1_1'
path1 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.0/0_0'
path2 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.3/0_0'
path3 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.6/0_0'
path4 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-3.0/0_0' #THIS IS STILL RUNNING

paths = [path0, path1, path2, path3, path4]

# load and format data
IncludeLIGO = True

alpha1_6 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0, 'Z': 0}
alpha2_0 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0, 'Z': 0}
alpha2_3 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0, 'Z': 0}
alpha2_6 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0, 'Z': 0}
alpha3_0 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0, 'Z': 0}
alpha_data = [alpha1_6, alpha2_0, alpha2_3, alpha2_6, alpha3_0]
alpha_values = ['1.6', '2.0', '2.3', '2.6', '3.0']
alphas_used = []
DataDir = './data/'
PlotDir = './plots'

print('Loading data...')

for i in range(len(paths)):
	if paths[i] != '':
		
		# download initial simulation data
		M1, M2 = pd.read_csv(paths[i]+'/initial.bhmerger.dat',skiprows=1,header=0,delim_whitespace=True,dtype='str',usecols=(5,6)).values.T
		r_c, r_h, time = pd.read_csv(paths[i]+'/initial.dyn.dat',skiprows=1,header=0,delim_whitespace=True,dtype='str',usecols=(7,20,0)).values.T
		
		# get metallicities
		Z_str = paths[i][50:56]
		print(f'Z_str = {Z_str}')
		Z_float = float(Z_str)

		# add in ejected BHs (masses are in M_SUN)
		ejected_data = pd.read_csv(f"{DataDir}bbh_ejected_{alpha_values[i]}.csv")
		ejected_M1 = ejected_data['M1']
		ejected_M2 = ejected_data['M2']
		ejected_M1 = np.array(ejected_M1)
		ejected_M2 = np.array(ejected_M2)
		# print('ejected_M1 =', ejected_M1)
		# print('ejected_M2 =', ejected_M2)
		# print('type(ejected_M1) =', type(ejected_M1))
		# print('type(M1) =', type(M1))
		M1 = np.append(M1, ejected_M1)
		M2 = np.append(M2, ejected_M2)

		# change data to floats
		M1 = change_to_floats(M1)
		M2 = change_to_floats(M2)
		ejected_M1 = change_to_floats(ejected_M1)
		ejected_M2 = change_to_floats(ejected_M2)
		r_c = change_to_floats(r_c)
		r_h = change_to_floats(r_h)
		time = change_to_floats(time)

		(m_conv, mstar_conv, t_conv, tnbody_conv, l_conv) = conv_all(paths[i]+'/initial.conv.sh')
		# M1 *= m_conv    ## this is already in M_sun
		# M2 *= m_conv    ## this is already in M_sun
		r_c *= l_conv
		r_h *= l_conv
		time *= t_conv

		for j in range(len(M1)): # Ensures that M1 is the larger body
			# if 75 < M1[j] and M1[j] < 85:
			# 	print(f'M1 = {M1[j]}')
			if M2[j] > M1[j]:
				larger = M2[j]
				M2[j] = M1[j]
				M1[j] = larger
		for k in range(len(ejected_M1)):
			# if 75 < M1[k] and M1[k] < 85:
			# 	print(f'M1 = {M1[k]}')
			if ejected_M2[k] > ejected_M1[k]:
				larger = ejected_M2[k]
				ejected_M2[k] = ejected_M1[k]
				ejected_M1[k] = larger

		Z_list = [Z_float for i in range(len(M1))]
		Z_list = np.array(Z_list)
		ejected_Z_list = [Z_float for i in range(len(ejected_M1))]
		ejected_Z_list = np.array(ejected_Z_list)

		print(f'alpha = {alpha_values[i]}')
		print(f'shape(M1) = {np.shape(M1)}')
		print(f'shape(M2) = {np.shape(M2)}')

		print(f'shape(Z_list) = {np.shape(Z_list)}')
		print(f'shape(ejected_Z_list) = {np.shape(ejected_Z_list)}')

		alpha_data[i]['M1'] = M1 						## M1, M2, and Z contain the values for all BHs, whether or not they were ejected
		alpha_data[i]['M2'] = M2
		alpha_data[i]['ejected_M1'] = ejected_M1 		## The ejected keys contain ONLY the ejected cluster values
		alpha_data[i]['ejected_M2'] = ejected_M2
		alpha_data[i]['r_c'] = r_c
		alpha_data[i]['r_h'] = r_h
		alpha_data[i]['time'] = time
		alpha_data[i]['Z'] = Z_list
		alpha_data[i]['ejected_Z'] = ejected_Z_list

		# for 
		# 	print()

		# print(f"type(alpha_data[{i}]['M2']) = {type(alpha_data[i]['M2'])}")
		# for l in range(len(alpha_data[i]['M2'])):
		# 	m2 = alpha_data[i]['M2'][l]
		# 	if m2 > 80:
		# 		print(f'M2 = {m2}')

		# print(f"For alpha={alpha_values[i]}, the average r_c is {np.mean(r_c)}")
		# print(f"For alpha={alpha_values[i]}, the average r_h is {np.mean(r_h)}")
		# print('len(r_c) =', len(r_c))
		# print('len(r_h) =', len(r_h))

		# print(alpha_data[i])

		alphas_used.append(alpha_values[i])

print('Simulation data loaded...')

if IncludeLIGO:
	LIGO_data = pd.read_csv(DataDir + 'LIGO_data.csv')
	# print(LIGO_data.columns)
	# BBH_idx = np.where(M2 > 3)
	# print("BBH_idx:", BBH_idx)

	M1 = LIGO_data['mass_1_source']
	M1_low = -1.0 * LIGO_data['mass_1_source_lower']
	M1_up = LIGO_data['mass_1_source_upper']
	M2_low = -1.0 * LIGO_data['mass_2_source_lower']
	M2_up = LIGO_data['mass_2_source_upper']
	M2_unfiltered = LIGO_data['mass_2_source']
	
	# Remove events that aren't BBH mergers
	M1 = M1[M2_unfiltered > 3]
	M1_low = M1_low[M2_unfiltered > 3]
	M1_up = M1_up[M2_unfiltered > 3]
	M2_low = M2_low[M2_unfiltered > 3]
	M2_up = M2_up[M2_unfiltered > 3]
	M2 = M2_unfiltered[M2_unfiltered > 3]

	# print("len(M1) =", len(M1))
	# print("len(M1_low) =", len(M1_low))
	# print("len(M1_up) =", len(M1_up))
	# print("len(M2) =", len(M2))
	# print("len(M2_low) =", len(M2_low))
	# print("len(M2_up) =", len(M2_up))

	LIGO_data = {'M1': M1,
			'M1_error': [M1_low, M1_up],
			'M2': M2,
			'M2_error': [M2_low, M2_up]}

	print('LIGO data loaded...')

print('All data loaded...')
print('Beginning plotting...')


# plot data
######################################################################### Mass Histograms

f, ax = plt.subplots(1,2,figsize=(14,7))
w = 10
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i in range(len(alpha_data)):
	masses = [alpha_data[i]['M1'], alpha_data[i]['M2']]
	LIGO_masses = [LIGO_data['M1'], LIGO_data['M2']]
	LIGO_err = [LIGO_data['M1_error'], LIGO_data['M2_error']]
	for j in range(len(masses)):
		if (i==0) and (j==0) and (IncludeLIGO): # only do this once
			for LIGO_idx in range(len(LIGO_masses)):
				y0 = LIGO_masses[LIGO_idx]
				y0_up = y0 - LIGO_err[LIGO_idx][0]
				y0_low = y0 + LIGO_err[LIGO_idx][1]
				y00=np.sort(y0)
				y00_up=np.sort(y0_up)
				y00_low=np.sort(y0_low)
				p00=1.*np.arange(len(y0))/(len(y0)-1)
				ax[LIGO_idx].fill_betweenx(p00, y00_low, y00_up, alpha=0.3, color='black', label='LIGO error')
				# ax[LIGO_idx].plot(y00_up, p00, linestyle=':', drawstyle='steps', linewidth=3.5, c='grey', label='LIGO data error')
				# ax[LIGO_idx].plot(y00_low, p00, linestyle=':', drawstyle='steps', linewidth=3.5, c='grey', label='LIGO data error')
				ax[LIGO_idx].plot(y00, p00, linestyle='-', drawstyle='steps', linewidth=3.5, c='black', label='LIGO data')
		y0 = masses[j]
		# weights = np.ones_like(y0)/float(len(y0))
		y00=np.sort(y0)
		p00=1.*np.arange(len(y0))/(len(y0)-1)
		plot00=ax[j].plot(y00, p00, linestyle='-', drawstyle='steps', linewidth=3.5, label=f'alpha={alpha_values[i]}')

	# fill
	# ax[0].hist(alpha_data[i]['M1'], fill=True, lw=1, alpha=1, cumulative=False, histtype='step', bins=np.arange(0, 90+w, w), label=f'alpha={alpha_values[i]}')
	# ax[1].hist(alpha_data[i]['M2'], fill=True, lw=1, alpha=1, cumulative=False, histtype='step', bins=np.arange(0, 90+w, w))	

	# lines
	# ax[0].hist(alpha_data[i]['M1'], fill=False, color=colors[i], lw=2, histtype='step', bins=np.arange(0, 90, w))
	# ax[1].hist(alpha_data[i]['M2'], fill=False, color=colors[i], lw=2, histtype='step', bins=np.arange(0, 90, w))

plot_label = ''
for label in alphas_used:
	if label != alphas_used[0]:
		plot_label += ','
	plot_label += label

# print(f'Alpha values plotted: {plot_label}')

ax[0].set_xlabel(r'Mass of Primary ($M_{SUN}$)', fontsize='x-large')
ax[1].set_xlabel(r'Mass of Secondary ($M_{SUN}$)', fontsize='x-large')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax[1].legend(by_label.values(), by_label.keys(), loc='lower right')

f.suptitle('Primary and Secondary Masses of Black Hole Mergers', fontsize='xx-large')
f.savefig(PlotDir + f'/MassDistributions_WithEjected_alpha={plot_label}.pdf', bbox_inches='tight')
print('Mass Histograms complete')

######################################################################### M2 vs M1

f, ax = plt.subplots(2,3,figsize=(3*len(alphas_used),14))
ax[1][2].set_visible(False)

# ax_idx = 0

# calculate axes limits
x_max = 0
y_max = 0
for i in range(len(paths)):
	xi_max = np.amax(alpha_data[i]['M1'])
	yi_max = np.amax(alpha_data[i]['M2'])
	if xi_max > x_max:
		x_max = xi_max
	if yi_max > y_max:
		y_max = yi_max
print(f'x_max = {x_max}, {type(x_max)}')
print(f'y_max = {y_max}, {type(y_max)}')
x_lim = x_max + 10
y_lim = y_max + 30
x = np.linspace(0,x_lim,100)

for i in range(len(paths)):
	if paths[i] != '':
		# print('i =', i)
		if i > 2:
			row = 1
		else:
			row = 0
		col = i % 3
		if IncludeLIGO:
			# ax[ax_idx].errorbar(LIGO_data['M1'], LIGO_data['M2'], yerr=LIGO_data['M2_error'], xerr=LIGO_data['M1_error'], fmt='o', color='black', ecolor='grey', zorder=1, label='LIGO data')
			ax[row][col].errorbar(LIGO_data['M1'], LIGO_data['M2'], yerr=LIGO_data['M2_error'], xerr=LIGO_data['M1_error'], fmt='o', color='black', ecolor='grey', zorder=1, label='LIGO data')
		M1 = alpha_data[i]['M1']
		M2 = alpha_data[i]['M2']
		ejected_M1 = alpha_data[i]['ejected_M1']
		ejected_M2 = alpha_data[i]['ejected_M2']
		# print(f"For alpha={alpha_values[i]}, M1 = {M1}")
		# print(f"For alpha={alpha_values[i]}, M2 = {M2}")
		ax[row][col].scatter(M1, M2, color='red', zorder=2, label='Simulation data')
		# ax[ax_idx].scatter(ejected_M1, ejected_M2, color='blue', zorder=3, label="Ejected black holes")
		ax[row][col].plot(x,x,'--', color='red', zorder=0)
		ax[row][col].set_ylim(0,y_lim)
		ax[row][col].set_xlim(0,x_lim)
		ax[row][col].set_xlabel(r'Mass of Primary ($M_{SUN}$)', fontsize='x-large')
		ax[row][col].set_ylabel(r'Mass of Secondary ($M_{SUN}$)', fontsize='x-large')
		ax[row][col].set_title(f'alpha = {alpha_values[i]}')
		# ax_idx += 1

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
# ax[0].legend(by_label.values(), by_label.keys(), loc='upper left')
f.legend(by_label.values(), by_label.keys())

f.suptitle(f'Mass of Secondary vs. Mass of Primary', fontsize='xx-large')		
f.savefig(PlotDir + f'/M2_vs_M1_WithEjectedSameColor_alpha={plot_label}.pdf', bbox_inches='tight')
print('M2 vs M1 complete')

######################################################################### M2 vs M1 (single plot)
'''
f, ax = plt.subplots(1,1,figsize=(8,7))

if IncludeLIGO:
	ax.errorbar(LIGO_data['M1'], LIGO_data['M2'], yerr=LIGO_data['M2_error'], xerr=LIGO_data['M1_error'], fmt='o', color='black', ecolor='gray', zorder=0, label='LIGO data')

for i in range(len(paths)):
	if paths[i] != '':
		M1 = alpha_data[i]['M1']
		M2 = alpha_data[i]['M2']
		ax.scatter(M1, M2, label=f'alpha={alpha_values[i]}', zorder=i+1)


x = np.linspace(0,x_lim,100)
ax.plot(x,x,'--')
ax.set_ylim(0,y_lim)
ax.set_xlim(0,x_lim)
ax.set_xlabel(r'Mass of Primary ($M_{SUN}$)', fontsize='x-large')
ax.set_ylabel(r'Mass of Secondary ($M_{SUN}$)', fontsize='x-large')
ax.set_title(f'Mass of Secondary vs. Mass of Primary', fontsize='xx-large')
ax.legend(loc='upper left')
f.savefig(PlotDir + f'/M2_vs_M1_singleplot_alpha={plot_label}.pdf', bbox_inches='tight')
print('M2 vs M1 (single plot) complete')

######################################################################### Mass Ratio vs M1

f, ax = plt.subplots(1,len(alphas_used),figsize=(6*len(alphas_used),7))
ax_idx = 0
for i in range(len(paths)):
	if paths[i] != '':
		if IncludeLIGO:
			M1 = LIGO_data['M1']
			M2 = LIGO_data['M2']
			# How to calculate y-error?
			ax[ax_idx].errorbar(M1, M2/M1, xerr=LIGO_data['M1_error'], fmt='o', color='black', ecolor='gray', zorder=0, label='LIGO data')
		M1 = alpha_data[i]['M1']
		M2 = alpha_data[i]['M2']
		ax[ax_idx].scatter(M1, M2/M1, color='red', zorder=1)
		ax[ax_idx].set_xlabel(r'Mass of Primary ($M_{SUN}$)', fontsize='x-large')
		ax[ax_idx].set_ylabel(r'$M2/M1$', fontsize='x-large')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}', fontsize='xx-large')
		ax_idx += 1
f.suptitle(f'Mass Ratio vs. Mass of Primary')		
f.savefig(PlotDir + f'/Mass_Ratio_vs_Primary_alpha={plot_label}.pdf', bbox_inches='tight')
print('Mass Ratio vs M1 complete')

######################################################################### Core Radius vs Half Mass Radius

# PLEASE NOTE that this code only takes every 100th value of both r_h and r_c to save time

# print('Now plotting Core Radius vs Half Mass Radius...')
f, ax = plt.subplots(1,len(alphas_used),figsize=(6*len(alphas_used),7))
ax_idx = 0
for i in range(len(paths)):
	if paths[i] != '':
		# print(f'alpha value = {alpha_values[i]}')
		r_c = alpha_data[i]['r_c'][::100]
		r_h = alpha_data[i]['r_h'][::100]
		# print(f"For alpha={alpha_values[i]}, M1 = {M1}")
		# print(f"For alpha={alpha_values[i]}, M2 = {M2}")
		ax[ax_idx].scatter(r_h, r_c)
		ax[ax_idx].set_xlabel('Half Mass Radius (pc)', fontsize='x-large')
		ax[ax_idx].set_ylabel('Core Radius (pc)', fontsize='x-large')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}', fontsize='xx-large')
		ax_idx += 1
f.suptitle(f'Core Radius vs. Half Mass Radius')		
f.savefig(PlotDir + f'/CoreRadius_vs_HalfMassRadius_scatter_alpha={plot_label}.pdf', bbox_inches='tight')
print('Core Radius vs Half Mass Radius (scatter) complete')

######################################################################### M2 vs M1

f, ax = plt.subplots(1,len(alphas_used),figsize=(6*len(alphas_used),7))
ax_idx = 0
for i in range(len(paths)):
	if paths[i] != '':
		# print('i =', i)
		r_c = alpha_data[i]['r_c'][::100]
		r_h = alpha_data[i]['r_h'][::100]
		# time = alpha_data[i]['time']
		# print(f"For alpha={alpha_values[i]}, M1 = {M1}")
		# print(f"For alpha={alpha_values[i]}, M2 = {M2}")
		ax[ax_idx].hist2d(r_h, r_c, bins=50, cmap='plasma')
		ax[ax_idx].set_xlabel('Half Mass Radius (pc)', fontsize='x-large')
		ax[ax_idx].set_ylabel('Core Radius (pc)', fontsize='x-large')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}', fontsize='xx-large')
		ax_idx += 1
f.suptitle(f'Core Radius vs. Half Mass Radius')		
f.savefig(PlotDir + f'/CoreRadius_vs_HalfMassRadius_2dhist_alpha={plot_label}.pdf', bbox_inches='tight')
print('Core Radius vs Half Mass Radius (2d histogram) complete')
'''