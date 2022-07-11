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
# path0 = ''
path1 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.0/0_0'
path2 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.3/0_0'
# path2 = ''
path3 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-2.6/0_0'
# path4 = ''
path4 = '/projects/b1095/newlin/cmc/IMFgrid2/rundir/rv1rg8z0.0017n8e5/w5fb0.05fbh0.05alpha3-3.0/0_0' #THIS IS STILL RUNNING

paths = [path0, path1, path2, path3, path4] # Add more when I receive them
PlotDir = './plots'

# load and format data
IncludeLIGO = True

alpha1_6 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0} # correct these alpha values if need be
alpha2_0 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0}
alpha2_3 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0}
alpha2_6 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0}
alpha3_0 = {'M1': 0, 'M2': 0, 'r_c': 0, 'r_h': 0, 'time': 0}
alpha_data = [alpha1_6, alpha2_0, alpha2_3, alpha2_6, alpha3_0]
alpha_values = ['1.6', '2.0', '2.3', '2.6', '3.0']
alphas_used = []

print('Loading data...')

for i in range(len(paths)):
	if paths[i] != '':
		M1, M2 = pd.read_csv(paths[i]+'/initial.bhmerger.dat',skiprows=1,header=0,delim_whitespace=True,dtype='str',usecols=(5,6)).values.T
		r_c, r_h, time = pd.read_csv(paths[i]+'/initial.dyn.dat',skiprows=1,header=0,delim_whitespace=True,dtype='str',usecols=(7,20,0)).values.T
		M1 = change_to_floats(M1)
		M2 = change_to_floats(M2)
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
			if M2[j] > M1[j]:
				larger = M2[j]
				M2[j] = M1[j]
				M1[j] = larger

		alpha_data[i]['M1'] = M1
		alpha_data[i]['M2'] = M2
		alpha_data[i]['r_c'] = r_c
		alpha_data[i]['r_h'] = r_h
		alpha_data[i]['time'] = time

		# print(f"For alpha={alpha_values[i]}, the average r_c is {np.mean(r_c)}")
		# print(f"For alpha={alpha_values[i]}, the average r_h is {np.mean(r_h)}")
		# print('len(r_c) =', len(r_c))
		# print('len(r_h) =', len(r_h))

		# print(alpha_data[i])

		alphas_used.append(alpha_values[i])

print('Simulation data loaded...')

if IncludeLIGO:
	LIGO_data = pd.read_csv('LIGO_data.csv')
	# print(LIGO_data.columns)

	M1 = LIGO_data['mass_1_source']
	M1_low = -1.0 * LIGO_data['mass_1_source_lower']
	M1_up = LIGO_data['mass_1_source_upper']
	M2 = LIGO_data['mass_2_source']
	M2_low = -1.0 * LIGO_data['mass_2_source_lower']
	M2_up = LIGO_data['mass_2_source_upper']

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

ax[0].set_xlabel(r'Mass of Primary ($M_{SUN}$)')
ax[1].set_xlabel(r'Mass of Secondary ($M_{SUN}$)')

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
f.legend(by_label.values(), by_label.keys())

f.suptitle('Primary and Secondary Masses of Black Hole Mergers')
f.savefig(PlotDir + f'/MassDistributions_alpha={plot_label}.pdf', bbox_inches='tight')
print('Mass Histograms complete')

######################################################################### M2 vs M1

f, ax = plt.subplots(1,len(alphas_used),figsize=(6*len(alphas_used),7))
ax_idx = 0
x = np.linspace(0,90,90)
for i in range(len(paths)):
	if paths[i] != '':
		# print('i =', i)
		if IncludeLIGO:
			ax[ax_idx].errorbar(LIGO_data['M1'], LIGO_data['M2'], yerr=LIGO_data['M2_error'], xerr=LIGO_data['M1_error'], fmt='o', color='black', ecolor='grey', zorder=1, label='LIGO data')
		M1 = alpha_data[i]['M1']
		M2 = alpha_data[i]['M2']
		# print(f"For alpha={alpha_values[i]}, M1 = {M1}")
		# print(f"For alpha={alpha_values[i]}, M2 = {M2}")
		ax[ax_idx].scatter(M1, M2, color='red', zorder=2, label='Simulation data')
		ax[ax_idx].plot(x,x,'--', color='red', zorder=0)
		ax[ax_idx].set_ylim(0,90)
		ax[ax_idx].set_xlim(0,140)
		ax[ax_idx].set_xlabel(r'Mass of Primary ($M_{SUN}$)')
		ax[ax_idx].set_ylabel(r'Mass of Secondary ($M_{SUN}$)')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}')
		ax_idx += 1

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
f.legend(by_label.values(), by_label.keys())

f.suptitle(f'Mass of Secondary vs. Mass of Primary')		
f.savefig(PlotDir + f'/M2_vs_M1_alpha={plot_label}.pdf', bbox_inches='tight')
print('M2 vs M1 complete')

######################################################################### M2 vs M1 (single plot)

f, ax = plt.subplots(1,1,figsize=(8,7))

if IncludeLIGO:
	ax.errorbar(LIGO_data['M1'], LIGO_data['M2'], yerr=LIGO_data['M2_error'], xerr=LIGO_data['M1_error'], fmt='o', color='black', ecolor='gray', zorder=0, label='LIGO data')

for i in range(len(paths)):
	if paths[i] != '':
		M1 = alpha_data[i]['M1']
		M2 = alpha_data[i]['M2']
		ax.scatter(M1, M2, label=f'alpha={alpha_values[i]}', zorder=i+1)


x = np.linspace(0,90,90)
ax.plot(x,x,'--')
ax.set_ylim(0,90)
ax.set_xlim(0,140)
ax.set_xlabel(r'Mass of Primary ($M_{SUN}$)')
ax.set_ylabel(r'Mass of Secondary ($M_{SUN}$)')
ax.set_title(f'Mass of Secondary vs. Mass of Primary')
ax.legend()
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
		ax[ax_idx].set_xlabel(r'Mass of Primary ($M_{SUN}$)')
		ax[ax_idx].set_ylabel(r'$M2/M1$')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}')
		ax_idx += 1
f.suptitle(f'Mass Ratio vs. Mass of Primary')		
f.savefig(PlotDir + f'/Mass_Ratio_vs_Primary_alpha={plot_label}.pdf', bbox_inches='tight')
print('Mass Ratio vs M1 complete')

######################################################################### Core Radius vs Half Mass Radius
'''
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
		ax[ax_idx].set_xlabel('Half Mass Radius (pc)')
		ax[ax_idx].set_ylabel('Core Radius (pc)')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}')
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
		ax[ax_idx].set_xlabel('Half Mass Radius (pc)')
		ax[ax_idx].set_ylabel('Core Radius (pc)')
		ax[ax_idx].set_title(f'alpha = {alpha_values[i]}')
		ax_idx += 1
f.suptitle(f'Core Radius vs. Half Mass Radius')		
f.savefig(PlotDir + f'/CoreRadius_vs_HalfMassRadius_2dhist_alpha={plot_label}.pdf', bbox_inches='tight')
print('Core Radius vs Half Mass Radius (2d histogram) complete')
'''