import numpy as np
import matplotlib.pyplot as plt

DataDir = './data/'
PlotDir = './plots/'

alpha_values = ['1.6', '2.0', '2.3', '2.6', '3.0']
rates = [[] for i in range(len(alpha_values))]  		## list of rates of all mergers
in_cluster = [[] for i in range(len(alpha_values))]  	## list of in-cluster merger rates
out_cluster = [[] for i in range(len(alpha_values))]  	## list of out-of-cluster merger rates

in_redshifts = [[] for i in range(len(alpha_values))]
out_redshifts = [[] for i in range(len(alpha_values))]
redshifts = [[] for i in range(len(alpha_values))]

# # Reorder the combined redshifts and rates
# for i in range(len(alpha_values)):
# 	idx = np.argsort(redshifts[i])
# 	redshifts[i] = np.array(redshifts[i])[idx]
# 	rates[i] = np.array(rates[i])[idx]

for i in range(len(alpha_values)):
	for j in range(len(redshifts[i])):
		new_idx = 0
		keep_going = True
		while keep_going:
			if redshifts[new_idx] <= redshifts[j] and redshifts[new_idx+1] >= redshifts[j]:
				## Change location of redshift
				z = redshifts[j]
				redshifts.insert(new_idx, z)
				redshifts.remove(z)

				## Change location of rate
				rate = rates[i][j]
				rates.insert(new_idx, rate)
				redshifts.remove(rate)

				keep_going = False
			else:
				new_idx +=1

# file_types = []
# for i in range(len(alpha_values)):
#     file_types.append(f"in_cluster_alpha={alpha_values[i]}")
#     file_types.append(f"out_cluster_alpha={alpha_values[i]}")

for i in range(len(alpha_values)):
	in_cluster[i] = np.load(f'{DataDir}cumulative_rates/in_cluster_alpha={alpha_values[i]}.npy')
	# print(f'len(in_cluster) = {len(in_cluster[i])}')
	out_cluster[i] = np.load(f'{DataDir}cumulative_rates/out_cluster_alpha={alpha_values[i]}.npy')
	# print(f'len(out_cluster) = {len(out_cluster[i])}')
	rates[i] = np.append(in_cluster[i], out_cluster[i])
	# print(f'len(rates[i]) = {len(rates[i])}')
	# print(f'rates[i] = {rates[i]}\n')
	in_redshifts[i] = np.load(f'{DataDir}redshifts/in_cluster_alpha={alpha_values[i]}.npy')
	out_redshifts[i] = np.load(f'{DataDir}redshifts/out_cluster_alpha={alpha_values[i]}.npy')
	redshifts[i] = np.append(in_redshifts[i], out_redshifts[i])

	# print(f'\nalpha = {alpha_values[i]}')
	# print('The following two lengths should be the same')
	# print(f'len(rates) = {len(rates[i])}')
	# print(f'len(redshifts) = {len(redshifts[i])}')
	
	# print('\nin_cluster[i]:')
	# print(in_cluster[i])
	# print('\nout_cluster[i]:')
	# print(out_cluster[i])

	# print('\nin_redshifts[i]:')
	# print(in_redshifts[i])
	# print('\nout_redshifts[i]:')
	# print(out_redshifts[i])

# Plotting
#########################################################################################
# f, ax = plt.subplots(1,5,figsize=(35,7))
# for i in range(len(rates)):
# 	ax[i].hist(rates[i])
# 	ax[i].set_xlabel('Cumulative Rate of BBH Mergers [yr^-1]', fontsize='x-large')
# 	ax[i].set_title(f'alpha={alpha_values[i]}', fontsize='x-large')

# f.savefig(PlotDir + f'/CumulativeRate.pdf', bbox_inches='tight')

# #########################################################################################
# f2, ax2 = plt.subplots(2,5,figsize=(35,14))
# f2.suptitle('Cumulative Rate of BBH Mergers', fontsize='xx-large')
# for i in range(len(rates)):
# 	ax2[0][i].hist(in_cluster[i], color='mediumblue')
# 	ax2[1][i].hist(out_cluster[i], color='red')
# 	ax2[0][i].set_xlabel('Merger Rate in Clusters [yr^-1]', fontsize='x-large')
# 	ax2[1][i].set_xlabel('Merger Rate of Ejected BBHs [yr^-1]', fontsize='x-large')
# 	ax2[0][i].set_title(f'alpha={alpha_values[i]}', fontsize='x-large')

# f2.savefig(PlotDir + f'/CumulativeRate_SeparateEjected.pdf', bbox_inches='tight')

#########################################################################################
# for the cumulative rates I would try to mimic the one I sent you, where as you move in
# z, you keep adding the rates ( += ) and then just plot that number as a function of z

in_rate_sums = [[] for i in range(len(alpha_values))]
out_rate_sums = [[] for i in range(len(alpha_values))]

# Reverse the order of rate and redshift arrays
# This will allow us to add rates as z increases (z is currently decreasing)
# for i in range(len(alpha_values)):
# 	in_cluster = np.flip(in_cluster)
# 	out_cluster = np.flip(out_cluster)
# 	in_redshifts = np.flip(in_redshifts)
# 	out_redshifts = np.flip(out_redshifts)

for i in range(len(alpha_values)):
	current_in_sum = 0
	current_out_sum = 0
	for j in range(len(in_redshifts[i])):
		# print(f'i={i}, j={j}')
		current_in_sum += in_cluster[i][j]
		in_rate_sums[i] = np.append(in_rate_sums[i], current_in_sum)
	for k in range(len(out_redshifts[i])):
		current_out_sum += out_cluster[i][j]
		out_rate_sums[i] = np.append(out_rate_sums[i], current_out_sum)
	# print(f'\nrate_sums[i] = {rate_sums[i]}')

f3, ax3 = plt.subplots(2,5,figsize=(42,14))
for i in range(len(alpha_values)):
	ax3[0][i].plot(np.flip(in_redshifts[i]), in_rate_sums[i])
	ax3[1][i].plot(np.flip(out_redshifts[i]), out_rate_sums[i])
	ax3[0][i].set_ylabel(f'Cumulative Rate of BBH Mergers [yr^-1]', fontsize='x-large')
	ax3[1][i].set_ylabel(f'Cumulative Rate of BBH Mergers [yr^-1]', fontsize='x-large')
	ax3[0][i].set_xlabel('z', fontsize='x-large')
	ax3[1][i].set_xlabel('z', fontsize='x-large')
	ax3[0][i].set_title(f'In Cluster, alpha={alpha_values[i]}', fontsize='x-large')
	ax3[1][i].set_title(f'Out Of Cluster, alpha={alpha_values[i]}', fontsize='x-large')
	ax3[0][i].set_yscale('log')
	ax3[1][i].set_yscale('log')
	ax3[0][i].set_ylim(top=10**7)
	ax3[1][i].set_ylim(top=10**7)

# f3.yscale('log')
f3.suptitle('Merger Rate vs Redshift', fontsize='xx-large')
f3.savefig(PlotDir + f'/SeparatedRateVsRedshift.pdf', bbox_inches='tight')

#########################################################################################
# for the cumulative rates I would try to mimic the one I sent you, where as you move in
# z, you keep adding the rates ( += ) and then just plot that number as a function of z

# rate_sums = [[] for i in range(len(alpha_values))]
# for i in range(len(alpha_values)):
# 	rate_sums[i] = np.append(in_rate_sums[i], out_rate_sums[i])
# 	idx = np.argsort(redshifts[i])
# 	redshifts[i] = np.array(redshifts[i])[idx]
# 	rate_sums[i] = np.array(rate_sums[i])[idx]

rate_sums = [[] for i in range(len(alpha_values))]
for i in range(len(alpha_values)):
	current_sum = 0
	for j in range(len(redshifts[i])):
		current_sum += rates[i][j]
		rate_sums[i] = np.append(rate_sums[i], current_sum)

f4, ax4 = plt.subplots(1,5,figsize=(35,7))
for i in range(len(alpha_values)):
	ax4[i].plot(redshifts[i], rate_sums[i])
	ax4[i].set_yscale('log')
	ax4[i].set_ylabel(f'Cumulative Rate of BBH Mergers [yr^-1]', fontsize='x-large')
	ax4[i].set_xlabel('z', fontsize='x-large')
	ax4[i].set_title(f'alpha={alpha_values[i]}', fontsize='x-large')

f4.suptitle('Merger Rate vs Redshift')
f4.savefig(PlotDir + f'/RateVsRedshift.pdf', bbox_inches='tight')
