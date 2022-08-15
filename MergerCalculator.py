import numpy as np
import scipy.integrate as integrate

def get_merger_rate(z_lim):
	'''
	Returns the cumulative BBH merger rate.
	Input:
		z_lim - given redshift value (float)
	'''
	dVcdz = ## comoving volume at redshift z_lim
	rho_CG = 2.31 ## volumetric number density of clusters (units: Mpc^-3)
	dNdt = ## number of events in a unit time at a given redshift z_lim
	
	## Run these two lines to get the script to work before the actual values are determined
	# dVcdz = 2
	# dNdt = 3
	
	R = rho_CG * dNdt
	rate_tuple = integrate.quad(lambda z: R * (dVcdz) * pow(1+z, -1), 0, z_lim)
	rate = rate_tuple[0] ## the integrate.quad() function returns a tuple in the form: (estimated result, upper error bound)
	return rate

print(get_merger_rate(4))