import numpy as np
from scipy.stats import lognorm

def get_mass_weights(masses):
	'''
	Takes a numpy array of cluster masses (masses) and returns an array of mass weights (mass_weights).
	'''
	mass_weights = masses**(-2)
	return mass_weights

def get_metallicity_weights(metallicities): ## Finish this! CHECK IT!!!
	'''
	Takes a numpy array of cluster metallicities in units of solar metallicities (metallicities)
	and returns an array of metallicity weights based on a log normal distribution (Z_weights).
	'''
	sigma = 10**0.5 							## CONFIRM THIS VALUE
	log_mean = 0.153 - 0.074*(z**1.34)			## CONFIRM THIS VALUE
	mean = 10 ** mean
	Z_weights = lognorm.pdf(metallicities, sigma, scale=exp(mean))
	return Z_weights

def get_M_Z_weights(masses, metallicities):
	'''
	Takes a numpy array of cluster masses (masses) and returns an array of weights based on the mass and metallicity (M_Z_weights).
	'''
	mass_weights = get_mass_weights(masses)
	metallicity_weights = get_metallicity_weights(masses)
	M_Z_weights = mass_weights * metallicity_weights
	return M_Z_weights


def get_detection_weights(): ## Finish this!
	'''
	'''