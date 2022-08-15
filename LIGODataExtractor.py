import numpy as np
import pandas as pd

# This was just a place to write code
# This code is now a part of BHMergerPlotter.py

data = pd.read_csv('LIGO_data.csv')
print(data.columns)

M1 = data['mass_1_source']
M1_low = data['mass_1_source_lower']
M1_up = data['mass_1_source_upper']
M2 = data['mass_2_source']
M2_low = data['mass_2_source_lower']
M2_up = data['mass_2_source_upper']

data = {'M1': M1
		'M1_error': [M1_low, M1_up],
		'M2': M2
		'M2_error': [M2_low, M2_up]}