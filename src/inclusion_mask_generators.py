"""description of this module, interface.py(interface).
please update this docstring as you develop.

interface.py should contain functions that corresponde to each keyword in the "weighted_by" columns of the model_specification.csv. 
Each keyword needs to map to a function here that will return the right vector of n weights for n spines
"""

import numpy as np



def all_spines(param_array):
    weight_array = binary_weights(param_array, 100)
    return weight_array

def only(param_array):
	#the param arrays used with "only" should already be binary... 
	return param_array.astype(bool)

def top_20_percentile(param_array):
    weight_array = binary_weights(param_array, 20)
    return weight_array

def bottom_20_percentile(param_array):
    weight_array = binary_weights(param_array, 20, 'lowest')
    return weight_array

def top_50_percentile(param_array):
    weight_array = binary_weights(param_array, 50)
    return weight_array

def bottom_50_percentile(param_array):
    weight_array = binary_weights(param_array, 50, 'lowest')
    return weight_array


def binary_weights(param_array, threshold_percentage, direction = 'highest'):
    #if there are 100 spines and threshold_percentage, this will take the 20 spines with the highest params. NOT all spines with a param within 20% of the highest param.
    total_spines = len(param_array)
    threshold_n = round(total_spines*threshold_percentage/100)
    if direction == 'lowest':
        ind = np.argpartition(param_array, threshold_n)[:threshold_n]
    elif direction == 'highest':
        ind = np.argpartition(param_array, -threshold_n)[-threshold_n:]
    weights = np.zeros(len(param_array))
    weights[ind] = 1
    return weights