"""description of this module, interface.py(interface).
please update this docstring as you develop.

interface.py should contain functions that corresponde to each keyword in the "weighted_by" columns of the model_specification.csv. 
Each keyword needs to map to a function here that will return the right vector of n weights for n spines
"""
import numpy as np

def summation_integration(traces, metadata_df):
    return np.sum(traces, axis=0) #should be directions x samples



def two_layer_by_fov_integration(traces, spines_per_fov_list, threshold=1):
    num_spines, num_directions, num_samples = traces.shape

    input_to_soma_trace = np.zeros([num_directions, num_samples])
    current_spine = 0
    for fov_num, num_spines in enumerate(spines_per_fov_list):
        activation_of_fov = np.sum(traces[current_spine:+current_spine+num_spines,:,:], axis=0)/num_spines
        output_of_fov = activation_of_fov - threshold
        output_of_fov = np.maximum(output_of_fov, 0)
        input_to_soma_trace += output_of_fov
        current_spine += num_spines
    return input_to_soma_trace #should be directions x samples


def two_layer_by_fov_1_integration(traces, spines_per_fov_list):
    return two_layer_by_fov_integration(traces, spines_per_fov_list, threshold=1)

def two_layer_by_fov_half_integration(traces, spines_per_fov_list):
    return two_layer_by_fov_integration(traces, spines_per_fov_list, threshold=.5)

def two_layer_by_fov_quarter_integration(traces, spines_per_fov_list):
    return two_layer_by_fov_integration(traces, spines_per_fov_list, threshold=.25)

def two_layer_by_fov_tenth_integration(traces, spines_per_fov_list):
    return two_layer_by_fov_integration(traces, spines_per_fov_list, threshold=.1)