"""description of this module, interface.py(interface).
please update this docstring as you develop.

interface.py should contain functions that corresponde to each keyword in the "weighted_by" columns of the model_specification.csv. 
Each keyword needs to map to a function here that will return the right vector of n weights for n spines
"""



weight_matricies_dict['size_AND_dist'] = weight_matricies_dict['spine_size'] * weight_matricies_dict['distance_from_soma']
#weight_matricies_dict['size_dist_resp_only'] = weight_matricies_dict['size_AND_dist']*weight_matricies_dict['resp_only']
weight_matricies_dict['nn_resp_spine_dist_AND_corr'] = weight_matricies_dict['corr_to_nn_resp_spine']*weight_matricies_dict['dist_to_nn_resp_spine']
