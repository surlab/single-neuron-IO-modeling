"""description of this module, weight_generators.py(interface).
please update this docstring as you develop.

weight_generators.py should contain functions that corresponde to each keyword in the "weighted_by" columns of the model_specification.csv. 
Each keyword needs to map to a function here that will return the right vector of N weights for N spines
Each weight function takes the spine metadata as an input and returns this length N vector
The functions are in the order they appear in the model_specifications.csv (weighted_by column first and then inclusion_based_on column)
"""
from src import helper_functions as hf
from src import config as cfg
import numpy as np
from src import trial_selection

def weight_identity(array):
    return 1.0*np.ones(array.shape)

def weight_direct_from_value(array):
    return array

def weight_from_normalized_nonzero_value(weight_array):

    #Need to replace NANs with 0
    nonzero_weight_array = weight_array - np.nanmin(weight_array)
    nonnan_weight_array = np.nan_to_num(nonzero_weight_array, nan=0.0, posinf=0.0, neginf=0.0)
    #print('nonnan_weight_array', nonnan_weight_array)
    #we only want to normalize if its linear (and makes sense to normalize within cell, not within FOV)
    #but actually this won't affect binary weights - just dividing by 1.
    normalized_weight_array = (nonnan_weight_array)/(np.max(nonnan_weight_array))
    return normalized_weight_array

def get_weight_array(spine_data, param_func, weight_func=weight_from_normalized_nonzero_value, fov_based=False):
    param_list = []
    weight_iterable = []
    for i, (fov_activity, fov_metadata )in enumerate(hf.fov_generator(spine_data)):
        fov_params = param_func(fov_activity, fov_metadata)
        try:
            fov_params = fov_params.tolist() #to use extend reliably later this needs to be a list
        except AttributeError:
            pass #if we can't make it a list this way its probably already a list, otherwise let it fail on the extend

        if fov_based:
            weight_array = weight_func(np.array(fov_params))
            weight_iterable.extend(weight_array.astype(list))
        else:
            param_list.extend(list(fov_params))

    
    if not(fov_based):
        weight_iterable = weight_func(np.array(param_list).astype(np.float32))
    weight_array = np.array(weight_iterable)
    assert(np.max(weight_array) <= 1.0, 'The weights should be max 1.0')
    assert(np.min(weight_array) >= 0.0, 'The weights should be min 0.0')
    return weight_array


def weights_from_democratic_weights(spine_data):
    return get_weight_array(spine_data, hf.get_fovs_spine_distance_from_soma, weight_identity)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#these ones need to be pplied with the responsive spine mask

def weights_from_dist_to_nearest_resp_spine(spine_data):
    return get_weight_array(spine_data, hf.get_fovs_dist_to_nearest_resp_spine, weight_from_distance) #spines closer to other responsive spines may be amplified so we need to invert this. 





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#these ones are applied to all spines


##################################################################################
def weights_from_spine_size(spine_data):
    return get_weight_array(spine_data, hf.get_fovs_spines_size, weight_from_size, fov_based=True)#fov based normalization of max is supported from kasai and also accoutns for varied imaging depth/quality

def weight_from_size(size_array):
    #y = mx+b
    b = .2
    m = .8
    weights = size_array/np.max(size_array)*m+b
    return weights


##################################################################################
def weights_from_spine_distance_from_soma(spine_data):
    weight_array = get_weight_array(spine_data, hf.get_fovs_spine_distance_from_soma, weight_from_distance)
    nonnan_weight_array = np.nan_to_num(weight_array, nan=0.0, posinf=0.0, neginf=0.0)
    return nonnan_weight_array
                 
#def weight_from_distance_lin(dist_array):
#    #y = mx+b
#    b = 1
#    m = -.002
#    return dist_array*m+b

def weight_from_distance(dist_array):
    #y = span1*e**(-x/k_1) + span2*e**(-x/k_2)
    span1 = 1/7
    k_1 = 20
    span2 = 6/7
    k_2 = 250
    return span1*np.exp(-dist_array/k_1) + span2*np.exp(-dist_array/k_2)


##################################################################################
def weights_from_size_AND_dist(spine_data):
	dist_weights = weights_from_distance_from_soma(spine_data)
	size_weights = weights_from_spine_size(spine_data)
	return dist_weights*size_weights #(elementwise multiplication)



##################################################################################
#def weights_from_neck_len_lin(spine_data):
#    return get_weight_array(spine_data, hf.get_fovs_spine_neck_length, weight_from_neck_len)
#
#def weight_from_neck_len(size_array):
#    #y = mx+b
#    b = 1
#    m = -.8
#    return size_array*m+b
##################################################################################


def weights_from_one_over_peak_height(spine_data):
    spine_data_xr, spines_per_fov_list = trial_selection.compile_spine_traces(spine_data)
    weight_array = 1/compute_peak_heights(spine_data_xr)
    return weight_from_normalized_nonzero_value(weight_array)

def compute_peak_heights(spine_data_xr):
    #Exclude the first 4 seconds so we don't double count those
    #skip this for now... may actually not be necessary. I think I'm picking the traces so I have 1 second before and 3 seconds after, although might just be the plotting

    #could smooth the traces first...this would help with shot noise. but they are already highly averaged...
    #Take the max amplitude for each presentation
    maxs = spine_data_xr.max(dim='samples')
   
    
    peak_height_list = []
    #find the top 1%, this needs to be done for each spine and take the mean
    for one_spines_maxs in maxs.transpose("spines", ...):
        #print(one_spines_maxs.spines.item()) #use .item to get the label
        #print(one_spines_maxs.shape)
        
        
        maxes_1_dim = np.array(one_spines_maxs).flatten()

        num_to_average = int(maxes_1_dim.size*cfg.norm_peaks_percent/100)
        
        top_ind = np.argpartition(maxes_1_dim, -num_to_average)[-num_to_average:]
        
        top_maxs = maxes_1_dim[top_ind]
        peak_height_list.append(np.median(top_maxs))

    return np.array(peak_height_list)












