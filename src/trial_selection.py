import numpy as np
import xarray as xr


from src import helper_functions as hf
from src import config as cfg

def all_trials(fov_activity_meta):
    bap_trials = baps_trials_only(fov_activity_meta)
    return np.ones(bap_trials.shape)

def baps_trials_only(fov_activity_meta):
    return np.array(hf.get_bap_trials_meta(fov_activity_meta)).astype(int)

def no_bap_trials(fov_activity_meta):
    bap_trials = baps_trials_only(fov_activity_meta)    
    return -1*bap_trials




def get_subset_mask(fov_activity_meta, mask_func):
    trial_param = mask_func(fov_activity_meta)
    trials = hf.get_stim_num(fov_activity_meta)
    presentations = hf.get_presentation_num(fov_activity_meta)
    reshaped_trial_params = np.reshape(trial_param, (trials,presentations))
    #print(reshaped_trial_params.shape)
    #Should be directions x presntations - verified. As long as downstream uses correctly
    return reshaped_trial_params


def mask_traces(weights, traces):
    #weights is stims x presentations
    #traces is coming in as spines x stims x presentations x samples
    #and its just a numpy array, not an Xarray.
    tile_dims = (traces.shape[3], traces.shape[0], 1,1) #( samples x spines)
    tiled_weights = np.tile(weights, tile_dims )
    tiled_weights = np.moveaxis(tiled_weights, 0,3)
    masked_traces = traces*tiled_weights
    return masked_traces


def compile_spine_traces(spine_data, mask_func = all_trials):
    #idea is to put them all in a single array that can be sliced quickly
    #dimensison are spines x stims x trials x timestamps

    #we want to be able to slice this so that for a certain spine (all spines in an FOV)
    #we grab the Nth trial - all timstams and all stims
    #ideally we could do this in a single slice...
    #we will be able to do this using this notation
    #######
    #a = np.zeros((5, 16, 10, 91))
    #print(a.shape)
    ##lets try to grab the 1st presentation for spine 1 and the second for spine 2
    #sliced = a[np.arange(0, a.shape[0], 1), :,(0,2,4,6,7) , :]
    #print(sliced.shape)
    ########

    activity_list = []
    spines_per_fov_list = []
    spine_labels = []
    for i, (fov_activity_meta, fov_metadata )in enumerate(hf.fov_generator(spine_data)):
        fov_activity = np.array(fov_activity_meta['trial_traces'][:,:,0,:,:])

        fov_activity = fov_activity.swapaxes(0,-1)  #send samples to the last axes
        fov_activity = fov_activity.swapaxes(  1, 2) #bring directions in front of presentations
            #Now should be spines x directions x presentations x samples

        bap_trials = get_subset_mask(fov_activity_meta, mask_func)
        #Should be directions x presntations

        fov_activity_masked = mask_traces(bap_trials, fov_activity)
        #Make sure you check this

        #now grab only the frames of interest
        start_idx, end_idx = hf.get_sample_idxs(fov_activity_meta, cfg.start_s, cfg.end_s)
        #print(start_idx, end_idx)
        fov_activity_masked_subset = hf.select_timesteps(fov_activity_masked, start_idx, end_idx)
        #print('##')
        #print(np.shape(fov_activity_masked_subset))

        activity_list.extend(list(fov_activity_masked_subset))

        spines_in_fov = len(fov_activity_masked_subset)
        spines_per_fov_list.append(spines_in_fov)

        these_spine_labels = [f'fov_{i}_spine_{j}' for j in range(spines_in_fov)]

        spine_labels.extend(these_spine_labels)

    #print('####')
    #print(type(activity_list))
    #print(len(activity_list))
    all_spine_activity_array = np.array(activity_list) #this is spines x directions x presentations x samples


    #here seems like a good time to make this into an xarray as well...
    #except it looks like xarray won't let me slice this way. So would have to go to numpy and back again...
    #spine_labels = ['fov_name'+str(i) for i in range(a.shape[0])]
    direction_labels = hf.get_direction_labels(spine_data, fov_num = 0)
    sample_labels = hf.get_trial_time_labels(spine_data, fov_num = 0)[start_idx:end_idx]

    data_xr = xr.DataArray(all_spine_activity_array,
    coords={'spines': spine_labels,'directions': direction_labels,'samples': sample_labels},
    dims=["spines", "directions", "presentations", "samples"])
    ######
    return data_xr, spines_per_fov_list