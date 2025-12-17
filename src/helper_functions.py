#intended to be functions that are used in the other src fils (computation.py, plotting.py) but avoids putting them in config.py
#each other module should import helper_functions as hf

import numpy as np
from src import config as cfg
from src import data_io as io

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#For organizational ease - first are the helper functions utilized to produce different weightings
#The functions are in the order they appear in the model_specifications.csv (weighted_by column first and then inclusion_based_on column)

#def weights_from_corr_to_nearest_resp_spine(spine_data):
#    return get_weight_matrix(spine_data, hf.nn_spine_stim_corr)

#def weights_from_spine_amp(spine_data):
#    return get_weight_matrix(spine_data, hf.spine_amplitude)

#def weights_from_spine_reliability(spine_data):
#    return get_weight_matrix(spine_data, hf.spine_reliability)

#def weights_from_spine_reliability(spine_data):
#    return get_weight_matrix(spine_data, hf.spine_reliability)

#def weights_from_spine_reliability(spine_data):
#    return get_weight_matrix(spine_data, hf.spine_reliability)

def get_soma_traces(soma_data):
    soma_field_2 = get_soma_activity_meta(soma_data)
    soma_traces = get_traces(soma_data)
    #should be spines x directions x presentations x samples
    start_idx, end_idx = get_sample_idxs(soma_field_2, cfg.start_s, cfg.end_s)
    soma_traces = select_timesteps(soma_traces, start_idx, end_idx)
    return soma_traces

def get_sample_idxs(fov_activity_meta, start_s, end_s):

    relative_frame_times = np.array(fov_activity_meta['trial_time'])

    def find_nearest(array, value): ## from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        array = np.asarray(array)
        idx = np.nanargmin(np.abs(array - value))
        return idx, array[idx]

    start_idx, actual_value = find_nearest(relative_frame_times, start_s)
    end_idx, actual_value = find_nearest(relative_frame_times, end_s)
    return start_idx, end_idx




def select_timesteps(traces, first_sample=0, last_sample=-1):
    if len(traces.shape) == 3:
        selected_timesteps = traces[:,:,first_sample:last_sample]
    else:
        selected_timesteps = traces[:,:,:,first_sample:last_sample]
    return selected_timesteps

def get_stim_on_traces(traces):
    return select_timesteps(traces, first_sample =cfg.stim_start, last_sample =cfg.stim_end )




def get_fovs_responsive_spines(fov_activity, fov_metadata):
    return list(fov_activity['responsive'])[0]


def get_fovs_unresponsive_spines(fov_activity, fov_metadata):
    return ~list(fov_activity['responsive'])[0]


def get_fovs_apical_spines(fov_activity, fov_metadata):
    ref = fov_metadata['structural_data']['Apical_Basal'][0][0]
    bool_apical = chr(fov_metadata[ref][0]).lower()=='a'
    num_spines = get_num_spines_in_fov(fov_metadata)
    bool_apical_list = [bool_apical]*num_spines
    return bool_apical_list


def get_fovs_basal_spines(fov_activity, fov_metadata):
    ref = fov_metadata['structural_data']['Apical_Basal'][0][0]
    bool_apical = chr(fov_metadata[ref][0]).lower()=='b'
    num_spines = get_num_spines_in_fov(fov_metadata)
    bool_apical_list = [bool_apical]*num_spines
    return bool_apical_list


def get_fovs_corr_to_nearest_resp_spine(fov_activity, fov_metadata):
    dist_corr_array = fov_activity['nn_dist_corr'][1]
    dist_corr_array = np.nan_to_num(dist_corr_array)
    return dist_corr_array


def get_fovs_dist_to_nearest_resp_spine(fov_activity, fov_metadata):
    dist_corr_array = fov_activity['nn_dist_corr'][0]
    dist_corr_array = np.nan_to_num(dist_corr_array)
    return dist_corr_array


def get_fovs_spine_amplitude(fov_activity, fov_metadata):
    return list(fov_activity['max'])[0]




def get_fovs_spine_preferred_orientation(fov_activity, fov_metadata):
    return list(fov_activity['pref_ori'])[0]

def get_fovs_spine_preferred_direction(fov_activity, fov_metadata):
    return list(fov_activity['pref_dir'])[0]

def get_fovs_spine_dsi(fov_activity, fov_metadata):
    return list(fov_activity['DSI'])[0]

def get_fovs_spine_median_resp(fov_activity, fov_metadata):
    return list(fov_activity['median_amp'])[0]




def get_fovs_spines_size(fov_activity, fov_metadata):
    return list(fov_metadata['spine_size'])[0] #strange syntax here because it needs to return a list... must be a better way but this works. Coudl we do the coersion to flat list on the other end? in the calling function?


def get_fovs_spine_distance_from_soma(fov_activity, fov_metadata):
    num_spines = get_num_spines_in_fov(fov_metadata)
    fov_dist = fov_metadata['structural_data']['DistanceFromRoot_um'][0][0]
    spines_dist = [fov_dist]*num_spines
    return spines_dist


def get_fovs_spine_neck_length(fov_activity, fov_metadata):
    return list(fov_metadata['neckLength'])[0]


def get_fovs_spine_osi(fov_activity, fov_metadata):
    return list(fov_activity['OSI'])[0]


def get_fovs_spine_dsi(fov_activity, fov_metadata):
    return list(fov_activity['DSI'])[0]


def get_fovs_spine_preferred_direction(fov_activity, fov_metadata):
    return list(fov_activity['pref_dir'])[0]


def get_fovs_spine_dsi(fov_activity, fov_metadata):
    return list(fov_activity['DSI'])[0]


def get_fovs_spine_amplitudes(fov_activity, fov_metadata):
    return np.array(fov_activity['median_amp'])



def get_fovs_spine_reliability(fov_activity, fov_metadata):
    #at its preferred direction
    all_reliabilities = np.array(fov_activity['reliability_index']) #spines x directions
    preferred_stim_inidces = list(fov_activity['max'])[1] - 1 #account for change in indexing!
    preferred_stim_inidces = np.array(preferred_stim_inidces).astype(int)
    spine_indices = np.arange(0,len(preferred_stim_inidces))
    pref_reliabilities = all_reliabilities[preferred_stim_inidces, spine_indices]
    return pref_reliabilities


def get_fovs_random_integers(fov_activity, fov_metadata):
    param_array = list(fov_activity['DSI'])[0]
    random_param_array = np.random.randint(0,high=len(param_array), size=len(param_array))
    return random_param_array


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






def nn_spine_trace_corr(fov_activity, fov_metadata):
    nearest_neighbor_list = get_spine_nn_idxs(fov_metadata)
    corr_mat = get_trace_corr(fov_activity)
    indicies = list(range(len(nearest_neighbor_list)))
    return corr_mat[indicies,nearest_neighbor_list]


def nn_spine_stim_corr(fov_activity, fov_metadata):
    nearest_neighbor_list = get_spine_nn_idxs(fov_metadata)
    corr_mat = get_stim_corr(fov_activity)
    indicies = list(range(len(nearest_neighbor_list)))
    return corr_mat[indicies,nearest_neighbor_list]

def get_spine_nn_idxs(fov_metadata):
    #find min for each spot in array
    pairwise_dist_array = np.array(fov_metadata['pairwise_dist'])

    #we need to replace the diagonal or it will always return the index of itself. 
    np.fill_diagonal(pairwise_dist_array, 5)#np.max(pairwise_dist_array))
    nearest_neighbor_list = np.argmin(pairwise_dist_array, axis=0)
    return nearest_neighbor_list


def get_trace_corr(fov_activity):
    return np.array(fov_activity['total_trial_pairwise_corr'])


def get_stim_corr(fov_activity):
    return np.array(fov_activity['corr_vis_periods_baps_excluded'])


def get_dist_from_root(spine_data, fov_num = 0):
    metadata_dict = get_spine_metadata(spine_data, fov_num)
    return metadata_dict['structural_data']['DistanceFromRoot_um'][0][0]

def get_dist_from_branch(spine_data, fov_num = 0):
    metadata_dict = get_spine_metadata(spine_data, fov_num)
    return metadata_dict['structural_data']['DistanceFromBranch_um'][0][0]


def get_spines_pixel_coords(spine_data, fov_num=0):
    metadata_dict = get_spine_metadata(spine_data, fov_num)
    return metadata_dict['spine_xy']


def get_direction_labels(spine_data, fov_num = 0):
    spine_activity = get_spine_activity(spine_data, fov_num)
    ref = spine_activity['OSI_DSI_angles'][0][1]
    return np.array(spine_data[ref])[:,0]


def get_trial_time_labels(spine_data, fov_num = 0):
    spine_activity = get_spine_activity(spine_data, fov_num)
    return np.array(spine_activity['trial_time'])[:,0]


def get_neck_length(spine_data, fov_num = 0):
    ref = spine_data['dend_cell'][3,fov_num]#['stem_stats']#['neckLength']
    spine_field_3 = spine_data[ref]
    return spine_field_3
    #spine_activity = get_spine_activity(spine_data, fov_num)
    #return np.array(spine_activity['trial_time'])[:,0]




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##################################################################################
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_model_names(model_spec, k=0):
    full_model_name = f"{model_spec['trials_included']}-{model_spec['spines_included']}-{model_spec['inclusion_based_on']}-{model_spec['weighted_by']}-{model_spec['integration']}-{model_spec['somatic_function']}"
    model_initials = ''.join([ s[0] for s in full_model_name.split('-') ])
    nickname =  f'model-{model_initials}-{k}'
    return full_model_name, nickname

def fovs_in_dataset(spine_data):
    num_fovs = spine_data['dend_cell'].shape[1]
    return num_fovs


def get_fov_zetas(fov_activity_meta):
    return fov_activity_meta['ZETA_test'][0]

def get_fov_idx_from_all_spine_idx(i, spine_data):
    mapping_dict = _all_idx_to_fov_idx_mapping(spine_data)
    return mapping_dict[i]


def _all_idx_to_fov_idx_mapping(spine_data):
    mapping_dict = {}
    spine_i = 0
    for fov_idx, (fov_activity, fov_meta) in enumerate(fov_generator(spine_data)):
        num_spines = get_num_spines_in_fov(fov_meta)
        for spine_fov_idx in range(num_spines):
            mapping_dict[spine_i] = (fov_idx, spine_fov_idx)
            spine_i +=1
    return mapping_dict

def get_num_spines_in_fov(fov_activity_meta):
    return len(list(fov_activity_meta['spine_size'])[0]) #spine_size, this is good and from kyle in imageJx - SpineArea is poorly computed by Gregg

def get_spine_activity(spine_data, fov_num = 0):
    ref = spine_data['dend_cell'][2,fov_num]
    spine_field_2 = spine_data[ref]
    return spine_field_2


def get_precomputed_tuning_curve(soma_data):
    soma_field_2 = get_soma_activity_meta(soma_data)
    return np.array(soma_field_2['median_amp'])

def include_soma(soma_data):
    return True
    #soma_field_3 = io._todict(soma_data[3])
    #return bool(soma_field_3['included'])


def get_responsive_status(soma_data):
    #get preferred orientation
    soma_field_2 = get_soma_activity_meta(soma_data)
    return bool(soma_field_2['responsive'])




def get_soma_responsive_status(soma_data):
    #get preferred orientation
    soma_field_2 = get_soma_activity_meta(soma_data)
    return bool(soma_field_2['responsive'])

def get_soma_preferred_orientation(soma_data):
    #get preferred orientation
    soma_field_2 = get_soma_activity_meta(soma_data)
    return soma_field_2['pref_ori']

def get_soma_preferred_direction(soma_data):
    #get preferred orientation
    soma_field_2 = get_soma_activity_meta(soma_data)
    return soma_field_2['pref_dir']

def get_soma_median_response_amplitudes(soma_data):
    #get preferred orientation
    soma_field_2 = get_soma_activity_meta(soma_data)
    return soma_field_2['median_amp']


def get_reliability_at_pref_stim(soma_data):
    med_resp_amps = get_soma_median_response_amplitudes(soma_data)
    preferred_stim_idx = np.argmax(med_resp_amps)
    
    soma_field_2 = get_soma_activity_meta(soma_data)
    reliabilities = soma_field_2['reliability_index']
    return reliabilities[preferred_stim_idx]


    #tuning_curve = get_precomputed_tuning_curve(soma_data)
    #preferred_stim_index = np.argmax(tuning_curve)
    #zeta_results = soma_field_2['ZETA_test_dir'][preferred_stim_index]
    #p_value = io._todict(zeta_results)['dblP']
    #print(p_value)
    #responsive_status = 'unresponsive'
    #if p_value<.05:
    #    responsive_status = 'responsive'
    #return responsive_status


def get_spine_metadata(spine_data, fov_num = 0):
    ref = spine_data['dend_cell'][3,fov_num]
    spine_field_3 = spine_data[ref]
    return spine_field_3


def fov_generator(spine_data):
    num_fovs = fovs_in_dataset(spine_data)
    for fov_num in range(num_fovs):
        yield get_spine_activity(spine_data, fov_num), get_spine_metadata(spine_data, fov_num)


def get_xyz_coords_of_fov(spine_data, fov_num = 0):
    metadata_dict = get_spine_metadata(spine_data, fov_num )
    coordinate_list = spine_data[metadata_dict['structural_data']['coordinate'][0][0]][:]
    new_list = [chr(i[0]) for i in coordinate_list][1:-1]
    string_nums = ''.join(new_list)[0:-1]

    num_list = []
    this_num = ''
    for char in string_nums:
        if char == ',':
            num_list.append(float(this_num))
            this_num = ''
        else:
            this_num += char
    num_list.append(float(this_num))
    return np.array(num_list)


def get_fov_name(spine_data, fov_num = 0):
    metadata_dict = get_spine_metadata(spine_data, fov_num)
    ascii_list = metadata_dict['seg_name']
    fov_name = ''.join(chr(i[0]) for i in ascii_list)
    return fov_name


def get_bap_trials(spine_data, fov_num = 0):
    metadata_dict = get_spine_metadata(spine_data, fov_num)
    return get_bap_trials_meta(metadata_dict)

def get_bap_trials_meta(metadata_dict):
    #print(metadata_dict.keys())
    return metadata_dict['BAP_trials']


def get_presentation_num(fov_activity_meta):
    return fov_activity_meta['trial_amp'].shape[0]


def get_stim_num(fov_activity_meta):
    return fov_activity_meta['trial_amp'].shape[1]


def get_branch_order(spine_data, fov_num = 0):
    metadata_dict = get_spine_metadata(spine_data, fov_num)
    return metadata_dict['structural_data']['order'][0][0]


def get_fov(activity_data_struct, fov=0):
    try:
        if fov%1==0:
            ref = activity_data_struct['dend_cell'][2,fov]
        else:
            ref = fov
    except Exception as E:
        ref = fov
    spine_field_2 = activity_data_struct[ref]#['DSI']
    return spine_field_2


def get_traces(soma_data, fov=0, spine_index=0):
    #tried to make it work for both soma and spine structure
    #Don't think this should be used for spines anymore so commenting out this part to avoid unexpected errors

    #try:
    #    this_fov = get_fov(activity_data_struct, fov=fov)
    #    try:
    #        spine_traces = np.array(this_fov['trial_traces'][:,:,0,:,spine_index].swapaxes(0,-1))
    #    except TypeError as E:
    #        spine_traces = np.array(io._todict(activity_data_struct[2])['trial_traces'][:,:,0,:,spine_index].swapaxes(0,-1))
    #        #spine_traces = np.array(this_fov['trial_traces'])#[:,:,0,:,spine_index].swapaxes(0,-1))
    #        #print(shape(spine_traces))
    #        #raise()
    #    traces = spine_traces

    #except IndexError as E:
    soma_field_2 = get_soma_activity_meta(soma_data)
    soma_traces = np.array(soma_field_2['trial_traces'])
    traces = soma_traces

    return traces

def get_soma_activity_meta(soma_data):
    return io._todict(soma_data[2])


def get_example_traces(spine_data):
    ref = spine_data['dend_cell'][2,0]
    fov_field_2 = spine_data[ref]
    ex_spine_trace = np.array(fov_field_2['trial_traces'][:,:,0,:,0].swapaxes(0,-1))
    return ex_spine_trace

def get_spine_global_coords(h, spine_data, fov_num, spine_idx, manual_adjstment=np.array([0,0,0])):
    xyz_coords = get_xyz_coords_of_fov(spine_data, fov_num = fov_num)
    #nearest_section, sec_coords, min_dist  = find_closest_section(h, xyz_coords)

    spine_pix_coords = get_spines_pixel_coords(spine_data, fov_num)[:,spine_idx]
    #delta_from_center =
    delta_pix_from_center = cfg.fov_pixel_dim/2 - spine_pix_coords
    delta_um_from_center = list(delta_pix_from_center*cfg.pixel_size)
    delta_um_from_center.append(0) #add the 3rd dimension
    spine_global_coords = -1*manual_adjstment+xyz_coords+np.array([-1,1,0])*np.array(delta_um_from_center)
    return spine_global_coords


def returnSegmentCoordinates(section):
    # Get section 3d coordinates and put in numpy array
    #from https://www.neuron.yale.edu/phpBB/viewtopic.php?t=1528
    #thank you landoscape
    n3d = section.n3d()
    x3d = np.empty(n3d)
    y3d = np.empty(n3d)
    z3d = np.empty(n3d)
    L = np.empty(n3d)
    for i in range(n3d):
        x3d[i]=section.x3d(i)
        y3d[i]=section.y3d(i)
        z3d[i]=section.z3d(i)

    # Compute length of each 3d segment
    for i in range(n3d):
        if i==0:
            L[i]=0
        else:
            L[i]=np.sqrt((x3d[i]-x3d[i-1])**2 + (y3d[i]-y3d[i-1])**2 + (z3d[i]-z3d[i-1])**2)

    # Get cumulative length of 3d segments
    cumLength = np.cumsum(L)

    N = section.nseg

    # Now upsample coordinates to segment locations
    xCoord = np.empty(N)
    yCoord = np.empty(N)
    zCoord = np.empty(N)
    dx = section.L / (N-1)
    for n in range(N):
        if n==N-1:
            xCoord[n]=x3d[-1]
            yCoord[n]=y3d[-1]
            zCoord[n]=z3d[-1]
        else:
            cIdxStart = np.where(n*dx >= cumLength)[0][-1] # which idx of 3d segments are we starting at
            cDistFrom3dStart = n*dx - cumLength[cIdxStart] # how far along that segment is this upsampled coordinate
            cFraction3dLength = cDistFrom3dStart / L[cIdxStart+1] # what's the fractional distance along this 3d segment
            # compute x and y positions
            xCoord[n] = x3d[cIdxStart] + cFraction3dLength*(x3d[cIdxStart+1] - x3d[cIdxStart])
            yCoord[n] = y3d[cIdxStart] + cFraction3dLength*(y3d[cIdxStart+1] - y3d[cIdxStart])
            zCoord[n] = z3d[cIdxStart] + cFraction3dLength*(z3d[cIdxStart+1] - z3d[cIdxStart])
    return xCoord, yCoord, zCoord


################################
#These ones are sorta doing computations...

def find_closest_section(h, xyz_coords):
    for section in h.allsec():
        for seg_coords in section.psection()['morphology']['pts3d']:
            dist = np.linalg.norm(xyz_coords - seg_coords[:-1]) #last point is radius
            try:
                if min_dist> dist:
                    min_dist = dist
                    nearest_section = section
                    sec_coords = seg_coords[:-1]
            except UnboundLocalError as E:
                min_dist = dist
                nearest_section = section
                sec_coords = seg_coords[:-1]

    return nearest_section, sec_coords, min_dist


def get_period_means(traces):
    selected_period_means = np.empty(traces.shape)

    #print(traces.shape)

    scale_width = 0
    #Reduce each time period trace to the mean in each period
    for i in range(cfg.num_tranges):
        mean_activity_in_trange = np.mean(traces[:,:,i*cfg.timepoints_per_period:(i+1)*cfg.timepoints_per_period], axis=2)
        #print(mean_activity_in_trange.shape)
        for j in range(cfg.timepoints_per_period):
            selected_period_means[:,:,j+i*cfg.timepoints_per_period] = mean_activity_in_trange #tried doing this with tile and ran into trouble/
    return selected_period_means


def get_bool_activity(traces):
    #Use this to determine the soma threshold
    #plt.hist(first_entry_after_stim)
    #spine_threshold = 2

    #then boolean over or under the soma threshold
    bool_activity = traces>cfg.soma_threshold
    return bool_activity




