
import numpy as np
from src import config as cfg
from src import helper_functions as hf
import pandas as pd
from scipy import stats 
import ipdb 
import xarray as xr

import statsmodels.api as sm
from statsmodels.formula.api import ols









def shuffle_along_axis(a, axis): #from https://stackoverflow.com/questions/5040797/shuffling-numpy-array-along-a-given-axis
    idx = np.random.rand(*a.shape).argsort(axis=axis)
    return np.take_along_axis(a,idx,axis=axis)


def shuffle_spine_traces(spine_traces):
    #right now this shuffles without regard to fov... may want to change this in the future?
    rng = np.random.default_rng()
    return rng.permuted(spine_traces, axis=1)





def init_traces_xarray(all_spine_activity_array, simulated_trials_per_stim):
    num_directions = len(all_spine_activity_array['directions'])
    num_samples= len(all_spine_activity_array['samples'])
    model_output = xr.DataArray(np.zeros((num_directions, simulated_trials_per_stim, num_samples)),
    coords={'directions': all_spine_activity_array['directions'],'samples': all_spine_activity_array['samples']},
    dims=["directions", "presentations", "samples"])
    return model_output



def sample_trial_from_fov(all_spine_activity_array, spines_per_fov_list):
    num_presentations = len(all_spine_activity_array['presentations'])
    #print(len(spines_per_fov_list))
    #draw a random integer for each fov
    rand_trial = np.random.randint(0, high=num_presentations, size=(len(spines_per_fov_list)))
    rand_trials_for_spines = list(np.repeat(rand_trial, spines_per_fov_list))
    spine_indicies = np.arange(0, len(all_spine_activity_array['spines']), 1)
    numpyfied_activity = np.array(all_spine_activity_array)
    simulated_trail_traces = numpyfied_activity[spine_indicies, :, rand_trials_for_spines, :]

    simulated_trail_traces = xr.DataArray(simulated_trail_traces,
        coords={'spines': all_spine_activity_array['spines'], 'directions': all_spine_activity_array['directions'],'samples': all_spine_activity_array['samples']},
        dims=["spines", "directions", "samples"])
    return simulated_trail_traces #should be spines x directions x samples



def sample_trials_from_fov(all_spine_activity_array, spines_per_fov_list, draw_trials):
    num_presentations = len(all_spine_activity_array['presentations'])
    #print(len(spines_per_fov_list))
    #draw a random integer for each fov
    random_trial_ints = np.random.randint(0, high=num_presentations, size=(len(spines_per_fov_list), draw_trials))
    rand_trials_for_spines = list(np.repeat(random_trial_ints, spines_per_fov_list))


    raise() #need to test that this works

    spine_indicies = np.arange(0, len(all_spine_activity_array['spines']), 1)
    numpyfied_activity = np.array(all_spine_activity_array)
    simulated_trail_traces = numpyfied_activity[spine_indicies, :, rand_trials_for_spines, :]

    simulated_trial_traces = xr.DataArray(simulated_trail_traces,
        coords={'spines': all_spine_activity_array['spines'], 'directions': all_spine_activity_array['directions'],'samples': all_spine_activity_array['samples']},
        dims=["spines", "directions", "samples"])
    return simulated_trial_traces #should be spines x directions x simulated_trials x samples



###############################################################################################################

def apply_weights(weights, traces):
    #you should not need to tile as long as one of them is the matching dimension. 
    # see https://stackoverflow.com/questions/32189190/numpy-array-multiplication-with-arrays-of-arbitrary-dimensions
    #tile_dims = (len(traces['directions']), len(traces['samples']), 1)
    weighted_traces = (weights*traces.T).T
    return weighted_traces




###############################################################################################################


###############################################################################################################
#Not particular to the model - mostly tuning curve and handling onset/offset/traces
###############################################################################################################








#def compare_tuning_curves_dot(means_1, means_2):
#    #first make both unit norm
#    means_1_unit_norm = means_1/np.dot(means_1, means_1)
#    means_2_unit_norm = means_2/np.dot(means_2, means_2)
#    return np.dot(means_1_unit_norm, means_2_unit_norm)


def compare_tuning_curves_anova(trial_amps_1_df, trial_amps_2_df):

    df = pd.concat([trial_amps_1_df, trial_amps_2_df])
    model = ols("amplitude ~ C(stim) + C(source) + C(source):C(stim)",data = df).fit()
        #C() indicates that the variable should be treated as categorical. They should also be strings at this point even though not strictly necessary
    anova_table = sm.stats.anova_lm(model, typ=2)
    p_value = anova_table.at[ 'C(source)' ,'PR(>F)']
    #the p value is the only thing we really need...
    return p_value

def convert_trial_amps_to_df(trial_amps, source=None):
    #each amplitude needs to be a row with columns: amplitude, model_name, stimulus
    num_presentations = trial_amps.shape[1]
    num_stims = trial_amps.shape[0]

    df_list = []
    for i in range(num_presentations):
        for j in range(num_stims):
            df_list.append({'source':source, 'stim': str(j), 'amplitude': trial_amps[j,i]})
    return pd.DataFrame(df_list)



#####################
#Neuron stuff

def compute_branch_order_and_dist(section, xyz_coords):
    order = 1
    dist = 0
    seg_fraction, seg_coords, _ = find_closest_segment(section, xyz_coords)
    dist_from_branch = section.L*seg_fraction
    dist+=dist_from_branch

    this_section = section
    for i in range(1000): #setting this to something reasonable for now, will avoid infinite loops
        parent  = this_section.parentseg().sec
        #print(parent.name())
        if 'soma' in parent.name():
            return order, dist, dist_from_branch
        order+=1
        dist += parent.L
        this_section = parent
    return order, dist, dist_from_branch


def get_branch_order_and_dist(h, spine_data, shifts_by_fov):
    stat_dict = {}
    for fov_num, fov in enumerate(spine_data['dend_cell'][2,:]):
        fov_name = hf.get_fov_name(spine_data, fov_num)
        print(fov_name)
        stat_dict[fov_name] = {}

        #current_input_dict[fov_num] = {}
        #ref = spine_data['dend_cell'][2,fov]
        fov_field_2 = spine_data[fov]#['DSI']
        spine_count = fov_field_2['trial_traces'].shape[-1]

        #should probably make this more elegant, I have copy and pasted this motif...
        spines_pixel_coords = hf.get_spines_pixel_coords(spine_data, fov_num)
        spines_global_coords = np.zeros((spines_pixel_coords.shape[1], 3))
        for i in range(spines_pixel_coords.shape[1]):
            spines_global_coords[i,:] = hf.get_spine_global_coords(h, spine_data, fov_num, i)

        manual_shift = shifts_by_fov[fov_num]
        #print(manual_shift)
        if (manual_shift == np.array([0,0,0])).all():
            print('using inferred shift')
            optimal_shift = list(estimate_offset(h, spine_data, fov_num, max_shift=10, iterations=2))
            #optimal_shift.append(0)
            shift = np.array(optimal_shift)
        else:
            print('using manual shift')
            shift = manual_shift



        test_spines_global_coords = shift_coords(spines_global_coords, shift)

        dist_list = []
        order_list = []
        for i in range(spine_count):
            nearest_section, sec_coords, min_sec_dist = hf.find_closest_section(h, test_spines_global_coords[i,:])
            #sec_fraction, seg_coords, min_seg_dist = find_closest_segment(nearest_section, test_spines_global_coords[i,:])
            order, dist, dist_from_branch = compute_branch_order_and_dist(nearest_section, test_spines_global_coords[i,:])
            print(f'order: {order}, dist: {dist}')
            dist_list.append(dist)
            order_list.append(order)
            #iclamp = h.IClamp(nearest_section(sec_fraction))
            #current_input_dict[fov_num][i] = iclamp
            #pp.pprint(nearest_section.psection()['point_processes'])
            #raise
        stat_dict[fov_name]['distance_to_soma'] = dist_list
        stat_dict[fov_name]['branch_order'] = order_list
    return stat_dict



def estimate_offset(h, spine_data, fov_num, max_shift = 10, iterations = 2, verbose = False):
    spines_pixel_coords = hf.get_spines_pixel_coords(spine_data, fov_num)

    coords_array, section_list = generate_section_mapped_3dpoints(h)


    spines_global_coords = np.zeros((spines_pixel_coords.shape[1], 3))
    for i in range(spines_pixel_coords.shape[1]):
        spines_global_coords[i,:] = hf.get_spine_global_coords(h, spine_data, fov_num, i)


    sections_to_search = []
    #I think its now running fast enough that we can just do it all again... but this is a bit cleaner and doesn't require passing back out the section names
    shift_idxs_to_test = [[-max_shift,-max_shift], [max_shift,-max_shift], [-max_shift,max_shift], [-max_shift,-max_shift], [0, 0]]
    for (x_shift, y_shift) in shift_idxs_to_test:
        shift = np.array([x_shift, y_shift, 0])
        test_spines_global_coords = spines_global_coords+np.tile(shift, (spines_global_coords.shape[0], 1))

        for i in range(test_spines_global_coords.shape[1]):
            spine_global_coords = test_spines_global_coords[i,:]
            nearest_section, sec_coords, min_dist  = find_closest_section_fast(coords_array, section_list, spine_global_coords)
            if not (nearest_section in sections_to_search):
                sections_to_search.append(nearest_section)

    print(f'Searching within these sections {sections_to_search}')


    # Did this to try to refine, make it use segment coordinates instead of section coordinates.
    full_coords_list = []
    full_section_list = []
    location_within_section_list = []
    for section in sections_to_search:
        coords_array, section_list, location_within_section = generate_segment_mapped_3dpoints(section)
        full_coords_list.extend(list(coords_array.T))
        full_section_list.extend(section_list)
        location_within_section_list.extend(location_within_section)
    all_segment_coords = np.array(full_coords_list)
    #print(f'####{all_segment_coords.shape}')
    #all_segment_coords = coords_array
    #full_section_list = section_list


    def find_optimal_shift(max_shift, iter_num, base_shift, all_segment_coords, full_section_list, verbose = verbose):
        if verbose:
            print(f'on iteration {iter_num}')
        cum_dist_at_shift = np.zeros((2*max_shift, 2*max_shift))

        # now implemented faster mechanism so don't have to regrab global spine coords

        for x_shift in range(-max_shift,max_shift):
            #print(x_shift)
            for y_shift in range(-max_shift,max_shift):
                this_x_shift =  x_shift/(10**iter_num)
                this_y_shift =  y_shift/(10**iter_num)
                manual_adjustment = np.array([this_x_shift+base_shift[0], this_y_shift+base_shift[1], 0])
                #print(manual_adjustment)

                cum_dist_at_shift[x_shift+max_shift, y_shift+max_shift] = find_distance_for_shift(all_segment_coords, full_section_list, spines_global_coords, manual_adjustment)
        #print(cum_dist_at_shift)

        #TODO if we wanted to make this better we could impose a cost here - don't want to move far distances for small gains
        #so scale the distance by the shift. But would require a lot of finnicky tweaking I think...

        optimal_shift = np.unravel_index(cum_dist_at_shift.argmin(), cum_dist_at_shift.shape)
        optimal_shift_relative = (np.array(optimal_shift) - np.array([max_shift, max_shift]))/(10**iter_num)+base_shift
        if verbose:
            print(f'Optimal shift is: {optimal_shift_relative}')
            print(f'Minimized cumulative distance is: {cum_dist_at_shift.min()}')
        return optimal_shift_relative

    base_shift = np.array([0,0])
    for i in range(iterations):
        base_shift = find_optimal_shift(max_shift, i, base_shift, all_segment_coords, full_section_list)
    optimal_shift = list(base_shift)
    optimal_shift.append(0)
    optimal_shift = np.array(optimal_shift)
    if verbose:
        check_min_dist = find_distance_for_shift(all_segment_coords, full_section_list, spines_global_coords, optimal_shift, verbose=verbose)
        print(f'Optimal shift is: {optimal_shift}')
        print(f'Minimized cumulative distance is: {check_min_dist}')
    return optimal_shift


def generate_section_mapped_3dpoints(h):
    coords_list = []
    section_list = []
    for section in h.allsec():
        points_list = section.psection()['morphology']['pts3d']
        coords_list.extend(points_list)
        this_sec_list = [section]*len(points_list)
        section_list.extend(this_sec_list)
    coords_array = np.array(coords_list)
    coords_array = coords_array[:, :3] #remove the radius

    return coords_array, section_list



def find_closest_section_fast(coords_array, section_list, xyz_coords):
    tiled_xyz_coords = np.tile(xyz_coords, (max(coords_array.shape),1))
    assert(tiled_xyz_coords.shape[1] == 3)
    diffs = coords_array - tiled_xyz_coords
    dists = np.linalg.norm(diffs, axis=1)

    nearest_section = section_list[dists.argmin()]
    nearest_coords = coords_array[dists.argmin(),:]
    return nearest_section, nearest_coords, dists.min()


def generate_segment_mapped_3dpoints(section):
    xCoord, yCoord, zCoord = hf.returnSegmentCoordinates(section)

    all_seg_coords = np.array([xCoord, yCoord, zCoord])
    section_list = [section]*len(xCoord)
    location_within_section = np.linspace(0,1,section.nseg)
    return all_seg_coords, section_list, location_within_section




def find_closest_segment(section, xyz_coords):
    xCoord, yCoord, zCoord = hf.returnSegmentCoordinates(section)

    all_seg_coords = np.array([xCoord, yCoord, zCoord])
    #print(seg_coords.shape)
    find_coords = np.tile(xyz_coords, (len(xCoord), 1)).T
    #print(find_coords.shape)
    dists = np.linalg.norm(all_seg_coords - find_coords, axis=0)
    #print(dists.shape)
    nearest_segment_i = dists.argmin()
    min_dist = dists.min()
    seg_coords = all_seg_coords[:,nearest_segment_i]

    #for i, (x,y,z) in enumerate(zip(xCoord, yCoord, zCoord)):
    #    dist = np.linalg.norm(xyz_coords - np.array([x,y,z])) #last point is radius
    #    try:
    #        if min_dist> dist:
    #            min_dist = dist
    #            nearest_segment_i = i
    #            seg_coords =  np.array([x,y,z])
    #    except UnboundLocalError as E:
    #        min_dist = dist
    #        nearest_segment_i = i
    #        seg_coords =  np.array([x,y,z])
    seg_fraction = nearest_segment_i/(section.nseg-1)
    return seg_fraction, seg_coords, min_dist,




def total_distance(coords_array, section_list, global_coords, verbose=False):
    total_distance = 0
    for i in range(global_coords.shape[0]):
        spine_global_coords = global_coords[i,:]
        nearest_section, sec_coords, min_dist  = find_closest_section_fast(coords_array, section_list, spine_global_coords)
        if verbose:
            print(f'spine num: {i}, spine coords: {spine_global_coords}')
            print(f'section: {nearest_section}, section_coords: {sec_coords}, min_sec_dist: {min_dist}')
        total_distance += min_dist
    #print(f'total_distance: {total_distance}')
    return total_distance


def find_distance_for_shift(coords_array, section_list, global_coords, apply_shift, verbose=False):
    if verbose:
        print(f'Calculating distance for all spines with a shift of {apply_shift}!!!!!')
    shifted_global_coords = global_coords - np.tile(apply_shift, (global_coords.shape[0], 1))
    return total_distance(coords_array, section_list, shifted_global_coords, verbose=verbose)


def find_distance_for_shift_slow(h, spine_data, fov_num, shift_coords, verbose=False):
    cumulative_distance=0
    spines_pixel_coords = hf.get_spines_pixel_coords(spine_data, fov_num)

    coords_array, section_list = generate_section_mapped_3dpoints(h)

    spines_global_coords = np.zeros((spines_pixel_coords.shape[1], 3))
    for i in range(spines_pixel_coords.shape[1]):
        spines_global_coords[i,:] = hf.get_spine_global_coords(h, spine_data, fov_num, i)

    cumulative_dist_fast_method = find_distance_for_shift(coords_array, section_list, spines_global_coords, shift_coords, verbose=verbose)
    if verbose:
        print(f'Cumulative distance from fast method: {cumulative_dist_fast_method}')


    for i in range(spines_pixel_coords.shape[1]):
        spines_global_coords[i,:] = hf.get_spine_global_coords(h, spine_data, fov_num, i, shift_coords)
    for i in range(spines_pixel_coords.shape[1]):
        nearest_section, sec_coords, min_dist = hf.find_closest_section(h, spines_global_coords[i,:])
        seg_fraction, seg_coords, min_dist = find_closest_segment(nearest_section, spines_global_coords[i,:])
        if verbose:
            print(f'Spine {i} is {min_dist} microns from section {nearest_section}[{seg_fraction}]')
        cumulative_distance+= min_dist
    return cumulative_distance



def shift_coords(coords_array, shift):
    return coords_array+np.tile(shift, (coords_array.shape[0], 1))


def create_current_sources(h, spine_data, shifts_by_fov):
    current_input_dict = {}
    for fov_num, fov in enumerate(spine_data['dend_cell'][2,:]):
        current_input_dict[fov_num] = {}
        #ref = spine_data['dend_cell'][2,fov]
        fov_field_2 = spine_data[fov]#['DSI']
        spine_count = fov_field_2['trial_traces'].shape[-1]

        #should probably make this more elegant, I have copy and pasted this motif...
        spines_pixel_coords = hf.get_spines_pixel_coords(spine_data, fov_num)
        spines_global_coords = np.zeros((spines_pixel_coords.shape[1], 3))
        for i in range(spines_pixel_coords.shape[1]):
            spines_global_coords[i,:] = hf.get_spine_global_coords(h, spine_data, fov_num, i)

        manual_shift = shifts_by_fov[fov_num]
        #print(manual_shift)
        if (manual_shift == np.array([0,0,0])).all():
            print('using inferred shift')
            optimal_shift = list(estimate_offset(h, spine_data, fov_num, max_shift=10, iterations=2))
            #optimal_shift.append(0)
            shift = np.array(optimal_shift)
        else:
            print('using manual shift')
            shift = manual_shift

        test_spines_global_coords = shift_coords(spines_global_coords, shift)

        for i in range(spine_count):
            nearest_section, sec_coords, min_sec_dist = hf.find_closest_section(h, test_spines_global_coords[i,:])
            sec_fraction, seg_coords, min_seg_dist = find_closest_segment(nearest_section, test_spines_global_coords[i,:])

            iclamp = h.IClamp(nearest_section(sec_fraction))
            current_input_dict[fov_num][i] = iclamp
            #pp.pprint(nearest_section.psection()['point_processes'])
            #raise
    return current_input_dict

