
##################
#DEPRECATED

#crete a list of 100 numbers


a = 1
b = 2
c = a+b+10
print(c)





#def run_model(spine_data, spine_activity_array, spines_per_fov_list, weight_function):
#    model_traces = compute_model_output_from_random_sampled_fovs(spine_data,
#                                                                      spine_activity_array,
#                                                                      spines_per_fov_list,
#                                                                      simulated_trials_per_stim=cfg.simulated_trials_per_stim)
#    model_tuning_curve_normalized, model_max_amplitude = compute_normalized_tuning_curves(model_traces)
#    return model_traces, model_tuning_curve_normalized, model_max_amplitude



#Now we need a function to slice this array meaningfully
#basically we need a different random integer for each FOV
#mutliply that out to be an array of the right length
#def compute_model_output_from_random_sampled_fovs(spine_data,
#                                                all_spine_activity_array,
#                                                spines_per_fov_list,
#                                                simulated_trials_per_stim = 10,
#                                                weight_function=democratic_weights,
#                                                integration_function = linear_integration,
#                                                somatic_function = somatic_identity,
#                                                ):#

    ##Need to break this into seperate modules. Bootsrapping will be a bit tricky? What if different nonlinearities are better for different runs?
    #I guess you run them all and then sort it out after.
    #we dont want to have to save all the model outputs, that will get kinda big
    #so its ok to bootstrap this function after we have extracted the rele
    #except a list for each - full prameters sets of weight, integration and somatic pairs. Also should pass in flags for unresponsive and exclude baps here


#    num_directions = len(all_spine_activity_array['directions'])
#    num_samples= len(all_spine_activity_array['samples'])
#    model_output = xr.DataArray(np.zeros((num_directions, simulated_trials_per_stim, num_samples)),
#                                coords={'directions': all_spine_activity_array['directions'],'samples': all_spine_activity_array['samples']},
#                                dims=["directions", "presentations", "samples"])#
#
#    for i in range(simulated_trials_per_stim):
#        simulated_trial_traces = sample_trial_from_fov(all_spine_activity_array, spines_per_fov_list)#
#
#        #multiply by weights here
#        weighted_simulated_trial_traces = weight_function(spine_data, simulated_trial_traces)##
#
#        #apply integration model here
#        simulated_input_to_soma = integration_function(weighted_simulated_trial_traces)#
#
#        #apply somatic nonlinearity here
#        simulated_output_of_soma = somatic_function(simulated_input_to_soma)#
#
#        model_output[:,i,:] = simulated_output_of_soma
##
#
#    return model_output #should be directions x simulated_trials x samples (like the soma)


def test_MC_Pitts_model( soma_data, soma_act_thresh):
    #For now just minumize the pairwise difference between the two?
    pass


#def plot_MC_Pitts_model():
#    #sort by activity sorting within each direction (based off counts)#

#
    #sort by same directional sorting as somas

    #plot both
#    fig, axs = plt.subplots(1,2)
#    axs[0].imshow(flatten_for_image(bool_soma_activity))
    #axs[1].imshow(

#def run_MC_Pitts_model(spine_data, spine_act_thresh, soma_count_thresh):#

    #loop through each spine
    #TODO could we make this a generator?
    count_active_spines = np.zeros(spine_traces.shape)
    for fov in spines['dend_cell'][2,:]:
        ref = spines['dend_cell'][2,fov]
        fov_field_2 = spines[ref]#['DSI']
        spine_count = fov_field_2['trial_traces'].shape[-1]
        for i in range(spine_count):
            this_spine_traces = np.array(fov_field_2['trial_traces'][:,:,0,:,i].swapaxes(0,-1))

            #select the time period means of interest
            this_spine_sub_traces = select_timesteps(soma_traces)
            this_spine_period_means = get_period_means(this_spine_sub_traces)

            #apply spine_act_threshold on each
            bool_input = int(this_spine_period_means>spine_act_thresh)

            #add to sum (bool will make it a count
            count_active_spines += bool_input

    #done looping
    #apply soma_count_thresh to get to bool
    simulated_soma_response =  int(count_active_spines>soma_count_thresh)
    return simulated_soma_response





def get_summed_spine_trace_depracated(spine_data):
    #Get all the traces from all the spines
    #TODO could we make this a generator?

    ex_spine_trace = hf.get_example_traces(spine_data)
    summed_spine_traces = np.zeros(ex_spine_trace.shape)

    for fov in spine_data['dend_cell'][2,:]:
        #ref = spine_data['dend_cell'][2,fov]
        fov_field_2 = spine_data[fov]#['DSI']
        spine_count = fov_field_2['trial_traces'].shape[-1]
        for i in range(spine_count):
            this_spine_traces = hf.get_traces(spine_data, fov=fov, spine_index=i)
            summed_spine_traces += this_spine_traces

    return summed_spine_traces


def get_summed_trial_sampled_spine_trace_depracated(spine_data):
    sampling_mat = trial_sampling(spine_data)
    #dimensions = stims x simulated_trials x FOVs.

    ex_spine_trace = hf.get_example_traces(spine_data)
    stims = ex_spine_trace.shape[0]
    samples = ex_spine_trace.shape[-1]
    simulated_trials = sampling_mat.shape[1]
    num_fovs = sampling_mat.shape[-1]
    assert stims == sampling_mat.shape[0]

    summed_traces = np.zeros((stims, simulated_trials, samples))
    for stim_num in range(stims):
        for trail_num in range(simulated_trials):
            for fov_num in range(num_fovs):
                this_fov = hf.get_fov(spine_data, fov=fov_num)
                spine_count = this_fov['trial_traces'].shape[-1]
                for i in range(spine_count):
                    #Might be faster not to go all the way back to the file here...
                    this_spine_traces = hf.get_traces(spine_data, fov=fov_num, spine_index=i)
                    og_trial_num = sampling_mat[stim_num, trail_num, fov_num]
                    summed_traces[stim_num, trail_num, :] += this_spine_traces[stim_num, og_trial_num, :]
    return summed_traces




def trial_sampling_deprecated(spine_data):
    simulated_trials_per_stim = 10

    ex_spine_trace = hf.get_example_traces(spine_data)
    stims = ex_spine_trace.shape[0]
    stim_repeats = ex_spine_trace.shape[1]
    fovs = spine_data['dend_cell'][2,:].shape[0]

    sampling_mat = np.random.randint(0,stim_repeats,(stims, simulated_trials_per_stim, fovs)) #dimensions = stims x simulated_trials x FOVs. each index is the og_trial - direction as assumed the same order
    return sampling_mat
