from src import data_io as io
from src import plotting as plot
from src import computation as comp
from src import helper_functions as hf
from src import config as cfg
from src import weight_generators
from src import inclusion_mask_generators
from src import integration_functions
from src import somatic_functions
from src import evaluation_functions
from src import trial_selection

import matplotlib.colors as cm
from datetime import datetime


import xarray as xr
#import pandas as pd
import numpy as np
import ipdb
from matplotlib import pyplot as plt
#import seaborn as sns
import pandas as pd

from PIL import Image #this needs to be after matplotlib??
from scipy import stats  
from collections import defaultdict
import os
import traceback
import logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)


def init_globals():
    globals = defaultdict(dict)
    globals['failed_list']= []
    globals['data_dict_list'] = [] #Best practice is to create a list of dictionaries and then convert to df
    #see https://stackoverflow.com/questions/28056171/how-to-build-and-fill-pandas-dataframe-from-for-loop
    #and https://stackoverflow.com/questions/57000903/what-is-the-fastest-and-most-efficient-way-to-append-rows-to-a-dataframe
    return globals

def main_loop():
    globals = init_globals()
    
    if os.path.isdir(cfg.data_path) and cfg.walk_dirs:
        print('Walking to multiple datasets')
        for current_data_dir, dirs, files in os.walk(cfg.data_path, topdown=False):
            for filename in files:
                filepath = os.path.join(current_data_dir, filename)
                unprocessed = all([not(skip_keyword in filename) for skip_keyword in cfg.skip_keywords])
                soma_df = ('ASC' in filename) or ('BM' in filename)
                if soma_df and (cfg.re_run or unprocessed):
                    main(filepath, globals)
    else:
        print('Running on single dataset')
        main(cfg.data_path, globals)
    
    print("###########################################################")
    print("###########################################################")
    print("Finished running all models for all cells, saving aggregate data")
    print("###########################################################")
    #io.save_summary_plots(globals)
    
    df = pd.DataFrame(globals['data_dict_list'])
    date_str = datetime.today().strftime('%Y%m%d%H%M%S')
    name_str = f"single_neuron_simulation_scores_{cfg.simulated_trials_per_stim}_{date_str}"
    io.save_csv(df, name_keywords=name_str)

    io.save_named_iterable_to_json(failed_dirs_list=globals['failed_list'])
    io.save_named_iterable_to_json(failed_dirs_errors=globals['errors'])


def tuning_curve_df_helper(model_keyword, means):
    result_list = [model_keyword]
    column_list = ['model_keyword (V), stimulus (->)']
    for i, xarray_mean in enumerate(means):

        #can add stim descriptions here
        #print(xarray_mean.Coordinates)
        try:
            column_list.append(float(xarray_mean.directions))#need to make this label more general
            result_list.append(float(xarray_mean.values))
        except AttributeError as E:
            column_list.append(i)
            result_list.append(xarray_mean)
    return column_list, result_list


def main(current_data_path, globals = None):
    if globals is None:
        globals = init_globals()
    print("###########################################################")
    print("###########################################################")
    print("Current data dir: " + str(current_data_path))
    print("###########################################################")
    print("###########################################################")
    print("")
    print("")
    try:
        #load soma data
        current_data_dir, filename = os.path.split(current_data_path)
        experiment_id, extension = os.path.splitext(filename)
        soma_path = current_data_path
        soma_data = io.loadmat(soma_path)
        model_spec_df = cfg.run_model_specs

        if not(hf.include_soma(soma_data)):
            print("Excluding: " + str(current_data_path))
        else:
            print("Attempting to run adult spine models for data dir: " + str(current_data_path))
            responsive_status = hf.get_responsive_status(soma_data)
            print(f'responsive status: {hf.get_responsive_status(soma_data)}')
        

    ############################################
    #Get everything ready to run all the models (lists for each parameter)
    ############################################

            print('Loading spine data and filtering traces')
            ####load spine data
            spines_path = cfg.get_spines_path(soma_path)
            spine_data = io.loadmat(spines_path)


            #the only two that we need to do the legwork for beforehand are the included trials and the weights
            #the rest are all fast, static functions that can be called when we need them. 
            #### Get diferenct trace matricies for each exclusion criteria
            included_trials_list = {
                'all_trials': trial_selection.all_trials,
                'baps_trials_only': trial_selection.baps_trials_only, 
                'no_bap_trials': trial_selection.no_bap_trials,
                }
            inclusion_type_spine_traces_dict = {}
            spine_counts_dict = {}
            for inclusion_criteria, mask_func in included_trials_list.items():
                inclusion_type_spine_traces_dict[inclusion_criteria], spine_counts_dict[inclusion_criteria] = trial_selection.compile_spine_traces(spine_data, mask_func=mask_func)
###Could save these as intermediate? does this step take time?

            def generate_weight_matricies(model_spec_df, spine_data):
            #### Get weight matricies for each model we want to test
                print("###########################################################")
                print('')
                print('Generating weight matricies for each model')

                weight_matricies_dict = {}
                all_weightings = list(model_spec_df['weighted_by'])+list(model_spec_df['inclusion_based_on'])
                
                #commenting these out to do them in order for now
                #unique_weightings = all_weightings
                unique_weightings = set(all_weightings)
                unique_weightings = {x for x in unique_weightings if x==x} #remove nans
                unique_weightings.remove('_')
                for weight_name in unique_weightings:
                    print(f'Generating weight matrix for {weight_name}')
                    func_name = f'weights_from_{weight_name}'
                    #print('func_name', func_name)
                    try:
                        weight_func = getattr(weight_generators, func_name)
                    except AttributeError as E:
                        #If we are not rescaling the values for the weighting (applying a weight function to the retrieved values)
                        #Then we can pull the values directly using the get_fovs_attr functions in hf and the default weight_identity
                        #Probably we should specify this in the model_specifications.csv, or even precompute and save in a seperate module. 
                        #we should definitely save the weightiings we used.
                        func_name = f'get_fovs_{weight_name}'
                        param_func = getattr(hf, func_name)
                        
                        def weight_func(spine_data):
                            return weight_generators.get_weight_array(spine_data, param_func)
                        #weight_func = weight_generators.get_weight_array(spine_data, param_func)
                    weights = weight_func(spine_data)
                    #print('weights', weights)
                    weight_matricies_dict[weight_name] = weights
                    print(weight_name)
                    print(weights)
                return weight_matricies_dict

            weight_matricies_dict = generate_weight_matricies(model_spec_df, spine_data)
            print(weight_matricies_dict.keys())
            score_function_for_fitting = getattr(evaluation_functions, cfg.score_for_fitting)

    ############################################
    #First we generate the input to the soma for each model as specified in cfg.model_specifications
    ############################################

            print("###########################################################")
            print('')
            print('Begining to run the models')
            #a dictionary to put all the results in
            model_outputs = {}
            model_details_dict = {}

            #intiialize the dictionaries and model traces outputs
            for i, (index, model_spec) in enumerate(cfg.run_model_specs.iterrows()):

                full_model_name, nickname = hf.get_model_names(model_spec, index)
                model_details_dict[full_model_name] = {
                    'experiment_id': experiment_id,
                    'soma_responsive_status': responsive_status,
                    'model_nickname': nickname,
                    'full_model_name': full_model_name,
                    'trials_included': model_spec.trials_included,
                    'weighted_by': model_spec.weighted_by,
                    'spines_included': model_spec.spines_included,
                    'integration_function': model_spec.integration,
                    'somatic_function': model_spec.somatic_function,
                }

                spine_traces = inclusion_type_spine_traces_dict[model_spec.trials_included]
                model_outputs[full_model_name] = comp.init_traces_xarray(spine_traces, cfg.simulated_trials_per_stim)

            included_trials_list = [model_spec.trials_included for index, model_spec in cfg.run_model_specs.iterrows()]
            for included_trial_type, spine_traces in inclusion_type_spine_traces_dict.items():
                spines_per_fov_list = spine_counts_dict[included_trial_type]


                if included_trial_type in included_trials_list:
                        print(f"Running models on included trials: {included_trial_type}")
                        #Don't even loop through the simulated trials if none of the models need it

                        for i in range(cfg.simulated_trials_per_stim):
                            simulated_trial_traces = comp.sample_trial_from_fov(spine_traces, spines_per_fov_list)
                            #print('input', simulated_trial_traces)

                            #TODO eventually should numpify this ^^^ loop.... the whole matrix of traces could be shared across all models of the same trial type
                            #instead of just sharing 1 trial at a time
                            #just added an s in the function name
                            #simulated_trial_traces = comp.sample_trials_from_fov(spine_traces, spines_per_fov_list)
                            #this will be tougher if not everything is matrix multiplicatin... have to make those functions copatible too. 
                            

                            #preallocate disk space for memmap
                            #may actually not even want to do this. If we leave it in memory - many of them are redrawn, pointing to the same original trials
                            #so won't necessarily take up any extra space? have to test and see, numpy may actually make the copy for speed
                            #put all the simulated trial traces in there (even this could be numpified probably...)

                            #then loop through the models and apply the functions to the whole matrix all at once. 
                            for k, model_spec in cfg.run_model_specs.iterrows():
                                #print(f'model_k: {k}')
                                #print(model_spec)
                            #for k, model_spec in enumerate(cfg.models_list):
                                if model_spec.trials_included == included_trial_type:
                                #only run the models that require this included trial

                                    spines_per_fov_list = spine_counts_dict[model_spec.trials_included]

                                    ############################
                                    weights = weight_matricies_dict[model_spec['weighted_by']]
                                    #print('weights')
                                    #print(weights)

                                    ############################
                                    mask_func = getattr(inclusion_mask_generators, model_spec['spines_included'])
                                    if (str(model_spec['inclusion_based_on']) == 'nan') or (str(model_spec['inclusion_based_on']) == '_'):
                                        #intentionally using if-else instead of try, because we don't want other missing keys to pass through with the all spines mask
                                        param_array_determining_inclusion = weight_matricies_dict['democratic_weights']
                                    else:
                                        param_array_determining_inclusion = weight_matricies_dict[model_spec['inclusion_based_on']]

                                    #print('param_array_determining_inclusion')
                                    #print(param_array_determining_inclusion)
                                    mask = mask_func(param_array_determining_inclusion)
                                    #print('mask')
                                    #print(mask)

                                    ############################
                                    masked_weights = weights*mask
                                    #print('masked weights')
                                    #print(masked_weights)
                                    weighted_simulated_trial_traces = comp.apply_weights(masked_weights, simulated_trial_traces)
                                    #print('weights', weights)
                                    #print('mask', mask)
                                    #### ^^ up to here could be numpified easily. 

                                    #############################
                                    integration_function_name = f"{model_spec['integration']}_integration"
                                    integration_function = getattr(integration_functions, integration_function_name)
                                    #apply integration model here
                                    simulated_input_to_soma = integration_function(weighted_simulated_trial_traces, spines_per_fov_list)

                                    #add the new trace to the initialized array
                                    full_model_name, nickname = hf.get_model_names(model_spec, k)
                                    #print('output', simulated_input_to_soma)
                                    model_outputs[full_model_name][:,i,:] = simulated_input_to_soma

                                    #except KeyError as E: #Don't need anymore, we initialize earlier
                                    #    model_outputs[full_model_name] = comp.init_traces_xarray(spine_traces, cfg.simulated_trials_per_stim)
                                    #    model_outputs[full_model_name][:,i,:] = masked_simulated_trial_traces

    ############################################
    #Second we apply and fit the somatic output function for each model as specified in model_specifications.csv this could actually come later after saving. Making it more an entirely different module. its only the integration function that needs to be applied so we don't keep a bunch of uneceesary traces. 
    ############################################
            
            soma_traces = hf.get_soma_traces(soma_data)
            match_list = []
            match_counter = 0
            for k, model_spec in cfg.run_model_specs.iterrows():
                full_model_name, nickname = hf.get_model_names(model_spec, k)
                simulated_input_to_soma = model_outputs[full_model_name] 
                somatic_function = getattr(somatic_functions, f"somatic_{model_spec['somatic_function']}")
                simulated_output_of_soma, score, fit_params = somatic_function(simulated_input_to_soma, soma_traces, score_function=score_function_for_fitting)
                #print('output2', simulated_input_to_soma)
                model_outputs[full_model_name] = simulated_output_of_soma #overwriting here. We could have an option to save it beforehand - more effciitne to test all different somatic functions on same inuput, but only apply the weightings once instead of applying them for each different somatic function

                #for debugging
                match_list_2 = []
                for j, model_spec_2 in cfg.run_model_specs.iterrows():
                    full_model_name_2, nickname = hf.get_model_names(model_spec_2, j)
                    simulated_input_to_soma_2 = model_outputs[full_model_name_2] 
                    match = int(np.all(simulated_input_to_soma_2 == simulated_input_to_soma))
                    match_list_2.append(match) 
                    if match:
                        match_counter += 1

                match_list.append(match_list_2)
            #print(f'num_models: {k}, {j}')
            #print(f'number_matches: {match_counter}')
            #print(f'number_matches: {np.sum(np.array(match_list))}')
            #print(f'number_comparisons: {np.array(match_list).shape}')
            #print('MAtch_list')
            #print(match_list)
            #plt.imshow(np.array(match_list))
            #raise()
                    


    ############################################
    #Compare the model outputs with the somatic traces - this could actually come later after saving. Making it more an entirely different module. its only the integration function that needs to be applied so we don't keep a bunch of uneceesary traces. 
    ############################################

            temporal_likelihoods_list = []
            print("###########################################################")
            print('')
            print(f"Computing similarity metrics for each model")
            for full_model_name, model_traces in model_outputs.items():

                #print('model_traces', model_traces)
                #Compute various similarity metrics and add to the details dictionary for that model
                ############
                scores_dictionary, err_str, temporal_likelihood = evaluation_functions.compute_model_scores(soma_data, model_traces)
                model_details_dict[full_model_name].update(scores_dictionary)
                #import pdb; pdb.set_trace() 
                df_row = [full_model_name]
                df_row.extend(list(temporal_likelihood))
                temporal_likelihoods_list.append(df_row)

            #save these
            #consider adding columns to this one - need to get the lables of samples from soma_traces
            column_list = ['model_name'].extend(list(model_traces.coords['samples'].values))
            likelihood_df = pd.DataFrame(temporal_likelihoods_list, columns = column_list)
            name_str = f'{experiment_id}_temporal_loglikelihood.csv'
            tuning_curve_csv_path = os.path.join(cfg.collect_summary_at_path, experiment_id, name_str)
            dirpath = os.path.split(tuning_curve_csv_path)[0]   
            if not (os.path.isdir(dirpath)):
                os.mkdir(dirpath)
            likelihood_df.to_csv(tuning_curve_csv_path)
            #io.save_csv(df, name_keywords=name_str)

    ############################################
    #Save the traces and tuning curves
    ############################################
            print("###########################################################")
            print('')
            print(f"Saving the tuning curves for {experiment_id}")


            #initialize a list of lists where each list of our dataframe is a row
            simulation_means_list = []

            #first add the soma
            soma_traces = hf.get_soma_traces(soma_data)
            soma_tuning_curve_normalized, soma_max_amplitude = evaluation_functions.compute_normalized_tuning_curves(soma_traces)
            column_list, result_list = tuning_curve_df_helper('soma', soma_tuning_curve_normalized)
            simulation_means_list.append(result_list)

            for k, model_spec in cfg.run_model_specs.iterrows():
                full_model_name, nickname = hf.get_model_names(model_spec, k)
                model_traces = model_outputs[full_model_name]
                #nickname = model_details_dict[full_model_name]['model_nickname']
                #then add a tuning curve from each simulation
                model_tuning_curve_normalized, model_max_amplitude = evaluation_functions.compute_normalized_tuning_curves(model_traces)
                column_list, result_list = tuning_curve_df_helper(full_model_name, model_tuning_curve_normalized)
                simulation_means_list.append(result_list)

                #save the traces image - could make this configureable to each model spec if we wanted to only save some of them...
                #also should be in a seperate module. Not combined with putting the results in the larger DF. traces should be treated seperately from images too. 
                if cfg.save_traces:
                    #print('model_traces', model_traces)
                    
                    
                    #name = f'{experiment_id}_{full_model_name}_{cfg.simulated_trials_per_stim}'
                    plot.save_trace_image(model_traces, experiment_id=experiment_id, full_model_name=full_model_name, prefix=f'simulated_traces')
                    subdir = os.path.join(experiment_id, "real_soma_2P_measurements")
                    plot.save_trace_image(soma_traces, experiment_id=experiment_id, prefix = 'measured_traces')

            tuning_curve_df = pd.DataFrame(simulation_means_list, columns = column_list)
            name_str = f'{experiment_id}_simulation_mean_stim_response.csv'
            tuning_curve_csv_path = os.path.join(cfg.collect_summary_at_path, experiment_id, name_str)
            tuning_curve_df.to_csv(tuning_curve_csv_path)
            #io.save_csv(df, name_keywords=name_str)

 

    ############################################
    #Finally we run the shuffles if specified in cfg.num_shuffles, alternatively we just load the results from previously computed shuffles
    ############################################

            shuffle_scores_dict = {}
            if cfg.run_shuffles:
                for k, model_spec in enumerate(cfg.models_list):
                    if model_spec.run_shuffles:

                        shuffle_scores_list = [] #this list will become a dataframe of the scores computed for each shuffle and be saved
                        full_model_name, nickname = hf.get_model_names(model_spec, k)
                        print('Running shuffles for {full_model_name}')
                        ################################
                        included_trial_types = model_spec.trials_included
                        spine_traces = inclusion_type_spine_traces_dict[included_trial_types]
                        spines_per_fov_list = spine_counts_dict[included_trial_types]
                        for j in range(cfg.num_shuffles):
                            if j % 100 == 0:
                                print(f'shuffle {j}')

                            #shuffle spine traces here
                            shuffled_spine_traces = comp.shuffle_spine_traces(spine_traces)
                            #need it to be an xarray for some functions I think
                            shuffled_spine_traces = xr.DataArray(shuffled_spine_traces,
                                coords={'spines': spine_traces['spines'], 'directions': spine_traces['directions'],'samples': spine_traces['samples']},
                                dims=['spines', "directions", "presentations", "samples"]
                                )
                            #this is spines x directions x presentations x samples
                            
                            this_shuffle_output_traces = comp.init_traces_xarray(spine_traces, cfg.simulated_trials_per_stim)

                            for i in range(cfg.simulated_trials_per_stim):
                                simulated_trial_traces = comp.sample_trial_from_fov(shuffled_spine_traces, spines_per_fov_list)
                                #TODO eventually should numpify this ^^^ loop....

                                ############################
                                weights = weight_matricies_dict[model_spec.weighted_by]
                                #multiply by weights here
                                weighted_simulated_trial_traces = comp.apply_weights(weights, simulated_trial_traces)

                                ############################
                                mask = getattr(inclusion_mask_generators, model_spec.spines_included)
                                #multiply by weights here
                                masked_simulated_trial_traces = comp.apply_weights(mask, weighted_simulated_trial_traces)

                                #############################
                                integration_function = integration_functions[model_spec.integration]
                                #apply integration model here
                                simulated_input_to_soma = integration_function(weighted_simulated_trial_traces)

                                ##############################
                                somatic_function  = somatic_functions[model_spec.somatic_function]
                                #apply somatic nonlinearity (if using) here
                                simulated_output_of_soma = somatic_function(simulated_input_to_soma, soma_traces, score_function=None)

                                this_shuffle_output_traces[:,i,:] = simulated_output_of_soma


                            #Compute various similarity metrics
                            ############
                            shuffle_scores_dict, err_str , tl = comp.compute_model_scores(soma_data, this_shuffle_output_traces)
                            if err_str:
                                globals['errors'][current_data_dir] = f'{err_str}{full_model_name}'
                            #now we need to put them in a dataframe so all scores can be saved.
                            shuffle_scores_list.append(shuffle_scores_dict)

                        #save the dataframe if we need to grab the distribution again... the shuffling really belongs in a seperate script with a seperate df. 
                        #then here we should intentionally reload it to compute the p value
                        df = pd.DataFrame(shuffle_scores_list)
                        io.save_csv(df, name_keywords=f'shuffle_scores_{cfg.num_shuffles}_{experiment_id}_{full_model_name}')
     

    ############################################
    #load in the shuffles and compute the p value from them, add to be saved in the df
    ############################################
            if False:
            
                for k, model_spec in enumerate(cfg.models_list):
                    if model_spec.run_shuffles:
                        full_model_name, nickname = hf.get_model_names(model_spec, k)
                        
                        df = io.from_csv(df, name_keywords=f'shuffle_scores_{experiment_id}_{full_model_name}')
                        io.save_csv(df, name_keywords=f'shuffle_scores_{cfg.num_shuffles}_{experiment_id}_{full_model_name}')
                        #df = io.from_csv(df, name_keywords=f'shuffle_scores_{cfg.num_shuffles}_{experiment_id}_{full_model_name}')

                        unshuffled_traces = model_outputs[full_model_name]
                        unshuffled_scores_dict, err_str, tl = comp.compute_model_scores(soma_data, unshuffled_traces)
                        for metric_name, unshuffled_model_score in unshuffled_scores_dict.items():
                            try:
                                shuffle_array = df.loc[:, metric_name]
                                num_lower = (shuffle_array > unshuffled_model_score).sum()
                                this_score_pvalue = num_lower/len(shuffle_array)
                                key = f'{metric_name}_p_value'
                                model_details_dict[full_model_name][key] = this_score_pvalue
                            except KeyError as E:
                                #shuffle not present for this metric
                                pass
                
                    #This doesn't get saved within this function... should it? I've actually never really had it crash where I would want to resurrect it. the whole computation is fast enough

    ############################################
    #Finally #2, add all the details for each model to the global dataframe 
    ############################################
            for full_model_name, model_details_dict in model_details_dict.items():
                globals['data_dict_list'].append(model_details_dict)


    except Exception as E:
        globals['failed_list'].append(current_data_dir)
        err_str = f"There was an error processing data file: {current_data_path}"
        logger.error(err_str)
        logger.warning(traceback.format_exc())
        globals['errors'][current_data_dir] = traceback.format_exc()
        raise(E)
    
    
if __name__ == "__main__":
    main_loop()
    
    