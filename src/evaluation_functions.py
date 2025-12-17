"""description of this module, interface.py(interface).
please update this docstring as you develop.

interface.py should contain functions that corresponde to each keyword in the "weighted_by" columns of the model_specification.csv. 
Each keyword needs to map to a function here that will return the right vector of n weights for n spines
"""
from src import helper_functions as hf
import numpy as np
from sklearn.neighbors import KernelDensity
from scipy import stats 

from src import somatic_functions as soma

#score_functions = {'stim_on_llikelihood': comp.stim_on_neg_llikelihood,
#                    'similarity_score': comp.dotprod_based_similarity_score_neg,
#                  }


#########
#Default functions must be defined first so they have something to point to
def linear_normalization(array_in, norm_by=None):
    if norm_by is None:
        norm_by = array_in
    zeroed_array = array_in-np.min(norm_by)
    zeroed_norm_by = norm_by-np.min(norm_by)
    this_max = np.max(zeroed_norm_by)
    if this_max > .1:
        return zeroed_array/this_max
    else:
        return zeroed_array
    

def compute_trial_means(traces):
    on_period = hf.get_stim_on_traces(traces)
    #This is for the anova - output should be very similar to kyles trial amps
    #should be shape: stims x presentations as opposed to the tuning curve which is just stims x 1
    #on_period = on_period.reshape(on_period.shape[0], on_period.shape[1]*on_period.shape[2])
    try:
        trial_means = on_period.mean(dim='samples')
    except TypeError as E:
        trial_means = on_period.mean(axis=1)
    return trial_means


def compute_normalized_trial_means(traces):
    trial_means = compute_trial_means(traces)
    normalized_trial_means = linear_normalization(trial_means)
    return normalized_trial_means

def compute_mean_tuning(traces):
    trial_means = compute_trial_means(traces)
    try:
        stim_means = trial_means.mean(dim='presentations')
    except TypeError as E:
        stim_means = trial_means.mean(axis=1)
    return stim_means

def compute_median_tuning(traces):
    trial_means = compute_trial_means(traces)
    try:
        stim_medians = trial_means.median(dim='presentations')
    except (TypeError, AttributeError) as E:
        stim_medians = np.median(trial_means, axis=1)
    return stim_medians

def normalize_traces_by_tuning_curve(traces):
    tuning_curve = compute_tuning_curves(traces)
    return linear_normalization(traces, norm_by = tuning_curve)





def compute_and_compare_tuning_curves(traces_1, traces_2):
    means_1 = compute_tuning_curves(traces_1)
    means_2 = compute_tuning_curves(traces_2)
    return compare_tuning_curves(means_1, means_2)

def neg_compare_tuning_curves(means_1, means_2):
    return -1*compare_tuning_curves(means_1, means_2)

def compare_tuning_curves(means_1, means_2):
    #first make both unit norm
    means_1_unit_norm = make_unit_norm(means_1)
    means_2_unit_norm = make_unit_norm(means_2)
    return np.dot(means_1_unit_norm, means_2_unit_norm)

def make_unit_norm(array_in):
  return array_in/np.linalg.norm(array_in)


def compute_tuning_curves(traces):
    return compute_median_tuning(traces)

def compute_normalized_tuning_curves(traces):
    tuning_curve = compute_tuning_curves(traces)
    normalized_tuning_curve = linear_normalization(tuning_curve)
    max_amp = np.max(tuning_curve)
    return normalized_tuning_curve, max_amp


#############

def compute_model_scores(soma_data, model_traces):

    #Correlation
    soma_traces = hf.get_soma_traces(soma_data)
    soma_traces_normalized = normalize_traces_by_tuning_curve(soma_traces)
    #soma_means_normalized, soma_max_amplitude = compute_normalized_tuning_curves(soma_traces)
    soma_tunning_curve = hf.get_precomputed_tuning_curve(soma_data)
    soma_means_normalized = linear_normalization(soma_tunning_curve)
    model_tuning_curve_normalized, model_max_amplitude = compute_normalized_tuning_curves(model_traces)

    err_str=None
    try:
        model_corr_to_soma = stats.pearsonr(soma_means_normalized, model_tuning_curve_normalized)
    except ValueError as E:
        err_str = f'One or more of the model tuning curves seems to have been all NaNs'
        model_tuning_curve_normalized = np.nan_to_num(model_tuning_curve_normalized, nan=0.0, posinf=0.0, neginf=0.0)
        model_corr_to_soma = stats.pearsonr(soma_means_normalized, model_tuning_curve_normalized)
    #Dot product based similarity score

    relu_applied_output, best_score, fit_params = soma.somatic_2_param_relu(model_tuning_curve_normalized, soma_means_normalized, neg_compare_tuning_curves)
    #print(model_tuning_curve_normalized)
    #print(soma_means_normalized)
    #print(relu_applied_output)
    #print(best_score)
    #print(fit_params)
    
    model_similarity_score = compare_tuning_curves(soma_means_normalized, relu_applied_output)
    #compare_tuning_curves(soma_means_normalized, model_tuning_curve_normalized)
    


    mss2 = dotprod_based_similarity_score(model_traces, soma_traces_normalized)

    mean_log_likelihood, temporal_likelihood = model_likelihood(model_traces, soma_traces_normalized)
    #temporal_likelihood = [0]

    stim_on_ll = stim_on_llikelihood(model_traces, soma_traces_normalized)

    scores_dictionary = {
        'model_correlation_to_soma_r': model_corr_to_soma[0],
        'model_correlation_to_soma_p': model_corr_to_soma[1],
        'model_soma_similarity_score': model_similarity_score,
        'model_soma_sim_score_recomputed_tuningcurve':mss2,
        'mean_amplitude_stim_on_ll': stim_on_ll,
        'mean_log_likelihood': mean_log_likelihood
        }
    return scores_dictionary, err_str, temporal_likelihood #needs to be saved differently


def dotprod_based_similarity_score_neg(model_traces, soma_traces):
    return -1*dotprod_based_similarity_score(model_traces, soma_traces)

def dotprod_based_similarity_score(model_traces, soma_traces):
    model_tuning_curve_normalized, model_max_amplitude = compute_normalized_tuning_curves(model_traces)
    soma_tuning_curve_normalized, soma_max_amplitude = compute_normalized_tuning_curves(soma_traces)
    return compare_tuning_curves(soma_tuning_curve_normalized, model_tuning_curve_normalized)

def stim_on_neg_llikelihood(model_traces, soma_traces):
    return -1*stim_on_llikelihood(model_traces, soma_traces)

def stim_on_llikelihood(model_traces, soma_traces):
    #import pdb; pdb.set_trace()
    model_trial_amps = compute_trial_means(model_traces)
    soma_trial_amps = compute_trial_means(soma_traces)
    stim_on_ll = 0
    for j, stim in enumerate(model_traces.coords['directions']):
        #timestep = float(t)
        #for stim in 
        model_samples = model_trial_amps.isel(directions=j) #could also use .sel with the actual values
        #TODO want to use xarray this with the soma too eventually
        soma_samples = soma_trial_amps[j,:]
        this_stim_llikelihood = compute_loglikelihood(model_samples, soma_samples, method='kde')
        stim_on_ll += this_stim_llikelihood
    return stim_on_ll

def model_likelihood(model_traces, soma_traces):
    #samples should always be the last axis
    #import pdb; pdb.set_trace()
    temporal_likelihood = []
    for i, t in enumerate(model_traces.coords['samples']):
        this_sample_llikelihood = 0
        for j, stim in enumerate(model_traces.coords['directions']):
            #timestep = float(t)
            #for stim in 
            model_samples = model_traces.isel(samples=i, directions=j) #could also use .sel with the actual values
            #TODO want to use xarray this with the soma too eventually
            soma_samples = soma_traces[j,:,i]
            this_stim_llikelihood = compute_loglikelihood(model_samples, soma_samples, method='kde')
            this_sample_llikelihood += this_stim_llikelihood
        temporal_likelihood.append(this_sample_llikelihood)
    return np.mean(np.array(temporal_likelihood)), temporal_likelihood



def ll_kde(model_samples, actual_samples):
    #import pdb; pdb.set_trace()
    model_range = np.abs(max(model_samples) - min(model_samples))
    estimate_badwidth = (10*model_range/(len(model_samples)))+.01 #this should probably be based on the 25-75% not the full range to account for outliers

    # instantiate and fit the KDE model
    kde = KernelDensity(bandwidth=estimate_badwidth, kernel='gaussian')
    kde.fit(model_samples[:, None])

    # score_samples returns the log of the probability density
    logprob = kde.score_samples(actual_samples[:, None])
    return logprob



def ll_hist(model_samples, actual_samples):
    num_bins = int(np.floor(len(model_samples)/10))
    if (num_bins<1):
        raise AssertionError('Very few model samples, unlikely the model will have a sufficiently dense distribution estimated discretely.')
    
    global_min = float(min(np.min(model_samples), np.min(actual_samples))) #need to make these float, can't pass xarray object
    global_max = float(max(np.max(model_samples), np.max(actual_samples))) #need to make these float, can't pass xarray object
    #import pdb; pdb.set_trace()
    hist, bin_edges = np.histogram(np.array(model_samples), bins=num_bins, range=(float(global_min), float(global_max)), density=True)
    
    #to avoid indexing errors we remove the outer edges from the histgraom. Anything equal to those (or above but setting the min and max ensures this wont't happen)
    inner_edges = bin_edges[1:-1]
    indices = np.searchsorted(inner_edges, actual_samples, side='left') 
    likelihoods = hist[indices]
    return np.log(likelihoods)

def compute_loglikelihood_discrete_distribution(model_samples, actual_samples):
    return compute_loglikelihood(model_samples, actual_samples, method='hist')

def compute_loglikelihood(model_samples, actual_samples, method='kde'):
    if len(model_samples) < len(actual_samples):
        raise AssertionError('more test samples than model samples. unlikely the model will have a sufficiently dense distribution estimated discretely. Maybe the inputs to the function are switched?')
    #import pdb; pdb.set_trace()
    if method == 'kde':
        lls = ll_kde(np.array(model_samples), np.array(actual_samples))
    if method == 'hist':
        lls = ll_hist(model_samples, actual_samples)
        #gonna have a problem with 0s here aren't we...
    return np.sum(lls)