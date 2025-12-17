"""description of this module, interface.py(interface).
please update this docstring as you develop.

interface.py should contain functions that corresponde to each keyword in the "weighted_by" columns of the model_specification.csv. 
Each keyword needs to map to a function here that will return the right vector of n weights for n spines
"""
from scipy.optimize import minimize
import numpy as np

def somatic_identity(input_traces, actual_soma_traces, score_function):
    #we apply the same normalization that we would eventually on the somatic traces - its just an estimate of m and b but it is enough     
    #norm_traces =  normalize_traces_by_tuning_curve(input_traces)
    #score = score_function(norm_traces, actual_soma_traces)
    return input_traces, None, None








def somatic_1_param(input_traces, actual_soma_traces, score_function):
    m = 1
    guessed_params = [m]
    #bounds = [(0,100)]
    result = minimize(apply_error_func, guessed_params,
                  args=(somatic_1_fit, input_traces, actual_soma_traces, score_function),
                  #bounds = bounds,
                  )
    fit_params = result.x
    best_score = result.fun
    return best_score, fit_params

def apply_error_func(guessed_params, func_to_fit, input_traces, actual_soma_traces, score_function):
    #print('guessed_params', guessed_params)
    #print('input traces', input_traces)
    output_traces = func_to_fit(guessed_params, input_traces)
    #print('output traces', output_traces)
    #print('soma traces', actual_soma_traces)
    score = score_function(output_traces, actual_soma_traces)
    #print('score', score)
    return score

def somatic_1_fit(guessed_params, input_traces):
    [m] = guessed_params
    output_traces = m*input_traces
    return output_traces

def somatic_2_fit(guessed_params, input_traces):
    #unpack_guessed_params
    [m, x_intercept] = guessed_params
    output_traces = m*(input_traces - x_intercept)
    return output_traces

def somatic_relu_fit(guessed_params, input_traces):
    #unpack_guessed_params
    [m, x_intercept] = guessed_params
    output_traces = m*(input_traces - x_intercept)
    output_traces = np.maximum(output_traces, 0)
    return output_traces

def somatic_2_param_relu(input_traces, actual_soma_traces, score_function):
    return somatic_2_param(input_traces, actual_soma_traces, score_function,  fit_func = somatic_relu_fit)

def somatic_2_param_lin(input_traces, actual_soma_traces, score_function):
    return somatic_2_param(input_traces, actual_soma_traces, score_function,  fit_func = somatic_2_fit)

def somatic_2_param(input_traces, actual_soma_traces, score_function, fit_func = somatic_2_fit):
    guessed_params = [1,np.mean(input_traces)]
    result = minimize(apply_error_func, guessed_params,
                  args=(fit_func, input_traces, actual_soma_traces, score_function),
                  #bounds = bounds,
                  )
    fit_params = result.x
    simulated_output_of_soma = fit_func(fit_params, input_traces)
    best_score = result.fun
    return simulated_output_of_soma, best_score, fit_params

def somatic_doublerelu_fit(guessed_params, input_traces):
    #unpack_guessed_params
    [m, x_intercept, x2] = guessed_params
    y_at_x2 = m*(x2 - x_intercept)
    output_traces = m*(input_traces - x_intercept)
    output_traces[output_traces<0] = 0
    output_traces[input_traces>(x2)] = y_at_x2
    return output_traces


def somatic_sigmoid_fit(guessed_params, input_traces):
    #unpack_guessed_params
    [lateral_shift, flattness, basline_shift] = guessed_params
    x = input_traces
    y = 1/(1+np.exp(-(x-lateral_shift)/flattness))
    output_traces = (y+10*basline_shift)/(1+10*basline_shift)
    return output_traces

def somatic_3_param_sigmoid(input_traces, actual_soma_traces, score_function):
    return somatic_3_param(input_traces, actual_soma_traces, score_function, fit_func = somatic_sigmoid_fit)

#def somatic_3_param_doublerelu(input_traces, actual_soma_traces, score_function, fit_func):
#    return somatic_3_param(input_traces, actual_soma_traces, score_function, fit_func = somatic_doublerelu_fit)

def somatic_3_param(input_traces, actual_soma_traces, score_function, fit_func = somatic_doublerelu_fit):
    guessed_params = [.5,.1,.02]
    bounds = ([0,1],[0,1], [0,1] )
    result = minimize(apply_error_func, guessed_params,
                  args=(fit_func, input_traces, actual_soma_traces, score_function),
                  bounds = bounds,
                  )
    fit_params = result.x
    simulated_output_of_soma = fit_func(fit_params, input_traces)
    best_score = result.fun
    return simulated_output_of_soma, best_score, fit_params



#want to use some of this to fit a sigmoid eventually
#def somatic_3_param(input_traces, actual_soma_traces, score_function):
#    lateral_shift = 0.9 #@param {type:"slider", min:0, max:1, step:0.1}
#    flattness = 0.13 #@param {type:"slider", min:0, max:0.3, step:0.01}
#    basline_shift = 0.06 #@param {type:"slider", min:0, max:1, step:0.01}#
#
#    guessed_params = [lateral_shift, flattness, basline_shift]
#    bounds = [(0,100)]
#    result = minimize(apply_error_func, guessed_params,
#                  args=(somatic_3_fit, input_traces, actual_soma_traces, score_function),
#                  #bounds = bounds,
#                  )
#    fit_params = result.x
#    best_score = result.fun
#    return best_score, fit_params#

#def error_func_3_param(guessed_para#ms, input_traces, actual_soma_traces, score_function):
#    #unpack_guessed_params
#    #will probably take some debugging to make sure this is fit properly
#    [lateral_shift, flattness, basline_shift] = guessed_params
#    x = input_traces
#    y = 1/(1+np.exp(-(x-lateral_shift)/flattness))
#    output_traces = (y+10*basline_shift)/(1+10*basline_shift)
#    score = score_function(output_traces, actual_soma_traces)
#    return score
