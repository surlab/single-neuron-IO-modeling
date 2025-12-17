from src import data_io as io
from src import plotting as plot
from src import computation as comp
from src import config as cfg
from src import helper_functions as hf

#imports
import xarray as xr
#import pandas as pd
import numpy as np
import ipdb


from matplotlib import pyplot as plt

from PIL import Image #this needs to be after matplotlib??
from scipy.stats import stats   


import os
#from neuron import h, gui



#define paths    
test_path = os.path.join('demo_data', 'test.txt')
print(test_path)
print(os.getcwd())
print(os.path.exists(test_path))


#soma_path = "/Users/User/code/Adult-Spine-Models/scripts/demo_data/ASC26_cell_3_soma.mat"
#spines_path = "/Users/User/code/Adult-Spine-Models/scripts/demo_data/ASC26_cell_3_spines.mat"
spines_path = r"C:\Users\User\Dropbox (MIT)\2021 Gregg Sur rotation\ASC_experimental_data\2023-08 Spine Data\ASC10.mat"
#print(os.path.exists(soma_path))
print(os.path.exists(spines_path))


#data inputs
#io.readfile(test_path)
#soma = io.loadmat(soma_path)
spine_data = io.loadmat(spines_path)
#soma_data = soma


def one_over_peak_height(spine_data):
    spine_data_xr, spines_per_fov_list = comp.compile_spine_traces(spine_data)
    return 1/compute_peak_heights(spine_data_xr)

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


weights = one_over_peak_height(spine_data)
print(weights.shape)
print(weights)