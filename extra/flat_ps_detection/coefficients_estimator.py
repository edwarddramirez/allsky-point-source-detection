# For mexh2, we loaded functions.py -> functions_2.py and saved coefficients in ".../skymap_sps_..." -> ".../skymap_sps_2_..."
# and scale = 0.462 degree to be consistent with Zhong choice

import timeit as timeit
import sys
import numpy as np
import functions as fcts

ti = timeit.default_timer()

# load local directory
username="ramirez"
local_dir = "/het/p4/"+username+"/gcewavelets/cwt_v0.1_play/"

# load data
data_dir = local_dir + 'data/skymap/'
x_data, y_data = np.load(data_dir + 'sample_skymap_upper_wps_2.npy', allow_pickle = True)
arr_data = np.vstack((x_data, y_data)).T

# define 2d wavelet 
mexh = fcts._2d_wavelet('mexh')

# define range of translations (b) and scales (a) to compute cwt
results_dir = local_dir + 'results/'

## grid of b-values
# Due to limited memory, we compute the coefficients corresponding to a single translation
# v0.1: Due to farm priority, we compute coefficients on a horizontal line of constant by

# step_size = 0.005
# arr_bx_plot = np.arange(0,1+step_size,step_size)
# arr_bx = 0.5 * (arr_bx_plot[:-1] + arr_bx_plot[1:])
# Nx = len(arr_bx)

step_size = 0.005
arr_bx_plot = np.arange(-0.2,1.2+step_size,step_size, dtype = float)
arr_bx = 0.5 * (arr_bx_plot[:-1] + arr_bx_plot[1:])
Nx = len(arr_bx)

by = float(sys.argv[1]) 
arr_by = by * np.ones(Nx)

arr_b_line = np.stack((arr_bx, arr_by), axis = -1) # convert to datapoints

## grid of a-values
# for larger datasets, will only be able to compute one value of a
# arr_a = np.load(results_dir + 'arr_a.npy', allow_pickle = True)
log_arr_a_edge = np.linspace(-3,1,1000+1)
arr_a_edge = 10**log_arr_a_edge
arr_a = 0.5 * (arr_a_edge[:-1] + arr_a_edge[1:])

Nx = len(arr_bx)
Na = len(arr_a)
buf_arr_coefficients = np.zeros([Nx, 1, Na])

for nx in range(Nx):
    arr_b = arr_b_line[nx]
    # Step 1. Compute the arg(\psi) = (x - b) / a
    ## prepare broadcasting
    buf_arr_data = arr_data[np.newaxis,np.newaxis,np.newaxis]
    buf_arr_b = arr_b[np.newaxis, np.newaxis, np.newaxis, np.newaxis]
    buf_arr_a = arr_a[np.newaxis,np.newaxis,:,np.newaxis,np.newaxis]

    # broadcasting
    buf_arr_arg =  (buf_arr_data - buf_arr_b) / buf_arr_a

    # Note that indices of buf_f_output will be
    #  1st: bx index (axis)
    #  2nd: by index
    #  3rd: a index
    #  4th: data index
    #  5th: vector index of b or data

    # Step 2. Calculate \psi(x-b/a)
    buf_mexh_output = mexh.base_fct(buf_arr_arg)

    # Step 3. Estimate wavelet coefficient through sum over data
    ## F(b,a) ~ \sum_i(\psi^{b,a}(x_i)) / a / N

    # remove two dimensions of a-array to divide mexh
    buf_arr_a_sq = np.squeeze(buf_arr_a, axis = -1)
    buf_arr_a_sq = np.squeeze(buf_arr_a_sq, axis = -1)

    # estimate wavelet coefficient through sum
    coefficient_estimate = np.sum(buf_mexh_output, axis = -1) / buf_arr_a_sq / len(arr_data)
    
    # update coefficient bx-array
    buf_arr_coefficients[nx,0,:] = coefficient_estimate[0,0,:]

# save coefficient estimate
# unprocessed_data_dir = results_dir + 'preprocessed/'
str_by = str.format('{0:.5f}',by)
file_name = 'coefficient' + '_' + str_by

np.save(file_name, buf_arr_coefficients)

tf = timeit.default_timer()
print('Run Time (min): ' + str( (tf - ti) / 60 ) )