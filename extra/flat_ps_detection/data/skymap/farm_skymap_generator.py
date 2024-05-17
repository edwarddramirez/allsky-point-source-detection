import numpy as np
from scipy.stats import norm, uniform
import sys, os, time, fileinput

# load local directory
username="ramirez"
local_dir = "/het/p4/"+username+"/gcewavelets/cwt_v0.1_play/"

# load data
data_dir = local_dir + 'data/skymap/'
phi_events, theta_events = np.load(data_dir + 'sample_skymap_upper_angular.npy', allow_pickle=True)

# center offset data to [-15,15] x [-15,15]
phi_centered = phi_events + 10 * np.pi / 180
theta_centered = theta_events - 25 * np.pi / 180

x_data = ( phi_centered * 180 / np.pi / 15 + 1 ) / 2
y_data = ( (theta_centered - np.pi/2) * 180 / np.pi / 15 + 1 ) / 2
data = np.vstack((x_data, y_data)).T

# generate 1000 maps with 1 point source injected
## locations are chosen so that the edges of the point sources do not "touch"
## the region within 0.05 units from the boundary
## [1,2,3,4,5,6] = [100, 50, 30, 20, 15, 10]

N_counts = int(sys.argv[1]) 

file_name = 'map_'

# results_dir
results_dir = "skymaps_sps_" + str(N_counts) + '/'
os.system("mkdir -p "+ results_dir)

Nps = 1000

xloc_list = 0.7 * np.random.random(size = Nps) + 0.15
yloc_list = 0.7 * np.random.random(size = Nps) + 0.15

psf_degree = 0.4
psf_scale = psf_degree / 15  / 2

for nx in range(Nps):
    xloc = xloc_list[nx]
    nm_x = norm(scale = psf_scale, loc = xloc)

    yloc = yloc_list[nx]
    nm_y = norm(scale = psf_scale, loc = yloc)

    x_ps = nm_x.rvs(size = N_counts)
    y_ps = nm_y.rvs(size = N_counts)
        
    x_data_sps = np.concatenate((x_data, x_ps))
    y_data_sps = np.concatenate((y_data, y_ps))
    
    np.save(results_dir + file_name + str(nx), [x_data_sps, y_data_sps])
np.save(results_dir + 'ps_locs', [xloc_list, yloc_list])