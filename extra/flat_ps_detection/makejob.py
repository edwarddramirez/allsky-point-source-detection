#/usr/bin/env python
#

import sys, os, time, fileinput
import numpy as np
import glob

#username
username="ramirez"

#make folder for python scripts
executedir = "batch_exec"
os.system("rm -rf "+executedir)
os.system("mkdir "+executedir)

# data file for iteration
# step_size = 0.001
# arr_bx_plot = np.arange(0,1+step_size,step_size, dtype = float)
# arr_by_plot = np.arange(0,1+step_size,step_size, dtype = float)
# arr_bx = 0.5 * (arr_bx_plot[:-1] + arr_bx_plot[1:])
# arr_by = 0.5 * (arr_by_plot[:-1] + arr_by_plot[1:])

# define mesh grid
# mesh_bx, mesh_by = np.meshgrid(arr_bx,arr_by) # each output array (NxN shaped) contains x or y value at given (i,j)-th position
# mesh_bxby = np.stack((mesh_bx, mesh_by), axis=-1)

# Nx = len(arr_bx)
# arr_b = mesh_bxby.reshape(Nx*Nx,2) # flatten to 2D array

# arr_b = []
# for ny in range(Nx):
#     by = arr_bx[ny]
#     for nx in range(Nx):
#         bx = arr_by[nx]
#         b = np.array([bx,by])
#         print(ny)
#         arr_b.append(arr_b)

# arr_b = np.array(arr_b)

# step_size = 0.005
# arr_by_plot = np.arange(0,1+step_size,step_size, dtype = float)
# arr_by = 0.5 * (arr_by_plot[:-1] + arr_by_plot[1:])

step_size = 0.005
arr_by_plot = np.arange(-0.2,1.2+step_size,step_size, dtype = float)
arr_by = 0.5 * (arr_by_plot[:-1] + arr_by_plot[1:])

np.savetxt('by_list', arr_by, fmt = '%.5f', delimiter = ', ') # fmt specifies integer format
os.system("mv by_list batch_exec") # send to batch file for loading

#open file to contain submission commands for jdl files
dofile = open(executedir+"/do_all.src",'w')

basename = 'skymap_'

# identify directory containing results produced by main script
RESULTS_DIR = "/het/p4/"+username+"/gcewavelets/cwt_v0.1_play/results/preprocessed/skymap_up_wps_2"

# make directory if it doesn't already exist
os.system("mkdir -p "+RESULTS_DIR)
    
#define the local scratch folder
localdir = "$_CONDOR_SCRATCH_DIR"
gen_command = ( 'python '+'/het/p4/'+username+
               '/gcewavelets/cwt_v0.1_play/coefficients_estimator.py')

runname = basename

# key difference: iterate over exec.sh file, not the python file itself!
execfilename = "exec"+".sh"
executefile = open(executedir+"/"+execfilename,'w')
executefile.write("#!/bin/bash\n")
executefile.write("export VO_CMS_SW_DIR=\"/cms/base/cmssoft\"\n")
executefile.write("export COIN_FULL_INDIRECT_RENDERING=1\n")
executefile.write("export SCRAM_ARCH=\"slc5_amd64_gcc434\"\n")
executefile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")

#copy template directory to new location, and update its random number seed and run name
executefile.write("cd "+localdir+"\n")
#    executefile.write("mkdir -p results\n")
executefile.write("by=$1\n")
executefile.write(gen_command+' '+'$by'+'\n')
#    executefile.write('tar -czvf '+runname+'.tar.gz results/*\n')
#    executefile.write('cp *tar.gz '+RESULTS_DIR+'\n')
executefile.write('cp *npy '+RESULTS_DIR+'\n')
executefile.close()
os.system("chmod u+x "+executedir+"/"+execfilename)

#write out jdl script for job submission
jdlfilename = "exec"+".jdl.base"
jdlfile = open(executedir+"/"+jdlfilename,'w')
jdlfile.write("universe = vanilla\n")
jdlfile.write("+AccountingGroup = \"group_rutgers."+username+"\"\n")
jdlfile.write("Arguments = $(by)\n")
jdlfile.write("Executable = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/"+executedir+"/"+execfilename+"\n")
jdlfile.write("getenv = True\n")
jdlfile.write("should_transfer_files = NO\n")
jdlfile.write("priority = 0\n")
jdlfile.write("Output = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/condor/"+runname+'$(by)'+".out\n")
jdlfile.write("Error = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/condor/"+runname+'$(by)'+".err\n")
jdlfile.write("Log = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/condor/script.condor\n")
jdlfile.write("max_materialize = 500\n") # needs to be placed before queue
jdlfile.write("queue by from by_list\n") # data file should not have format
jdlfile.close()

dofile.write("condor_submit "+jdlfilename+"\n")

print('Done!')

dofile.close()