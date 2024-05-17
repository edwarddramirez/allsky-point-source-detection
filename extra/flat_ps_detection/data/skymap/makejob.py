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

step_size = 5
counts_list = np.arange(10,50+step_size,step_size, dtype = int)

np.savetxt('counts_list', counts_list, fmt = '%.0f', delimiter = ', ') # fmt specifies integer format
os.system("mv counts_list " + executedir) # send to batch file for loading

#open file to contain submission commands for jdl files
dofile = open(executedir+"/do_all.src",'w')

basename = 'newmaps_'

# identify directory containing results produced by main script
RESULTS_DIR = "/het/p4/"+username+"/gcewavelets/cwt_v0.1_play/data/skymap/maps_by_counts"

# make directory if it doesn't already exist
os.system("mkdir -p "+RESULTS_DIR)
    
#define the local scratch folder
localdir = "$_CONDOR_SCRATCH_DIR"
gen_command = ( 'python '+'/het/p4/'+username+
               '/gcewavelets/cwt_v0.1_play/data/skymap/farm_skymap_generator.py')

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
executefile.write("counts=$1\n")
executefile.write(gen_command+' '+'$counts'+'\n')
#    executefile.write('tar -czvf '+runname+'.tar.gz results/*\n')
#    executefile.write('cp *tar.gz '+RESULTS_DIR+'\n')
executefile.write('cp -r * '+RESULTS_DIR+'\n')
executefile.close()
os.system("chmod u+x "+executedir+"/"+execfilename)

#write out jdl script for job submission
jdlfilename = "exec"+".jdl.base"
jdlfile = open(executedir+"/"+jdlfilename,'w')
jdlfile.write("universe = vanilla\n")
jdlfile.write("+AccountingGroup = \"group_rutgers."+username+"\"\n")
jdlfile.write("Arguments = $(counts)\n")
jdlfile.write("Executable = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/data/skymap/"+executedir+"/"+execfilename+"\n")
jdlfile.write("getenv = True\n")
jdlfile.write("should_transfer_files = NO\n")
jdlfile.write("priority = 0\n")
jdlfile.write("Output = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/condor/"+runname+'$(counts)'+".out\n")
jdlfile.write("Error = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/condor/"+runname+'$(counts)'+".err\n")
jdlfile.write("Log = /het/p4/"+username+"/gcewavelets/cwt_v0.1_play/condor/script.condor\n")
# jdlfile.write("max_materialize = 500\n") # needs to be placed before queue
jdlfile.write("queue counts from counts_list\n") # data file should not have format
jdlfile.close()

dofile.write("condor_submit "+jdlfilename+"\n")

print('Done!')

dofile.close()