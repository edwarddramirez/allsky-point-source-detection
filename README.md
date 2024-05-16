# allsky-point-source-detection
This repo provides code for detecting point sources on a sphere (i.e., the sky) with the continuous wavelet transform. We use this tool to detect faint gamma-ray point sources from Fermi-LAT. See my skysearch repo for a pipeline. We use the technique proposed by https://arxiv.org/pdf/astro-ph/0212578.

# Notebooks
1. `01_ps_detection_cwt_flat.ipynb` - Point Source Detection in a Flat Region with the Continuous Wavelet Transform
   - See `/het/p4/ramirez/gcewavelets/cwt_v0.1_play/test_statistic_icwt.ipynb` and PPT files
3. `02_partition_sphere_into_flat_region.ipynb` - Partition Sky Into Approximately Flat Patches
   - See most of the notebooks in `/mnt/c/Users/Edwar/OneDrive - Rutgers University/projects/pswavelets/gce/notebooks`
5. `03.ipynb` - Generate Data in the Sphere (Background Photons + Point Sources)
   - See `/mnt/c/Users/Edwar/OneDrive - Rutgers University/projects/pswavelets/gce/notebooks/skymap_generation_.ipynb`
   - See `/het/p4/ramirez/gcewavelets/skysearch/data`
   - Or generate from scratch with `generate_skymap.py`
7. `04.ipynb` - All-Sky Point Source Detection
   - Consider just going through the contents of the python files
9. `05` - Results
   - Consider just copy and pasting `allsky_notes_plots` and comment
11. `06` - Next Steps
    - Go through next steps described in different branches and in later notes

# Progress (Delete Later)
1. It is better to create the first three notebooks with `pswavelets/gce/notebooks' directory and following the notes that you have already made.
2. For third notebook, just attach the python files from the `skysearch` file
3. For fourth notebook, attach allsky_notes plot and explain results
4. For 5th notebook, maybe attach the connection branch progress (connection between wavelet coefficient and signal strength)

# Installation
Run the `environment.yml` file by running the following command on the main repo directory:
```
conda env create
```
The installation works for `conda==22.9.0`. 
