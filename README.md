# allsky-point-source-detection
This repo provides code for detecting point sources on a sphere (i.e., the sky) with the continuous wavelet transform. We use this tool to detect faint gamma-ray point sources from Fermi-LAT. See my skysearch repo for a pipeline. We use the technique proposed by https://arxiv.org/pdf/astro-ph/0212578.

# Notebooks
`01` - Point Source Detection in a Flat Region with the Continuous Wavelet Transform
`02` - Partition Sky Into Approximately Flat Patches
`03` - All-Sky Point Source Detection
`04` - Results
`05` - Next Steps

# Progress (Delete Later)
1. It is better to create the first two notebooks with `pswavelets/gce/notebooks' directory and following the notes that you have already made.
2. For third notebook, just attach the python files from the `skysearch` file
3. For fourth notebook, attach allsky_notes plot and explain results
4. For 5th notebook, maybe attach the connection branch progress (connection between wavelet coefficient and signal strength)

# Installation
Run the `environment.yml` file by running the following command on the main repo directory:
```
conda env create
```
The installation works for `conda==4.14.0`. 
