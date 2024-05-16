# allsky-point-source-detection
This repo provides code for detecting point sources on a sphere (i.e., the sky) with the continuous wavelet transform. We use this tool to detect faint gamma-ray point sources from *Fermi*-LAT. See my skysearch repo for a pipeline. We use the technique originally proposed by [arXiv:0212578](https://arxiv.org/pdf/astro-ph/0212578).

# Notebooks
1. `01_ps_detection_cwt_flat.ipynb` - Point Source Detection in a Flat Region with the Continuous Wavelet Transform
2. `02_partition_sphere_into_flat_region.ipynb` - Partition Sky into Flat Patches and Project Data into Flat Patches
3. `03_ps_detection_cwt_spherical_patch.ipynb` - Point Source Detection in Spherical Patch with the Continuous Wavelet Transform (Bare-Bones Version of [skysearch](https://github.com/edwarddramirez/skysearch) repo)
4. `04_skysearch_and_next_steps.ipynb` - Discussion of [skysearch](https://github.com/edwarddramirez/skysearch) Repo and Future Steps to Take

# Installation
Run the `environment.yml` file by running the following command on the main repo directory:
```
conda env create
```
The installation works for `conda==22.9.0`. 
