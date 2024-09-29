# Instantaneous Hydrogen Bond Dynamics 
Source codes for the paper [Revisiting the Thickness of the Air-Water Interface from Two Extremes of Interface Hydrogen Bond Dynamics](https://doi.org/10.1021/acs.jctc.4c00457). 

<p align='center'>
<img src="4_plot/TOC.png" width="60%"/>
</p>

## Abstract
The air-water interface plays a crucial role in many aspects of science, because of its unique properties, such as a two-dimensional hydrogen bond (HB) network and completely different HB dynamics compared to bulk water. However, accurately determining the boundary of interfacial and bulk water, that is, the thickness of the air-water interface, still challenges experimentalists. Various simulation-based methods have been developed to estimate the thickness, converging on a range of approximately 3--10 (Å). In this study, we introduce a novel approach, grounded in density functional theory-based molecular dynamics and deep potential molecular dynamics simulations, to measure the air-water interface thickness, offering a different perspective based on prior research. To capture realistic HB dynamics in the air-water interface, two extreme scenarios of the interface HB dynamics are obtained: one underestimates the interface HB dynamics, while the other overestimates it. Surprisingly, our results suggest that the interface HB dynamics in both scenarios converges as the thickness of the air-water interface increases to 4 (Å). This convergence point, indicative of the realistic interface thickness, is also validated by our calculation of anisotropic decay of OH stretch and the free OH dynamics at the air-water interface.

## Overview
We have placed the AIMD simulation codes in `cp2k_code` and the MD simulation codes in `md_code`. The rest of the scripts and folders are used to analyze the trajectories obtained from AIMD and MD simulations. There are several folders that serve specific purposes:

- `m2_traj`: Folder to store JSON files describing the meta information of the trajectories. 
- `0_prepare`: Code to obtain the instantaneous density file (in cube file format) for a simulation trajectory.
- `1_case1`: Calculations for the first scenario.
- `2_case2`: Calculations for the second scenario. 
- `3_analyze`: Codes to analyze the two scenarios.
- `4_plot`: Code for visualizations.

## Installation
The codes are written in `bash`, `fortran` and `python`, and have been tested on Ubuntu.  To use the codes, several prerequisites need to be installed:
```bash
sudo apt install gfortran gnuplot jq
```
Here, `gfortran` is fortran compiler used to compile the Fortran code, `gnuplot` is used to plot the data, and `jq` is used to parse the JSON file.

Clone this repository:
```bash
git clone git@github.com:hg08/hb_ihb.git
cd hb_ihb
```

Create a conda environment:
```bash
conda env create -f environment.yml
conda activate ihb
```

## Run the code
Place the trajectory file and corresponding json files in `m2_traj` folder. And run these `bash` scripts in order:

1. `01_run_prepare.sh`
2. `01_run_case1.sh`
3. `01_run_case2.sh`
4. `02_run_analyse.sh`
5. `03_run_statistics.sh`
6. `04_plot.sh`

By default, the above scripts will process the trajectory corresponding to the latest JSON file in `m2_traj`. You can also specify the JSON file by adding the filename as an argument to these commands. If you have any questions about the code, feel free to [open an issue](https://github.com/hg08/hb_ihb/issues/new).

## Citation
Cited as:
> Gang Huang and Jie Huang, Revisiting the Thickness of the Air-Water Interface from Two Extremes of Interface Hydrogen Bond Dynamics, Journal of Chemical Theory and Computation,
DOI: 10.1021/acs.jctc.4c00457, 2024

Or

```
@article{huang2024,
  title = {Revisiting the Thickness of the Air-Water Interface from Two Extremes of Interface Hydrogen Bond Dynamics},
  url = {http://dx.doi.org/10.1021/acs.jctc.4c00457},
  DOI = {10.1021/acs.jctc.4c00457},
  journal = {Journal of Chemical Theory and Computation},
  publisher = {American Chemical Society (ACS)},
  author = {Huang, Gang and Huang, Jie},
  year = {2024},
  month = sep 
}
```
