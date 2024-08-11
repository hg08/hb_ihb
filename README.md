# Instantanoues hydrogen bond dynamics 

Source codes for the paper [Revisiting the Thickness of the Air-Water Interface from Two Extremes of Interface Hydrogen Bond Dynamics](https://arxiv.org/abs/2204.13941)

<p align='center'>
<img src="4_plot/TOC.png" width="60%"/>
</p>

## Abstract

The air-water interface plays a crucial role in many aspects of science, because of its unique properties, such as a two-dimensional hydrogen bond (HB) network and completely different HB dynamics compared to bulk water. However, accurately determining the boundary of interfacial and bulk water, that is, the thickness of the air-water interface, still challenges experimentalists. Various simulation-based methods have been developed to estimate the thickness, converging on a range of approximately 3 to 10 \AA. In this study, we introduce a novel approach, grounded in density functional theory-based molecular dynamics (DFTMD) and deep potential molecular dynamics (DeePMD) simulations, to measure the air-water interface thickness, offering a different perspective based on prior research. To capture realistic HB dynamics in the air-water interface, two extreme scenarios of the interface HB dynamics are obtained: one {underestimates} the interface HB dynamics, while the other {overestimates} it. Surprisingly, our results suggest that the interface HB dynamics in both scenarios converge as the thickness of the air-water interface increases to 4.0 \AA. This convergence point, indicative of the realistic interface thickness, is also validated by our calculation of anisotropic decay of OH stretch at the air-water interface.

## Installation
The codes are tested on Ubuntu. Install the following prerequisites:
```bash
sudo apt install gfortran gnuplot jq
```
Here, `gfortran` is used to compile the Fortran code, `gnuplot` is used to plot the data, and `jq` is used to parse the JSON file.

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
Place the trajectory file and corresponding json files in `m2_traj` folder.
