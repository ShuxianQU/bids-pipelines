# bids-pipelines
This is a matlab script for using BIDS (Brain Imaging Data Structure) apps to process BOLD functional MRI (fMRI) using offline recon such as those with amri_epi, where no dicoms are available. 
More info on BIDS can be found at https://bids.neuroimaging.io. Overall, this pipeline uses sMRIPrep and fMRIPrep for structural MRI faciliated fMRI preprocessing, followed by fMRI postprocessing using XCP-D. 

The pipeline fulfills the following steps:
1) create bids input data structure using dcm2bids (https://unfmontreal.github.io/Dcm2Bids);
2) perform structural MRI preprocessing using sMRIPrep (https://www.nipreps.org/smriprep);
3) conduct fMRI preprocessing using fMRIPrep (https://fmriprep.org/en/latest);
4) carry out functional connectivity related postprocessing using XCP-D (https://xcp-d.readthedocs.io/en/latest/index.html).

The pipeline assumes MP2RAGE to be used to define subject's native T1w space. 
For compatibility with the sMRIPrep pipeline which is optimized for MPRAGE, MP2RAGE is MPRAGEized (mostly to clean up the UNI background) using tools shared at: https://github.com/srikash/3dMPRAGEise

We have demonstrated the utility of this pipeline for resting state BOLD functional MRI (fMRI) at 10.5 T and reported our findings in the following paper: 

Advancing whole-brain BOLD fMRI in humans at 10.5 Tesla with motion-robust 3D EPI, parallel transmission and high-density RF receive coils. Shuxian Qu, Jiaen Liu, Peter van Gelderen, Jacco A. de Zwart, Jeff H. Duyn, Matt Waks, Russell Lagore, Alexander Bratch, Andrea Grant, Edward Auerbach, Lance Delabarre, Alireza Sadeghi-Tarakameh, Yigitcan Eryaman, Gregor Adriany, Kamil Ugurbil, and Xiaoping Wu. MRM 2025.
 

## Demonstrations
To grab an idea of how the pipieline works in fMRI processing, you may run the demo script, `bids_pipelines.m`. 
Please create a folder 'bids' first. For this quick demonstration, you will need to download the low resolution example fMRI data (2 mm isotropic, 10 volumes) and the corresponding anatomical MP2RAGE images (1 mm isotropic) from the subfolder "bold-data" shared at [this google drive](https://drive.google.com/drive/folders/1cVI2BXiPV-lKmIz1KD7RiYVmy8S9kSTL?usp=drive_link). Note that the demo script assumes that the example fMRI and structural data are stored under "~/myData/bold-data". 
