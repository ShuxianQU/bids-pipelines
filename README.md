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

If you find the toolbox helpful, please consider citing the following ISMRM abstract:

Qu, S., et al. Advancing whole-brain BOLD fMRI in humans at 10.5 Tesla with motion-robust 3D EPI and RF parallel transmission: initial experience, ISMRM 2025
 
