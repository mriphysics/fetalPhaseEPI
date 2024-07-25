# fetalPhaseEPI
Tools for phase unwrapping based 3D+T distortion correction for EPI applied to fetal MRI

This repository provides tools to perform distortion correction of EPI as described in the manuscript ''The developing Human Connectome Project fetal functional MRI release: Methods and data structures'', VR Karolis, L Cordero-Grande, AN Price, E Hughes, SP Fitzgibbon, V Kyriakopoulou, A Uus, N Harper, D Prokopenko, D Bridglal, J Willers Moore, S Wilson, M Pietsch, D Christiaens, M Deprez, L Williams, E Robinson, A Makropoulos, S-R Farahibozorg, J O'Muircheartaigh, MA Rutherford, D Rueckert, AD Edwards, T Arichi, SM Smith, E Duff, and JV Hajnal, biorXiv:doi.org/10.1101/2024.06.13.598863

The code has been developed in MATLAB and has the following structure:

###### ./
contains scripts for an example of distortion correction in fetal fMRI data: *fetalPhaseEPI_Exp1.m*, and generate a graphical illustration: *plot_Exp1.m*, and functions to call the unwrapping method: *unwrapCG.m* and revert the distortion: *undistortSinc.m*.

###### ./Tools
contains functions for visualization: *extractOrthogonalPlanes.m*, *visSegment.m*.

###### ./Tools/subtightplot
from  https://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot

###### ./Utils
contains functions that replace, extend or adapt some MATLAB built-in functions and implement generic methods: *aplGPU.m*, *blockGPU.m*, *build1DFTM.m*, *buildDFTM.m*, *buildFilter.m*, *buildStandardDFTM.m*, *dynInd.m*, *fctGPU.m*, *fftGPU.m*, *filtering.m*, *generateGrid.m*, *ifctGPU.m*, *ifftGPU.m*, *mapMat.m*, *matfun.m*, *mirroring.m*, *normm.m*, *numDims.m*, *parUnaFun.m*, *resampling.m*, *resPop.m*, *resSub.m*.

NOTE 1: Example data is provided in the dataset *x.mat* which contains the variables *x*: the distorted fMRI data, *Sequence*: one of the following 'FE' or 'SAFE', *EffectiveES*: effective inter-echo spacing (ms) for the provided FOV, *TE*: echo time or difference of echo times (ms), respectively for FE and SAFE, *voxsiz*: voxel size (mm) and factor to harmonize temporal dimension, *M* (optional): mask with area of interest for distortion correction. For runs without changing the paths, they should be placed in a folder
###### ../fetalPhaseEPIData
Data generated when running the example script is also stored in this folder as *resExp01.mat* which contains the variables *B*: the estimated field in Hz and *x*: the undistorted fMRI data.

NOTE 2: Computation times in a system with a NVIDIA RTX A6000 limited to 100W power have been below 1'.

##### Contact

Lucilio Cordero-Grande - lucilio.cordero@upm.es

##### Acknowledgments

The following projects and institutions have provided funding to support this release: The Developing Human Connectome Project, funded by the European Research Council under the European Union Seventh Framework Programme (FP/20072013)/ERC Grant Agreement no. 319456; MCIN/AEI/10.13039/501100011033/FEDER, EU under Project PID2021-911129022OA-I00; Universidad Politécnica de Madrid, providing computing resources on Magerit Supercomputer.
