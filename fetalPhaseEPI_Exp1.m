%FETALPHASEEPI_EXP1 script performs distortion correction of EPI as 
%described in the manuscript ''The developing Human Connectome Project 
%fetal functional MRI release: Methods and data structures'', V.R. Karolis, 
%L. Cordero-Grande, A.N. Price, E. Hughes, S.P. Fitzgibbon, 
%V. Kyriakopoulou, A. Uus, N. Harper, D. Prokopenko, D. Bridglal, 
%J. Willers Moore, S. Wilson, M. Pietsch, D. Christiaens, M. Deprez, 
%L. Williams, E. Robinson, A. Makropoulos, S.-R. Farahibozorg, 
%J. O'Muircheartaigh, M.A. Rutherford, D. Rueckert, A.D. Edwards, 
%T. Arichi, S.M. Smith, E. Duff, and J.V. Hajnal

clearvars
%EXPERIMENT PARAMETERS
parE.Plot=0;%0 to save results / 1 to save and plot results
parE.Verbosity=2;%Level of verbosity, from 0 to 3

curFolder=fileparts(mfilename('fullpath'));
addpath(genpath(curFolder));%Add code
pathData=strcat(curFolder,'/../fetalPhaseEPIData');%Data path

%READ DATA
if parE.Verbosity>0
    fprintf('Illustration of distortion correction for fetal fMRI data\n');
    fprintf('Reading input data...\n');tsta=tic;
end
load(fullfile(pathData,'x.mat'));%SMS FE fetal fMRI
if parE.Verbosity>0;fprintf('Finished reading input data in %.3fs\n',toc(tsta));end

%BO ESTIMATION PARAMETERS
parB.weightZ=1/10;%Weight of the slice information for unwrapping
parB.weightT=1/100;%Weight of the temporal information for unwrapping
parB.weightType='MagnitudeGradient2';%Type of weighting
parB.solverType='LSNonIt';%Type of solver
if exist('M','var');parB.M=M;else parB.M=[];end
parB.Verbosity=parE.Verbosity;%Level of verbosity, from 0 to 3

gpu=single(gpuDeviceCount && ~blockGPU);%Detects whether gpu computations are possible
if gpu;x=gpuArray(x);parB.M=gpuArray(parB.M);end

%SOLVE FOR B0
if parE.Verbosity>0;fprintf('Estimating B0...\n');tsta=tic;end  
B=unwrapCG(x,voxsiz,parB,Sequence);
if parE.Verbosity>0;fprintf('Finished estimating B0 in %.3fs\n',toc(tsta));end

%DISTORTION CORRECTION PARAMETERS
parU.gibbsRingi=1;%Fraction parameter of Gibbs ringing Tukey filter for smoothing the field
parU.oversampling=2;%Oversampling in resolution for reversing the gradient
parU.Verbosity=parE.Verbosity;%Level of verbosity, from 0 to 3

%UNDISTORT
if parE.Verbosity>0;fprintf('Reversing distortion...\n');tsta=tic;end  
[xu,B]=undistortSinc(x,B,TE,EffectiveES,parU,Sequence);
if parE.Verbosity>0;fprintf('Finished reversing distortion in %.3fs\n',toc(tsta));end

%WRITE RESULTS
if parE.Verbosity>0;fprintf('Saving results...\n');end
xu=gather(xu);B=gather(B);
save(fullfile(pathData,'resExp01.mat'),'xu','B');
if parE.Verbosity>0;fprintf('Finished saving results in %.3fs\n',toc(tsta));end

%PLOT RESULTS
if parE.Plot;plot_Exp1;end
