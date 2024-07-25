%GENERATES AN ANIMATION WITH RESULTS

clearvars
curFolder=fileparts(mfilename('fullpath'));
addpath(genpath(curFolder));%Add code
pathData=strcat(curFolder,'/../fetalPhaseEPIData');%Data path
load(fullfile(pathData,'x.mat'),'x','Sequence');%MB fMRI
load(fullfile(pathData,'resExp01.mat'),'xu','B');

NX=size(x);
if strcmp(Sequence,'SAFE');x=dynInd(x,1:2:NX(4),4);xu=dynInd(xu,1:2:NX(4),4);end%If SAFE we expect 4 dimensions at least, this extracts the SE
NX=size(x);
sl=[24 30 36];%Three slices arranged vertically
x=dynInd(x,sl,3);xu=dynInd(xu,sl,3);B=dynInd(B,sl,3);
B=B-min(B(:));B=B/max(B(:));B=permute(B,[1 3 2 4]);B=reshape(B,[NX(1)*3 NX(2) NX(4)]);
x=abs(x);x=x/max(x(:));x=permute(x,[1 3 2 4]);x=reshape(x,[NX(1)*3 NX(2) NX(4)]);
xu=abs(xu);xu=xu/max(xu(:));xu=permute(xu,[1 3 2 4]);xu=reshape(xu,[NX(1)*3 NX(2) NX(4)]);

figure
s=1;NT=size(x,3);
while 1
    subtightplot(1,3,1)
    visSegment(abs(x(:,:,s)),[],0,0);title('Original')
    subtightplot(1,3,2)
    visSegment(abs(xu(:,:,s)),[],0,0);title('Undistorted')
    subtightplot(1,3,3)
    visSegment(B(:,:,s),[],0,1);title('Estimated field')
    drawnow
    pause(2)
    s=mod(s,NT)+1;
end