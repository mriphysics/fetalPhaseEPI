function [x,B]=undistortSinc(x,B,TE,ESFOV,parU,Sequence)

%UNDISTORT applies an estiamted B0 field to undistort an image
%   [X,B]=CGUNWRAP(X,B,TE,ESFOV,SEQUENCE,PARU)
%   * X is the data to be undistorted
%   * B is the estimated B0 field (in radians)
%   * TE is the echo time
%   * ESFOV is the effective inter-echo spacing at this FOV
%   * {PARU} are the parameters for undistorting
%   * {SEQUENCE} is the sequence, one of the following: 'SAFE' for Spin And
%   Field Echo (dMRI) and (default) 'FE' for Field Echo (fMRI)
%   ** X is the undistorted data
%   ** B is the field converted to Hz
%

if nargin<5;parU=[];end
if ~isfield(parU,'gibbsRingingFactor');parU.gibbsRingingFactor=1;end
if ~isfield(parU,'oversampling');parU.oversampling=2;end
if ~isfield(parU,'Verbosity');parU.Verbosity=0;end
if nargin<6 || isempty(Sequence);Sequence='FE';end

gpu=isa(x,'gpuArray');
NX=size(x);
if parU.gibbsRingingFactor>0;B=real(filtering(B,buildFilter(NX(1:2),'tukeyIso',ones(1:2),gpu,parU.gibbsRingingFactor)));end
B=convertB0Field(B,TE,ESFOV,'rad','res');
NRes=NX;NRes(2)=parU.oversampling*NRes(2);
if strcmp(Sequence,'SAFE');B=repmat(B,[1 1 2]);B=reshape(B,NX);end
x=resampling(reverseDistortion(resampling(x,NRes),real(resampling(B,NRes))),NX);
if strcmp(Sequence,'SAFE');B=dynInd(B,1:2:NX(4),4);end
B=convertB0Field(B,TE,ESFOV,'res','Hz');

end

function x=reverseDistortion(x,B)%Reverses the distortion and intensity weighting effect of a B0 field based on spectral manipulation
    gpu=isa(x,'gpuArray');

    NPE=size(x,2);
    kGrid=generateGrid(NPE,gpu,NPE,ceil((NPE+1)/2));rGrid=generateGrid(NPE,gpu);
    [DFTM,DFTMH]=buildDFTM(rGrid{1}(:),kGrid{1}(:));

    ND=numDims(x);ND=max(ND,4);
    NX=size(x);NB=size(B);
    NX(end+1:4)=1;NB(end+1:max(4,ND))=1;
    B=repmat(B,NX./NB);
    [x,NN]=resSub(x,3:ND);B=resSub(B,3:ND);NN(end+1:3)=1;%We operate per 2D images
    blSzN=50;%Block size, reduce in case of memory issues
    for n=1:blSzN:NN(3);vN=n:min(n+blSzN-1,NN(3));
        xn=dynInd(x,vN,3);Bn=dynInd(B,vN,3);
        DFTMB0=bsxfun(@times,permute(DFTM,[3 2 4 1]),exp(bsxfun(@times,-2*pi*1i*permute(kGrid{1},[2 3 4 1]),Bn)));
        xn=sum(bsxfun(@times,xn,DFTMB0),2);
        xn=aplGPU(DFTMH,permute(xn,[1 4 3 2]),2);
        x=dynInd(x,vN,3,xn);
    end
    x=reshape(x,NX);
end

function x=convertB0Field(x,TE,ES,inp,out,N)

%CONVERTB0FIELD   Converts units of B0 field
%   X=CONVERTB0FIELD(X,TE,ES,INP,OUT,{N})
%   * X is a B0 field in units given by INP
%   * TE is the echo time at which a given field-induced dephasing has been observed
%   * ES is the inter-echo spacing which scales the distortion produced by the field
%   * INP are the input units of the field
%   * OUT are the output units of the field
%   * {N} is the grid size in the phase encoding direction
%   ** X is a B0 field in units given by OUT
%

if nargin<6;N=size(x,2);end

if strcmp(inp,'Hz')
    if strcmp(out,'pix');x=x*N*ES/1000;elseif strcmp(out,'rad');x=bsxfun(@times,x,2*pi*TE/1000);elseif strcmp(out,'res');x=x*ES/1000;end
elseif strcmp(inp,'pix')
    if strcmp(out,'Hz');x=x*1000/(ES*N);elseif strcmp(out,'rad');x=bsxfun(@times,x,2*pi*TE/(ES*N));elseif strcmp(out,'res');x=x/N;end
elseif strcmp(inp,'rad')
    if strcmp(out,'Hz');x=bsxfun(@times,x,1000./(2*pi*TE));elseif strcmp(out,'pix');x=bsxfun(@times,x,N*ES./(2*pi*TE));elseif strcmp(out,'res');x=bsxfun(@times,x,ES./(2*pi*TE));end
elseif strcmp(inp,'res')
    if strcmp(out,'Hz');x=x*1000/ES;elseif strcmp(out,'pix');x=x*N;elseif strcmp(out,'rad');x=bsxfun(@times,x,(2*pi*TE)./ES);end
end

end