function [z,w]=extractOrthogonalPlanes(x,y,par)

%EXTRACTORTHOGONALPLANES   Gets the orthogonal planes from a volume for 
%segmentation visualization
%   [Z,W]=EXTRACTORTHOGONALPLANES(X,{Y},{PAR})
%   * X is the image
%   * {Y} is the segmentation
%   * {PAR} serves to slice
%   ** Z contains the extracted planes for the image
%   ** W contains the extracted planes for the segmentation
%

x=x(:,:,:,1);
gpu=isa(x,'gpuArray');

N=size(x);
x=resampling(x,max(N)*ones(1,3),3);
M=size(x);
if nargin<3 || isempty(par)
    if ~isempty(y)
        y=y(:,:,:,1);
        y=resampling(y,max(N)*ones(1,3),3);  
        C=gridv(generateGrid(M,gpu,M,zeros(1,3)));   
        yy=y>0.5;
        par=sum(C.*yy,1:3)./sum(yy,1:3);par=par(:)';
    else
        par=mod(ceil(M/2)+1,M)+1;
    end
else
    par=par+ceil((M-N)/2);
end
z=zeros([M(1) M(2) 3],'like',x);
if ~isempty(y);w=zeros([M(1) M(2) 3],'like',y);end
for n=1:3
    x=shiftdim(x,1);
    if par(n)<=0 || par(n)>size(x,3);par(n)=ceil(size(x,3)/2);end
    z(:,:,n)=dynInd(x,round(par(n)),3);
    if ~isempty(y)
        y=shiftdim(y,1);
        w(:,:,n)=dynInd(y,round(par(n)),3);
    end
end
z=z(:,:);
if ~isempty(y);w=w(:,:);else w=[];end
