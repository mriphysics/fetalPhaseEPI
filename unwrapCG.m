function [x,B]=unwrapCG(x,voxsiz,parB,Sequence)

%UNWRAPCG calls an 2-norm phase unwrapping (Convex Norm Conjugate 
%Gradient) based on [1] DC Ghiglia, LA Romero. "Minimum Lp-norm 
%two-dimensional phase  unwrapping," J Opt Soc Am A, 13(10):1999-2013, Oct 
%1999, which is pre-weighted using the magnitude data as proposed in [2] DC 
%Ghiglia and LA Romero. "Robust two-dimensional weighted and unweighted
%phase unwrapping that uses fast transforms and iterative methods," J Opt
%Soc Am A, 11(1):107-117, Jan 1994, for dynamic field mapping.
%Non-iterative (free) and iterative (multiple of 2*pi) unwrapping are 
%considered based on Z Zhao et al. "Robust 2D phase unwrapping algorithm 
%based on the transport of intensity equation", Meas Sci Technol, 
%30(015201):1-8, Nov 2018
%   [X,B]=CGUNWRAP(X,{PARB},{VOXSIZ})
%   * X is the complex data used to estimate the field
%   * {VOXSIZ} are the voxel sizes of the data
%   * {PARB} are the parameters for unwrapping
%   * {SEQUENCE} is the sequence, one of the following: 'SAFE' for Spin And
%   Field Echo (dMRI) and (default) 'FE' for Field Echo (fMRI)
%   ** X is the estimation of the unwrapped field
%   ** B is the distance between the unwrapped and the wrapped field
%

ND=numDims(x);
if nargin<2 || isempty(voxsiz);voxsiz=ones(1,ND);end
ND=length(voxsiz);%Is used to establish dimensions to unwrap

if nargin<3;parB=[];end
if ~isfield(parB,'weightType');parB.weightType='Magnitude';end
if ~isfield(parB,'solverType');parB.solverType='LSIt';end
if ~isfield(parB,'weightZ');parB.weightZ=1;end
if ~isfield(parB,'weightT');parB.weightT=1;end
if ~isfield(parB,'M');parB.M=[];end
if ~isfield(parB,'Verbosity');parB.Verbosity=0;end

if nargin<4 || isempty(Sequence);Sequence='FE';end

gpu=isa(x,'gpuArray');

%INFO FOR UNWRAPPING
NX=size(x);N=NX(1:ND);NNX=prod(NX,1:ND);
if strcmp(Sequence,'SAFE');x=dynInd(x,2:2:NX(4),4).*conj(dynInd(x,1:2:NX(4),4));end%If SAFE we expect 4 dimensions at least
NX=size(x);N=NX(1:ND);NNX=prod(NX,1:ND);

%UNWRAPPING PARAMETERS
w=voxsiz;
if ND>=3;voxsiz(3)=voxsiz(3)/parB.weightZ;end
if ND>=4;voxsiz(4)=voxsiz(1)/parB.weightT;end
w=voxsiz;w(end+1:ND)=1;w=w(1:ND);
if parB.Verbosity>1;fprintf('Used weights:%s\n',sprintf(' %.2f',w));end

%PRECONDITIONER
Ps=-real(buildFilter(N*2,'2ndFiniteDiscrete',w,gpu,0,ones(1,ND)));%Second order finite difference with mirror boundary conditions
Ps(1)=1e9;%For inversion
Ps=Ps.^(-1);

x=resSub(x,ND+1:numDims(x));%Concatenate items
B=angle(x);x=abs(x);%Phase and magnitude
for l=1:size(x,ND+1)%Separately across items
    xl=dynInd(x,l,ND+1);Bl=dynInd(B,l,ND+1);
    for n=1:ND
        Bauxabs=ones(1,'like',xl);
        padEl=zeros(1,ND+1);padEl(n)=1;        
        %Few explored possibilities for weighting, more could be added, they get applied if a given string is used in weightFunc
        if contains(parB.weightType,'Magnitude');Bauxabs=bsxfun(@times,Bauxabs,sqrt(xl.*dynInd(xl,[2:NX(n) 1],n)));end%Weights given by harmonic mean of neighboring magnitudes     
        if contains(parB.weightType,'Gradient1');Bauxabs=Bauxabs.*padarray(abs(cos(wrapToPi(diff(Bl,1,n))/2)),padEl,0,'post');end%Weights by cos(dphi/2)
        if contains(parB.weightType,'Gradient2');Bauxabs=Bauxabs.*padarray(abs(cos(wrapToPi(diff(Bl,1,n))/2)).^2,padEl,0,'post');end%Weights by cos^2(dphi/2)
        if n==1;Mb=Bauxabs;else Mb=cat(ND+1,Mb,Bauxabs);end%Weights
    end

    %PHASE UNWRAPPING SOLVER    
    if strcmp(parB.solverType,'LSNonIt')%Non-iterative solver       
        Bn=real(CGsolver(inputPoisson(Bl,w),Mb,Ps,w));   
    else%Iterative solver
        [B1,K1,Bn,res1,resn1]=unwrapRes(Bl,Bl);%Residuals       
        nn=1;
        while 1
            [B1,K2,Bn,res2,resn2]=unwrapRes(res1,Bl,B1);
            if parB.Verbosity>2;fprintf('Iteration: %d / Max diff residual: %.2f / Mea diff residual: %.2f / Diff norm: %.2f / DiffPoints: %d\n',nn,max(abs(res1(:)-res2(:))),sqrt(normm(res2,res1)),abs(resn1-resn2),sum(K2(:)~=K1(:)));end
            nn=nn+1;
            if all(K2(:)==K1(:)) || abs(resn1-resn2)<1e-6*NNX || nn==100;break;end             
            K1=K2;resn1=resn2;res1=res2;
        end
    end
    x=dynInd(x,l,ND+1,Bn);
end
x=reshape(x,NX);B=reshape(B,NX);

%POSTPROCESSING B0
NDS=min(ND,3);
xref=resPop(x,1:NDS,prod(NX(1:NDS)),1);%Stack all spatial dimensions
if ~isempty(parB.M);xref=dynInd(xref,parB.M(:)~=0,1);end%Statistics over masked area
x=bsxfun(@minus,x,2*pi*round(median(xref(:))/(2*pi)));%Rounded to nearest 2Npi value
if ~isempty(parB.M);x=bsxfun(@times,x,parB.M);end%Masked

if nargout>=2
    if ~isempty(parB.M);B=bsxfun(@times,B,parB.M);end
    B=angle(exp(1i*(x-B)));
end


%UNWRAPPING RESIDUALS
function [B,K,x,res,resn]=unwrapRes(res,Bwrap,B0)               
    B=real(CGsolver(wrapToPi(encode(res,w)),Mb,Ps,w));
    if nargin>=3;B=B+B0;end
    K=round((B-Bwrap)/2/pi);
    x=Bwrap+2*K*pi;
    res=wrapToPi(x-B);
    resn=sqrt(normm(res));
end

end

%SOLVER
function x=CGsolver(y,Mb,Ps,w)    
    nX=200;%Number of iterations
    toler=2e-3;%Tolerance
    n=0;
    r=decode(y,w,Mb);n=n+1;
    
    NX=size(r);ND=numDims(r);
    x=zeros(NX,'like',r);
    z=filtering(r,Ps,1);
    p=z;
    rsold=sum(conj(z).*r,1:ND);
    l=true;
    if sqrt(min(abs(rsold(:))))<1e-9;l=false;end

    %ITERATIONS
    while l
        %SYSTEM MATRIX
        Ap=decode(encode(p,w),w,Mb);n=n+2;
    
        %UPDATES
        al=conj(rsold)./sum(conj(p).*Ap,1:ND);
        xup=bsxfun(@times,al,p);
        x=x+xup;
        xup=max(abs(xup(:)).^2);
        if xup<toler || n>=nX;break;end
        r=bsxfun(@minus,r,bsxfun(@times,al,Ap));
        z=filtering(r,Ps,1);

        rsnew=sum(conj(z).*r,1:ND);
        if sqrt(min(abs(rsnew(:))))<1e-9;break;end
        be=bsxfun(@times,rsnew,1./rsold);
        p=z+bsxfun(@times,be,p);
        rsold=rsnew;
    end
end

%DECODER
function xou=decode(x,w,Mb)
    if nargin<3;Mb=[];end
    ND=length(w);
    if ~isempty(Mb);x=bsxfun(@times,x,Mb);end    
    for n=1:ND
        padEl=zeros(1,ND+1);padEl(n)=1;
        xs=diff(padarray(dynInd(x,n,ND+1),padEl,0,'pre'),1,n)/w(n);
        if n==1;xou=xs;else xou=xou+xs;end
    end
end

%ENCODER
function xou=encode(x,w)
    ND=length(w);
    for n=1:ND
        padEl=zeros(1,ND);padEl(n)=1;
        xn=padarray(diff(x,1,n)/w(n),padEl,0,'post');
        if n==1;xou=xn;else xou=cat(ND+1,xou,xn);end
    end
end

%INPUT DATA
function xou=inputPoisson(x,w)
    ND=length(w);
    for n=1:ND   
        padEl=zeros(1,ND);padEl(n)=1;      
        xn=padarray(wrapToPi(diff(x,1,n))/w(n),padEl,0,'post');            
        if n==1;xou=xn;else xou=cat(ND+1,xou,xn);end
    end
end