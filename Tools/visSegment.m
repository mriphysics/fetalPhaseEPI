function visSegment(x,y,pau,show_ri,orthv,sl,dy,str,leg,folderName,fileName)

%VISSEGMENT   Visualizes a segmentation overlaid on a given image
%   VISSEGMENT(X,{Y},{PAU},{SHOW_RI},{ORTHV},{SL},{DY},{STR},{LEG},{FOLDERNAME},{FILENAME})
%   * X is the image
%   * {Y} is the segmentation
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {SHOW_RI} to show real and imaginary instead of modulus and phase
%   * {ORTHV} shows orthogonal views
%   * {SL} serves to select a given slice or a set of slices for orthogonal
%   views
%   * {DY} serves to select a given dynamic
%   * {STR} serves to write a title
%   * {LEG} serves to write a legend
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%

N=size(x);N(end+1:5)=1;
if nargin<2;y=[];end
if nargin<3 || isempty(pau);pau=1;end
if nargin<4 || isempty(show_ri);show_ri=0;end
if nargin<5 || isempty(orthv);orthv=single(N(3)>1);end
if nargin<6;sl=[];end
if nargin<7 || isempty(dy);dy=1;end
if nargin<8 || isempty(str);str='';end
if nargin<9;leg=[];end
if nargin<10;folderName=[];end
if nargin<11;fileName=[];end

assert(numel(dy)==1,'Dynamics are non-singleton');
%assert(numel(sl)==1,'Slices are non-singleton');%This may change for
N=[N(1:3) prod(N(4:end))];
x=dynInd(reshape(x,N),dy,4);
if ~isempty(y)
    y=dynInd(reshape(y,N),dy,4);
end
if orthv
    [x,y]=extractOrthogonalPlanes(x,y,sl);
else
    if isempty(sl);sl=mod(ceil(N(3)/2)+1,N(3))+1;end
    x=dynInd(x,sl,3);
    if ~isempty(y);y=dynInd(y,sl,3);end
end
x=double(gather(x));
if ~isempty(y)
    y=double(gather(y));
    y=imquantize(y,0.5:max(y(:)))-1;
end
if any(imag(x(:)))
    if show_ri;x=cat(1,(real(x)-min(real(x(:))))/(range(real(x(:)))+eps),(imag(x)-min(imag(x(:))))/(range(imag(x(:))))+eps);
    else x=cat(1,abs(x)/max(abs(x(:))),(angle(x)/pi+1)/2);
    end
    if ~isempty(y);y=repmat(y,[2 1]);end
end
%figure
imshow(x,[],'Border','tight')
hold on
if ~isempty(y)
    C=colororder;
    for n=1:max(y(:))
        z=single(y==n);
        contour(z,[0.95 0.95],'EdgeColor',C(mod(n-1,7)+1,:),'LineWidth',2.5);
    end
end
set(gcf, 'Position', get(0,'Screensize'))
text(1,10,str,'FontSize',24,'Color',[0.8500 0.3250 0.0980],'Interpreter','latex');
if ~isempty(leg);legend(leg,'FontSize',24/2,'Interpreter','none');end
if pau==1;pause;end
if pau==2 && ~isempty(folderName) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
    if ~exist(folderName,'dir');mkdir(folderName);end
    %print(strcat(folderName,filesep,fileName),'-dpng');
    warning('off','MATLAB:prnRenderer:opengl');export_fig(strcat(folderName,filesep,fileName,'.png'));warning('on','MATLAB:prnRenderer:opengl');
    close all
end
