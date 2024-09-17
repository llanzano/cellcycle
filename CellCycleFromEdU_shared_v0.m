%%   User-fiendly code based on the following article:
%$£   
%$£  Super-resolved analysis of colocalization between replication and transcription along the cell cycle
%$£  in a model of oncogene activation
%$£  Anna Provvidenza Privitera, Silvia Scalisi, Greta Paternò, Elena Cerutti, Morgana D’Amico, Pier Giuseppe Pelicci, 
%$£  Mario Faretta, Gaetano Ivan Dellino, Alberto Diaspro, Luca Lanzanò
%$£  Communications Biology 2024
%$£  
%$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£ 
%$£                                                                     %$£
%$£                              Luca Lanzanò                           %$£
%$£       University of Catania - Department of Physics and Astronomy   %$£
%$£           Istituto Italiano di Tecnologia - Nanoscopy               %$£
%$£                      User-Friendly Version (........)               %$£
%$£                                                                     %$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£


function [ ]=CellCycleFromEdU_shared_v0(filename );

% goal
% 
% count pixel/intensity in edu channel (must be first channel)
% load: count mask nuclei, edu mask, intensity image 
% open stacks
% output text file contains: 
% cell analyzed, cell class, Intensity of replication foci,  pixel density of replication foci , file index, stack frame, cell index in count mask)
% Cell CLASS  0=G1/G2  1=Early  2=Mid(…) 3=Late 4=unclassified
% the text file can be loaded into the excel template
%

% SORTING of cells 
%parameters
savefigure=0;
densfactor=0.4;
IntFactor1=3;
IntFactor2=3;
kslopeE=0;
kslopeL=0;
minpx=50;  % minimum pixels for EdU signal
SizeLimit=1.2;
celltest=0;  % cell # for test mode;  0 for 

DNAnormTHR=0.2;  % 0.2 confocal si SBG; 0.5 confocal no sBG;  0.3 STED con SBG

%green colormap
greenmap = zeros(51,3);
greenmap(:,2)=[0:0.02:1]';

cellfirst=1;
cellmax=100000; % max number of cells to analyze
AllResults=double.empty(20,0);

% nch=2;
chDNA=3; % channel for DNA (DAPI)
tag=[''];
mask=0;  % mask to use as ROI 
maskedu=1;
doublech=0; % flag for 1ch or 2ch analysis
maskinput=1;  % load n number of masks

%insert loop for multiple files
% if isfile('list.mat')
%     load('list.mat')
% else

% selection of files
[List, Names]=ICCS_SelectFiles(maskinput) ; % description: load Count mask nuclei, Count mask nuclei (again), edu Mask, intensity file 
save('list.mat', 'List' );
save('Files-Analyzed.mat', 'Names' );
% end

%read one file to check number of channels of intensity 
for ifile=1
filenamefull1=List{ifile,1};
MCstack=simple_ICCS_readfiletif(filenamefull1);  % read stack of count mask nuclei
nframes=size(MCstack,3);
if size(List,2)>2 &&~isempty(List{ifile,3})
filenamefulle2 = List{ifile,3};
Astack=simple_ICCS_readfiletif(filenamefulle2);  % order :  EdU, ch2, Dapi if 3 channels
nch=uint8(fix(size(Astack,3)./nframes)); % estimate number of channels
end
for istack=1
MC=MCstack(:,:,istack);
MinN=min(MC(MC>0),[], 'all');
MaxN=max(MC,[], 'all');  
Atest=Astack(:,:,(istack-1)*nch+1:(istack-1)*nch+nch);   

end
end
Nch=size(Atest,3);
ChEdu=1;
ChDNA=Nch;

% menu for channels
prompt = {'EdU channel', 'DNA channel'}; 
dlg_title = 'Press OK to confirm'  ; 
num_lines = 1;
def = { num2str(ChEdu), num2str(ChDNA) };
answer = inputdlg(prompt,dlg_title,num_lines,def);
ChEdu=str2num(answer{1});
ChDNA=str2num(answer{2});


cell=0;
% ResultSpot2=double.empty(0,8);
ResultCell=double.empty(0,11);
% ResultCellG1=double.empty(0,11);
% ResultCellEarly=double.empty(0,11);
% ResultCellLate=double.empty(0,11);
% ResultDist=double.empty(0,7);

% start reading and processing images
for ifile=1:size(List,1)
%open files ... 
% count mask nuclei
filenamefull1=List{ifile,1};
MCstack=simple_ICCS_readfiletif(filenamefull1);  % read stack of count mask nuclei
filenamefulle = List{ifile,2};
MCstack1=simple_ICCS_readfiletif(filenamefulle); % read stack of count mask spots (nuclei in this file!!!)
if size(List,2)>2 &&~isempty(List{ifile,3})
filenamefulle2 = List{ifile,3};
Astack=simple_ICCS_readfiletif(filenamefulle2);  % order :  EdU, ch2, Dapi
end
if maskedu>0 
filenamefulle4 = List{ifile,4};
MaskEduStack=simple_ICCS_readfiletif(filenamefulle4);  % load Mask for analysis of EdU pixels
end

for istack=1:size(MCstack,3)
MC=MCstack(:,:,istack);
MinN=min(MC(MC>0),[], 'all');
MaxN=max(MC,[], 'all');
% count mask CH1
MC1=MCstack1(:,:,istack);
% intensity file
if size(List,2)>2 &&~isempty(List{ifile,3})
A=Astack(:,:,(istack-1)*nch+1:(istack-1)*nch+nch);
end
if maskedu>0 
MaskEdu=MaskEduStack(:,:,istack);
else
MaskEdu=MC;  % 
end

if mask>0 && ~isempty(List{ifile,3+mask})  %mask>0?  &&   OK per ora
filenamefulle4 = List{ifile,3+mask};
MaskROI=simple_ICCS_readfiletif(filenamefulle4);  % load Mask for ROI of analysis
else
MaskROI=MC;  % 
end

X=size(MC,1);
Y=size(MC,2);
MinCell=max([cellfirst, MinN]);
MaxCell=min([cellmax, MaxN]);
for cellidx=MinCell:MaxCell
cell=cell+1;
% mask of single cell (value 1)    
MaskRaw=zeros(X,Y);
MaskRaw(MC==cellidx)=1;
% MaskRaw(MaskROI==0)=0;  % apply mask defining ROI, if loaded
%mask of spots in ch1 (values from k to k+nspots)
Mask1=MC1;
Mask1(MaskRaw==0)=0;
N1max=uint16(max(Mask1,[], 'all'));
%mask of spots in ch2 (values from k to k+nspots)
% Mask2=MC2;
% Mask2(MaskRaw==0)=0;
% N2max=uint16(max(Mask2,[], 'all'));

% order :  EdU, Damage, Dapi
Imgedu=A(:,:,ChEdu);
Imgspot=A(:,:,ChDNA); % not used 
Imgdapi=A(:,:,ChDNA);
Iedu=0;
Ispot=0;
%store cropped image
[row,col]=find(MaskRaw>0);  
Atemp=Imgedu;
Atemp(MaskRaw==0)=0;
Atemp(MaskEdu==0)=0; %visualizes only edu mask pixels
Acrop=Atemp(min(row):max(row),min(col):max(col)); %find crop region and crop
Imgeducrop{ cell } =Acrop ;
% Imgeducrop{ cell } = A(min(row):max(row),min(col):max(col),ChEdu); %find crop region and crop

%store dapi image channel
Atemp2=Imgdapi;
Atemp2(MaskRaw==0)=0;
Acrop2=Atemp2(min(row):max(row),min(col):max(col)); %find crop region and crop
Imgdapicrop{ cell } =Acrop2 ;

% operate on DNA image: must normalize for each cell
Imgdapi(MaskRaw==0)=0;
% Imgdapi=imtranslate(Imgdapi,[-120/45,-40/45]);     %translation
Imgdapinorm=double(Imgdapi)./max(double(Imgdapi),[],'all');  % normalized dna image for each cell
Idapinorm=mean(Imgdapinorm(MaskRaw>0));
Imgdapimask=Imgdapinorm; 
Imgdapimask(Imgdapinorm>=DNAnormTHR)=1; % mask for high DAPI signal 
Imgdapimask(Imgdapinorm<DNAnormTHR)=0;

if celltest==cell
Imgdapinormcrop = Imgdapinorm(min(row):max(row),min(col):max(col),1); %find crop region and crop
Imgdapimaskcrop = Imgdapimask(min(row):max(row),min(col):max(col),1); %find crop region and crop
figure
subplot(1,2,1)
imagesc(Imgdapinormcrop)
axis image
subplot(1,2,2)
imagesc(Imgdapimaskcrop)
axis image
end

Idapinormmask=0;
HCI=0;

Npxdapi = numel(find( MaskRaw>0 )) ;
Idapi=mean(Imgdapi(MaskRaw>0));

Npxedu = numel(find(MaskEdu>0 & MaskRaw>0 )) ;
if Npxedu>0
Iedu=mean(Imgedu(MaskEdu>0 & MaskRaw>0));

Idapinormmask=mean(Imgdapinorm(MaskEdu>0 & MaskRaw>0));  % Idapinormmask = av value of normalized DNA image in mask pixels !!!
NpxeduHC = numel(find(MaskEdu>0 & MaskRaw>0 & Imgdapimask>0 )) ;  % pixels with high DAPI signal 
HCI=NpxeduHC/Npxedu; % HC-index = percentage of pixels with high DAPI signal
end
% if N1>0
% Ispot=mean(Imgspot(Mask1>0));
% end

ResultsCellTemp=cat(2, cell , Idapinorm,  Idapinormmask, HCI, Npxdapi, Idapi, Npxedu, Iedu, ifile, istack, double(cellidx)   )  ;
% ResultsCellTemp=cat(2, cell , N1,  SizeAv1, N2,  SizeAv2, Ncoloc , ifile, double(cellidx)   )  ;

ResultCell=cat(1,ResultCell,ResultsCellTemp) ; % must initialize
% if Npxedu/Npxdapi < 0.005
% ResultCellG1=cat(1,ResultCellG1,ResultsCellTemp) ; % must initialize
% elseif ((Npxedu/Npxdapi>0.005) && (Npxedu/Npxdapi<0.4) && Iedu<(5000/0.4)*(Npxedu/Npxdapi))
%       ResultCellEarly=cat(1,ResultCellEarly,ResultsCellTemp) ; % must initialize  
%     else
%       ResultCellLate=cat(1,ResultCellLate,ResultsCellTemp) ; % must initialize
% end
if celltest==0
close all
end
% end
end
end
end
%end processing/storing


%Results cells at the end
% cell  , cell class,  Idapinormmask, dens edu , Npxdapi, Idapi, Npxedu, Iedu, ifile, istack, double(cellidx) 
% headerfile = 

% SORTING of cells 
%definition of variables
densedu=ResultCell(:,7)./ResultCell(:,5) ; 
ResultCell(:,4)=densedu ; 
Intedu=ResultCell(:,8);
% HCIpar=ResultCell(:,4);
sizenuclav=mean( ResultCell(:,5) ) ; 
maxdens=max( densedu ) ;
% HCIthr=median(HCIpar(HCIpar>0));
sizenucl=ResultCell(:,5);


% plot density edu edu vs I-edu
% figure
% semilogx(Intedu(HCIpar>HCIthr), densedu(HCIpar>HCIthr),'o',Intedu(HCIpar<HCIthr), densedu(HCIpar<HCIthr),'s');
% figure
% semilogx(Intedu, densedu,'o');
% figure; histogram(HCIpar(HCIpar>0),30)

% % SORTING of cells 
% %parameters
% minpx=50;  % minimum pixels for EdU signal
% densfactor=0.5;
% Imid = mean( Intedu(densedu>thrdens) ) ; 
% SizeLimit=1.2;
% IntFactor1=3;
% IntFactor2=3;
figcheck=1;
modesel='p';
figure
%plot dens edu vs int edu
while figcheck==1    
thrdens=densfactor*maxdens; 
mindens=minpx/sizenuclav;          % minimum pixel density for EdU signal
maxsizeearly=SizeLimit*sizenuclav;
Imin = min( Intedu(densedu>mindens) );
Imax = max( Intedu(densedu>mindens) );
SlopeScale=(Imax-Imin)/thrdens;
Ithearly = IntFactor1*Imin ; 
Ithlate = IntFactor2*Imin ; 
% kslopeE =  ( Ithlate - Ithearly  ) / thrdens ;

%conditions SIMPLIFIED
c1=densedu>mindens;
c2=densedu>thrdens;
c2m=densedu>1.0*thrdens;
c3=Intedu< ( Ithearly + kslopeE* SlopeScale * densedu ) ;
c4=Intedu>= ( Ithlate + kslopeL* SlopeScale *densedu ) ;
csize=sizenucl<maxsizeearly ;
cmid=c2m ;
cearly=(c1&~c2&c3 & csize) ; 
clate=(c1&~c2&c4) ;
cearlylarge=(c1&~c2&c3 & ~csize) ; 
% cunclass=c1&~c2 &(~cearly &~c4 )  ;
cunclass = ~ ( cearly | clate | cmid | (~c1) ) ;
% subplot(1,2,2)
% histogram(sizenucl( c1&~c2&c3 )./sizenuclav, 20) ;
% axis square
% subplot(1,2,1)
% order:  G1/G2, Mid, E, L
switch modesel
    case 'p'
subplot();
semilogx(Intedu(~c1), densedu(~c1),'*',Intedu(cmid),densedu(cmid),'s',Intedu(cearly),densedu((cearly)),'+',Intedu(clate),densedu((clate)),'d', Intedu(cunclass),densedu((cunclass)),'o');
axis square
% axis square
% Create legend
legend1 = legend('G1/G2', 'M', 'E', 'L' , 'none');
set(legend1,...
    'Position',[0.164880954030724 0.716269845811147 0.1464285697788 0.165476185934884]);
    case 'f'
        %border cells
if isempty(min(densedu(cearly>0)))==0 && isempty(min(Intedu(clate>0)))==0  && isempty(min(Intedu(cearly>0)))==0 && isempty(max(densedu(c1==0)))==0
celllateb=find(Intedu==min(Intedu(clate>0&cearlylarge==0)) & (clate>0&cearlylarge==0) ,1);
cellearlyb=find(Intedu==max(Intedu(cearly>0)) & cearly>0  ,1);
cellGb=find(densedu==max(densedu(c1==0))  & c1==0 ,1);
cellearlyb0=find(densedu==min(densedu(cearly>0))  & cearly>0 ,1);
cellearlybm=find(densedu==max(densedu(cearly>0)) & cearly>0  ,1);
cellmidb=find(densedu==min(densedu(cmid>0)) & cmid>0  ,1);
cellearlybs=find( sizenucl==max(sizenucl(cearly>0)) & cearly>0 ,1) ;

subplot(4,2,2)
% imagesc(  Imgeducrop {cellGb}  ) ;
imshow(  imfuse(Imgeducrop {cellGb} , Imgdapicrop {cellGb}  ))  ;
title('G1/G2')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,1)
% imagesc(  Imgeducrop {cellearlyb0}  ) ;
imshow(  imfuse(Imgeducrop {cellearlyb0} , Imgdapicrop {cellearlyb0}  ))  ;
title('EARLY S')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,3)
% imagesc(  Imgeducrop {cellearlybm}  ) ;
imshow(  imfuse(Imgeducrop {cellearlybm} , Imgdapicrop {cellearlybm}  ))  ;
title('EARLY S')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,4)
% imagesc(  Imgeducrop {cellmidb}  ) ;
imshow(  imfuse(Imgeducrop {cellmidb} , Imgdapicrop {cellmidb}  ))  ;
title('MID S')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,5)
% imagesc(  Imgeducrop {cellearlyb}  ) ;
imshow(  imfuse(Imgeducrop {cellearlyb} , Imgdapicrop {cellearlyb}  ))  ;
title('EARLY S')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,6)
% imagesc(  Imgeducrop {celllateb}  ) ;
imshow(  imfuse(Imgeducrop {celllateb} , Imgdapicrop {celllateb}  ))  ;
title('LATE S')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,7)
% imagesc(  Imgeducrop {cellearlybs}  ) ;
imshow(  imfuse(Imgeducrop {cellearlybs} , Imgdapicrop {cellearlybs}  ))  ;
title('EARLY S')
% colormap(greenmap);
colorbar('off')
axis image
axis off

subplot(4,2,8)
imagesc( zeros(1)   ) ;
if  isempty( min(sizenucl(cearlylarge>0)) )==0
cellearlylargeb=find( sizenucl==min(sizenucl(cearlylarge>0)) & cearlylarge>0 ,1);
% imagesc(  Imgeducrop {cellearlylargeb}  ) ;
imshow(  imfuse(Imgeducrop {cellearlylargeb} , Imgdapicrop {cellearlylargeb}  ))  ;
end
title('excluded EARLY (large nuclei)')
% colormap(greenmap);
colorbar('off')
axis image
axis off
  

end
end

prompt2 = {'Minimum number of EdU pixels for S phase:', 'Density Threshold for Mid (relative to maximum):','Intensity Threshold for Early (relative to minimum):', 'Intensity Threshold for Late (relative to minimum):','Maximum size for Early (relative to average):', 'file name','Plot(p) or Figure(f)'}; 
dlg_title2 = 'Input parameters'; 
num_lines = 1;
def2 = {num2str(minpx),num2str(densfactor),num2str([IntFactor1 kslopeE]),num2str([IntFactor2 kslopeL]), num2str(SizeLimit),tag, modesel};
answer2 = inputdlg(prompt2,dlg_title2,num_lines,def2);
figcheck=~isempty(answer2); 
if figcheck==1
minpx=str2num(answer2{1});
densfactor=str2num(answer2{2});

IntFactor1dual=str2num(answer2{3});
IntFactor1=IntFactor1dual(1);
kslopeE=IntFactor1dual(2);
IntFactor2dual=str2num(answer2{4});
IntFactor2=IntFactor2dual(1);
kslopeL=IntFactor2dual(2);

SizeLimit=str2num(answer2{5});
tag=answer2{6};
modesel=answer2{7};
end
end



for cell=1:size(ResultCell,1)
    
     
if densedu(cell)==0
       CellTag{ cell } = 'G1-G2';
       ResultCell(cell,2)=0;
else
%     if densedu(cell)>thrdens
    if cmid(cell)==1
 
       CellTag{ cell } = 'MID'; 
       ResultCell(cell,2)=2;
        
    else
       if cearly(cell)==1
       CellTag{ cell } = 'EARLY S'; 
       ResultCell(cell,2)=1;
       elseif clate(cell)==1
       CellTag{ cell } = 'LATE S'; 
       ResultCell(cell,2)=3;
       else
        
       CellTag{ cell } = 'none'; 
       ResultCell(cell,2)=4;
          
       end
    end
end
end

numE=numel(find((ResultCell(:,2)==1))) ;
numM=numel(find((ResultCell(:,2)==2))) ;
numL=numel(find((ResultCell(:,2)==3))) ;
numnone=numel(find((ResultCell(:,2)==4))) ;
numG=numel(find((ResultCell(:,2)==0))) ;

% show and save SORTED cells
% greenmap = zeros(51,3);
% greenmap(:,2)=[0:0.02:1]';
maxplot=3;
titleall={ 'G1-G2', 'EARLY S', 'MID', 'LATE S', 'none' } ;
figure; pie([numG numE, numM, numL, numnone],titleall);
for cellcat=1:5
nplot=min( ceil(sqrt(numel(find((ResultCell(:,2)==(cellcat-1)))))) , maxplot);
nplotall=  ceil(sqrt(numel(find((ResultCell(:,2)==(cellcat-1)))))) ;  

if savefigure==1
i=0;
figure
for cell=1: size(ResultCell,1) 
if strcmp(CellTag{ cell } , titleall{cellcat}  ) 
i=i+1;    
subplot(nplotall,nplotall,i)
imagesc(  Imgeducrop {cell}  ) ;
title( { CellTag{ cell }  ; num2str(ResultCell(cell,9:11))},'FontSize',4  );
colormap(greenmap);
colorbar('off')
axis image
axis off
end
end
fnameoutjpg = [tag,'-',titleall{cellcat},'.png'];
saveas(gcf,fnameoutjpg);
close
end

i=0;
figure
for cell=1: size(ResultCell,1) 
if strcmp(CellTag{ cell } , titleall{cellcat}  ) && i<maxplot*maxplot
i=i+1;    
subplot(nplot,nplot,i)
imagesc(  Imgeducrop {cell}  ) ;
title( { CellTag{ cell }  ; num2str(ResultCell(cell,9:11))},'FontSize',7  );
colormap(greenmap);
colorbar('off')
axis image
axis off
end
end

end

% cell  , cell class,  Idapinormmask, dens edu , Npxdapi, Idapi, Npxedu, Iedu, ifile, istack, double(cellidx) 
% selected: 
% cell  , cell class,  Iedu,  dens edu ,  ifile, istack, double(cellidx) 
ResultSelected (:,1:2) = ResultCell(:,1:2);
ResultSelected (:,3) = ResultCell(:,8);
ResultSelected (:,4) = ResultCell(:,4);
ResultSelected (:,5:7) = ResultCell(:,9:11);

filenameout=[tag,'-Results' ]; 
dlmwrite([filenameout,'.txt'],ResultSelected,'delimiter',';','precision',8);

% filenameout2=[tag,'-ResultsComplete' ]; 
% dlmwrite([filenameout2,'.txt'],ResultCell,'delimiter',';','precision',8);

filenameout1=[tag,'-par' ];
save([filenameout1,'.mat'], 'minpx', 'densfactor', 'IntFactor1', 'IntFactor2', 'Ithearly', 'Ithlate', 'maxdens', 'sizenuclav', 'Imin', 'thrdens', 'SizeLimit', 'kslopeE','kslopeL','Imax'  );


end


%required funtions

function A=simple_ICCS_readfiletif(fname)
info = imfinfo(fname);
nslice = numel(info);
A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end

end

function y=simpleICCS_smooth_simple(M,sm,n)
y=M;
if sm>0
filt = (1/(8+1/sm))*[1 1 1; 1 1/sm 1; 1 1 1]; % sm factor <=1 
    for i=1:n
    y = filter2(filt,y);
    end
end
    
end

function [y,Np, varargout]=simple_PadForICS_fromMask(x1,Extra,Mask)
[m,n,p]=size(x1);
Mask=double(Mask);
%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
end
%% adding average on zeros
for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 || isnan(Mask(i,j))
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function [y,Np, varargout]=simple_PadForICS_sm(x1,Extra,Thr, sm)
[m,n,p]=size(x1);
for kk=1:p
x1s(:,:,kk)=simpleICCS_smooth_simple(x1(:,:,kk),sm,2);
end

%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
ys=y;
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
    ys(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1s(:,:,k);
end

%% adding average on zeros
Mask=ys(:,:,1);
Mask2=ys(:,:,2);
Mask(Mask<=Thr(1)& Mask2<=Thr(2))=0;
Mask(Mask2>Thr(2)| Mask>Thr(1) )=1;

for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function  Output=simple_ICCS_CCFmean(x1,x2)

NumberOfAngles=180;
[X,Y]=size(x1);
%ACF=conv2(x1,x2,'same');
F1=fft2(x1);
F2=fft2(x2);
ACF= F1.*conj(F2);
G=((sum(sum(x1)))*(sum(sum(x2)))/X/Y);
ACF= ifft2(ACF);
ACF= fftshift(ACF)./G-1;

[R, C]=size(ACF);
if iseven(R)
r0=R/2+1;
else
r0=(R+1)/2;
end
if iseven(C)
c0=C/2+1;
% Radius=min(R/2,C/2);
else
c0=(C+1)/2;
% Radius=min((R-1)/2,(C-1)/2);
end
Radius=min(r0-1,c0-1);

if NumberOfAngles==1
    Output=ACF(r0,c0:end);
else
ACF1=flipud(ACF(1:r0-1,c0:end));
ACF2=ACF(r0:end,c0:end);
ProfMat=zeros(NumberOfAngles*2,Radius);

for j=1:2
    if j==1
        y=ACF1';
    else
        y=ACF2;
    end
    
% CALCULATION OF ROTATIONAL MEAN
% Definition of angles
t=(pi/NumberOfAngles/2:pi/NumberOfAngles/2:pi/2);
   
% Matrix
y=y(1:Radius,1:Radius);
% Cycle between the 2nd and 2nd to last angles
[~, y1y]=size(y);

for i=1:NumberOfAngles
   rt=ceil(cos(t(i))*(1:Radius));
   ct=ceil(sin(t(i))*(1:Radius));
   profile=y((rt-1).*y1y+ct);

   if j==1
   ProfMat(NumberOfAngles+i,:)=profile;
   else
%    ProfMat(i,:)=profile;
   ProfMat(i,:)=[profile(2:end),profile(end)];  % excluding the central ACF point (0,0)
   end   
end

end


Output=[double(ACF(r0,c0)) sum(ProfMat)./(2*NumberOfAngles)];
% 
% OrientedProfiles=min_fw_Profile;
% OrientedProfiles(2,:)=max_fw_Profile;
% Angles=[min_th,max_th];

end


end

function bool=iseven(x)

if mod(x,2) == 0
bool=1;
else
bool=0;
end
end

function [param, fval, chisq]=simple_ICCS_Fit_ICS_1Dsingle_dist_wfix(x,y,w0,delta0,Display, title1)
%fixed 

my=min(y);
My=max(y);
fun = @(Param) sum( (( (Param(1)+Param(2).*cosh(2*Param(3).*x/w0^2).*exp(-(((x-0).*(x-0)+Param(3)^2)./(w0*w0) ))) -y ).^2)./(abs(y))  ) ;
[param, chisqpar]=fminsearch( fun,[my My  delta0]);
w0=abs(w0);
param(3)=abs(param(3));
fval=(param(1)+param(2).*cosh(2*param(3).*x/w0^2).*exp(-(((x-0).*(x-0)+param(3)^2)./(w0*w0) )));
chisq=sum(  (  (param(1)+param(2).*cosh(2*param(3).*x/w0^2).*exp(-(((x-0).*(x-0)+param(3)^2)./(w0*w0) ))).^2)./(abs(y))  ) ; 
if Display==1
%     figure;
    plot(x,y,'o')
    hold on
    plot(x, fval, '--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat(title1, '  w=',num2str(w0,2),'   G0=',num2str(param(2),2),' d=',num2str(param(3),2)  ));
end
end
function [param, fval, chisq]=simple_ICCS_Fit_ICS_1Dsingle_dist(x,y,w0,delta0,Display, title1)
%fixed 

my=min(y);
My=max(y);
fun = @(Param) sum( (( (Param(1)+Param(2).*cosh(2*Param(4).*x/Param(3)^2).*exp(-(((x-0).*(x-0)+Param(4)^2)./(Param(3)*Param(3)) ))) -y ).^2)./(abs(y))  ) ;
[param, chisqpar]=fminsearch( fun,[my My w0 delta0]);
param(3)=abs(param(3));
param(4)=abs(param(4));
fval=(param(1)+param(2).*cosh(2*param(4).*x/param(3)^2).*exp(-(((x-0).*(x-0)+param(4)^2)./(param(3)*param(3)) )));
chisq=sum(  (  (param(1)+param(2).*cosh(2*param(4).*x/param(3)^2).*exp(-(((x-0).*(x-0)+param(4)^2)./(param(3)*param(3)) ))).^2)./(abs(y))  ) ; 
if Display==1
%     figure;
    plot(x,y,'o')
    hold on
    plot(x, fval, '--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat(title1, '  w=',num2str(param(3),2),'   G0=',num2str(param(2),2),' d=',num2str(param(4),2)  ));
end
end

function [param, fval, chisq]=simple_ICCS_Fit_ICS_1Dsingle(x,y,w0,Display, title1)
%fixed 

my=min(y);
My=max(y);
fun = @(Param) sum( (( (Param(1)+Param(2).*exp(-((x-0).*(x-0)./(Param(3)*Param(3)) ))) -y ).^2)./(abs(y))  ) ;
[param, chisqpar]=fminsearch( fun,[my My w0]);
param(3)=abs(param(3));
fval=(param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) )));
chisq=sum( (( (param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) ))) -y ).^2)./((param(2)^2)  )) ;

if Display==1
%     figure;
    plot(x,y,'o')
    hold on
    plot(x, fval, '--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat(title1, '  w=',num2str(param(3),2),'   G0=',num2str(param(2),2) ));
end
end

function B=simpleICCS_Threshold(A,thr)
  
if length(thr)==1
B=A;
B(B<=thr)=0;
B(B>thr)=1;
else
B=A;
B(B>thr(2))=0;
B(B<=thr(1))=0;
B(B>0)=1;
end

end

function const=simple_ICCS_Get_const(FWHM1,FWHM2)

FWHMval=2:0.05:5;

D = [0.2175    0.2203    0.2010    0.1863    0.1755    0.1714    0.1549;
    0.2203    0.1972    0.1850    0.1815    0.1737    0.1545    0.1400;
    0.2010    0.1850    0.1788    0.1772    0.1675    0.1525    0.1435;
    0.1863    0.1815    0.1772    0.1671    0.1539    0.1452    0.1347;
    0.1755    0.1737    0.1675    0.1539    0.1512    0.1389    0.1316;
    0.1714    0.1545    0.1525    0.1452    0.1389    0.1372    0.1309;
    0.15489	  0.13998	0.14346	  0.13471	0.13159   0.13092	0.12654] ;

[Xs,Ys] = meshgrid(double(1:7),double(1:7));
[Xq,Yq] = meshgrid(1:0.1:7,1:0.1:7);
Dint = interp2(Xs,Ys,D,Xq,Yq);
if FWHM1<2
    pos1=1;
elseif FWHM1>5
    pos1=61;
else
    [~, pos1]=min(abs((FWHMval-FWHM1)));
end

if FWHM2<2
    pos2=1;
elseif FWHM2>5
    pos2=61;
else
    [~, pos2]=min(abs((FWHMval-FWHM2)));
end
    
const=Dint(pos1,pos2);

end

function shift=ACF2Dshift(CCF,ACF1,sizeCCF)

[X,Y]=size(CCF);
zoom=1+round(0.5*(X-sizeCCF));
CCFzoom=CCF(zoom:X-zoom,zoom:Y-zoom);
ACF1zoom=ACF1(zoom:X-zoom,zoom:Y-zoom);
c0=round(sizeCCF/2)+1;
[X1,Y1]=meshgrid(0:sizeCCF,0:sizeCCF);
[xmax0, ymax0, r0, teta0]=Find_max_ICS(X1,Y1,ACF1zoom,'0');
[xmax, ymax, r, teta]=Find_max_ICS(X1,Y1,CCFzoom,'0');
rshift=sqrt( (xmax-xmax0)^2 + (ymax-ymax0)^2 ) ;
tetashift=angle( (xmax-xmax0) +1i* (ymax-ymax0) ) ;
% imagesc(Zfit)

shift = [ rshift  tetashift ] ;

end


function  ACF=simple_ICCS_CCF(x1,x2)
% 2D ACF

[X,Y]=size(x1);
%ACF=conv2(x1,x2,'same');
F1=fft2(x1);
F2=fft2(x2);
ACF= F1.*conj(F2);
G=((sum(sum(x1)))*(sum(sum(x2)))/X/Y);
ACF= ifft2(ACF);
ACF= fftshift(ACF)./G-1;

end


function [xmax, ymax, r, teta]=Find_max_ICS(X,Y,Z,Display)

M = max(Z,[],'all') ;
[imax,jmax] = find(Z==M) ;

xmax=X(imax,jmax);
ymax=Y(imax,jmax);
r=sqrt( X(imax,jmax)^2 + Y(imax,jmax)^2 ) ;
teta= angle(X(imax,jmax) +1i*Y(imax,jmax) ) ;

if Display==1
    figure;
    subplot(2,2,1)
    imagesc(Z)
  title(strcat('r=',num2str(r,2),' angle=',num2str(teta,2) ));
%    title(strcat('w=',num2str(param(3),2),'  r=',num2str(r,2),' angle=',num2str(teta,2) ));
%     title(strcat('w=',num2str(param(3),2),'   G0=',num2str(param(2),2) ));

end

if Display=='sub'
    imagesc(Z)
   title(strcat(' r=',num2str(r,2),' angle=',num2str(teta,2) ));
%     title(strcat('w=',num2str(param(3),2),'   G0=',num2str(param(2),2) ));

end

end



function [List, Names]=ICCS_SelectFiles(mask) ;

%open files ... 
filename='none';
c=1;
while filename~= 0
[filename,pathname, filterindex] = uigetfile({'*.tif'},['Select Count Mask Nuclei ',num2str(c)]);
if filename~= 0
filenamefull = [pathname, filename];  
List{c,1}=filenamefull;
Names{c,1}=filename(1:8);

% [filenamemask,pathnamemask, filterindex] = uigetfile({'*.tif'},['Select Count Mask Nuclei ',num2str(c)]);
% filenamefull3 = [pathnamemask, filenamemask]; 
% List{c,2}=filenamefull3;
List{c,2}=filenamefull;

[filenamemaske,pathnamemaske, filterindex] = uigetfile({'*.tif'},['Select Intensity file (EdU in ch1) ',num2str(c)]);
if filenamemaske ~= 0
filenamefulle = [pathnamemaske, filenamemaske];   
List{c,3}=filenamefulle ;
end
if mask>0
    for imask=1:mask
    [filenamemaskr,pathnamemaske, filterindex] = uigetfile({'*.tif'},['Select EdU Mask',num2str(imask),' of file ',num2str(c)]);
        if filenamemaskr ~= 0
        filenamefullr = [pathnamemaske, filenamemaskr];   
        List{c,3+imask}=filenamefullr ;
        end
    end
end

c=c+1;
end

end
end
