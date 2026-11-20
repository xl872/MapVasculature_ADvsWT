clc;close all; clear all;
% Vessel segmentation of each registered volume (Fig 2): Gaussian background
% (sigma 10) removed, sigmoid normalization exp(-(V-B)/1000+1) inside the brain
% mask, fast non-local means, then masks of the top 50 % ... 2 % of brain voxels
% in 2 % steps. Writes <name>_resultinrange50.mat (Vout, raw, percent) next to
% each input. First part: AD/; second part: WT/.
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
namelist = dir('./AD/*reg1+orig.HEAD'); % AD half
mask=BrikLoad ('ambmcmask+orig'); 
Vthr=0.98;
Sthr=0.5;
for index=1:length(namelist)
    
data=[namelist(index).folder,'/',namelist(index).name];%The data for segmentation (after reg by AFNI)
V = BrikLoad (data); 

'loaded'
%%
% volumeViewer(V);
%%

mask1=zeros(size(V));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask1=mask1.*mask;

Vout=maxfilt(V,3);
B = imgaussfilt3(Vout,10);
%
V1=exp(-(V-B)/1000+1).*mask1;   % sigmoid normalization, constant 1000
volumeViewer(V1);
% histogram(V1);

%%
V2=FastNonLocalMeans3D(V1,0.2);
raw=V2;
raw=raw(:,1:end-165,:);
% volumeViewer(V2);
%%
mask1=mask1(:,1:end-165,:);
[N,edges] = histcounts(raw(mask1==1));
N1=cumsum(N);
brainsize=sum(mask1(:),'all');

Vout={};
percent=Sthr:0.02:Vthr;
for perct= 1:length(percent)
[~,i]=find(N1>percent(perct)*brainsize,1);
V3=raw;
thr=edges(i)
V3(V3<thr)=0;
V3(V3~=0)=1;
Vout{perct}=V3;
end
save([data(1:end-5),'_resultinrange50.mat'],'Vout','raw',"percent",'-v7.3');


 end
%%
clc;close all; clear all;
% Second part: the same steps for WT/.
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
namelist = dir('./WT/*reg1+orig.HEAD'); % WT half: the registered reg1 volumes
Vthr=0.98;
Sthr=0.5;
for index=1:length(namelist)
    
data=[namelist(index).folder,'/',namelist(index).name];%The data for segmentation (after reg by AFNI)
V = BrikLoad (data); 
mask=BrikLoad ('ambmcmask+orig'); 
'loaded'
%%
% volumeViewer(V);
%%

mask1=zeros(size(V));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask1=mask1.*mask;


Vout=maxfilt(V,3);
B = imgaussfilt3(Vout,10);
%
V1=exp(-(V-B)/1000+1).*mask1;
% volumeViewer(V1);
% histogram(V1);

%%
V2=FastNonLocalMeans3D(V1,0.2);
raw=V2;
raw=raw(:,1:end-165,:);
% volumeViewer(V2);

%%
mask1=mask1(:,1:end-165,:);
[N,edges] = histcounts(raw(mask1==1));
N1=cumsum(N);
brainsize=sum(mask1(:),'all');

Vout={};
percent=Sthr:0.02:Vthr;
for perct= 1:length(percent)
[~,i]=find(N1>percent(perct)*brainsize,1);
V3=raw;
thr=edges(i);
V3(V3<thr)=0;
V3(V3~=0)=1;
Vout{perct}=V3;
end
save([data(1:end-5),'_resultinrange50.mat'],'Vout','raw',"percent",'-v7.3');
end

%%
function Vout=maxfilt(V,r)
Vout=V;
for i =1+r:size(V,1)-r
    for j=1+r:size(V,2)-r
        for k=1+r:size(V,3)-r
            V(i,j,k)=max(V(i-r:i+r,j-r:j+r,k-r:k+r),[],'all');
        end
    end
end
end

%%
function Vout=extend(V,r,r0)
Vout=V;

c1=zeros(r*2+1,r*2+1,r*2+1);
for i=1:r*2+1
    c1(i,i,i)=1;
end
c2=flip(c1,1);
c3=flip(c1,2);
c4=flip(c1,3);
for i =1+r:size(V,1)-r
    for j=1+r:size(V,2)-r
        for k=1+r:size(V,3)-r
            
            tV=V(i-r:i+r,j-r:j+r,k-r:k+r);
            if sum(tV,'all')==0
                continue
            else
                tV=(1+sign(tV-mean(tV,"all")))/2;
                if sign(sum(tV.*c1,'all')-2*r+r0)+sign(sum(tV.*c2,'all')-2*r+r0)+...
                    sign(sum(tV.*c3,'all')-2*r+r0)+sign(sum(tV.*c4,'all')-2*r+r0)+...
                    sign(sum(tV(r+1,r+1,:))-2*r+r0)+sign(sum(tV(r+1,:,r+1))-2*r+r0)+sign(sum(tV(:,r+1,r+1))-2*r+r0)>=1
                V(i,j,k)=max(V(i-r:i+r,j-r:j+r,k-r:k+r),[],'all');
                end
            end
        end
    end
end
end