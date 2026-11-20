clc;close all; clear all;
% Control permutations of the vessel-density maps (Fig S4, Methods): 5 WT + 5 AD
% against another 5 WT + 5 AD, 1000 epochs, tag 1 (top 50 %).
% One indexall=randperm(11) per epoch: group A = AD indexall(1:5) + WT
% indexall(6:10), group B = WT indexall(1:5) + AD indexall(6:10). Each animal's
% 450-um block VD map and its left-right mirror go into consecutive slots
% (20 maps per group); the slot array is cleared before each group. Per epoch:
% voxel-wise ttest2 inside the brain (p = 1 outside) and the difference map
% A - B. The seed is fixed (rng(SEED,'twister')).
% Inputs (WORK_DIR): AD/ and WT/ as for Random_group_result_double_density.m.
% Output (WORK_DIR): CON_Density<tag>nosr.mat (map, MVDV = p < 0.05, draws, SEED),
%   read by findROI.m for the Fig S4 null distributions.
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
%
SEED = 20261118;
PTHR=0.01;
epochs=1000;
mask= BrikLoad ('ambmcmask+orig.BRIK');



mask1=zeros(size(mask));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask=mask1.*mask;
mask=mask(:,1:end-165,:);
sr=15;
masksr=imresize3(mask,1/sr,"Method","nearest");
ANIMAL = struct('AD',{cell(1,11)},'WT',{cell(1,11)});                  % per-animal cache
% tag=11;
for tag=[1]%0:raw, 1~10:0.8~0.98 percent brain voxels, 11 <0.8 voxels, -1~-9:from 0.8~0.96 to 0.98 percent
% 6 11 16 21 25 -6 -11 -16 -21
    %50: 0:raw, 1~25:0.5~0.98 percent brain voxels, 26 <0.5 voxels, -1~-24:from 0.5~0.96 to 0.98 percent
rng(SEED,'twister');
fprintf('[%s] control VD tag %d: rng(%d,''twister''), %d epochs\n', datestr(now), tag, SEED, epochs);


close all;
map=[];MVDV=[];
sel=zeros(epochs,10);                                                 % record the draws
for epoch=1:epochs
    indexall=randperm(11);
    sel(epoch,:)=indexall(1:10);
namelist = dir('./AD/*resultinrange50.mat'); % or '.\AD\*.HEAD'
VD=[];VN=[];VDV=[];ALL=[];
slot=0;
sr=15;%size of average window: sr*sr*sr
for index=indexall(1:5)
[Vdv,~,ANIMAL.AD{index}]=animalVdv(ANIMAL.AD{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
namelist = dir('./WT/*resultinrange50.mat'); % or '.\AD\*.HEAD'
for index=indexall(6:10)
[Vdv,~,ANIMAL.WT{index}]=animalVdv(ANIMAL.WT{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
%

AD.VDV=VDV;

%indexall=randperm(11);
namelist = dir('./WT/*resultinrange50.mat'); % or '.\AD\*.HEAD'
VD=[];VN=[];VDV=[];ALL=[];
slot=0;
for index=indexall(1:5)
[Vdv,~,ANIMAL.WT{index}]=animalVdv(ANIMAL.WT{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
namelist = dir('./AD/*resultinrange50.mat'); % or '.\AD\*.HEAD'
for index=indexall(6:10)
[Vdv,~,ANIMAL.AD{index}]=animalVdv(ANIMAL.AD{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
WT.VDV=VDV;


%%
WT0=WT.VDV;%
AD0=AD.VDV;%
if size(WT0,4)~=20 || size(AD0,4)~=20 || any(all(WT0==0,[1 2 3])) || any(all(AD0==0,[1 2 3])) % guard: 20 non-empty maps per group
    error('epoch %d: expected 20 non-empty maps per group', epoch);
end

pVDV=ones(size(WT0(:,:,:,1)));
[~,pAll] = ttest2(AD0,WT0,'Dim',4);                                    % voxel-wise ttest2;
pVDV(masksr==1)=pAll(masksr==1);                                       %   p = 1 outside masksr

% MVDVn=imresize3((mean(WT0,4)-mean(AD0,4)),sr);
MVDVt=zeros(size(pVDV));
MVDVt(pVDV<0.05)=1;
MVDV(epoch,:,:,:)=MVDVt;

map(epoch,:,:,:)=(mean(WT0,4)-mean(AD0,4));
if mod(epoch,50)==0, fprintf('[%s] epoch %d\n', datestr(now), epoch); end
end
%
save(['CON_Density',num2str(tag),'nosr.mat'],  "map","MVDV",'sel','SEED','-v7.3')   % maps, p < 0.05 masks, draws, seed
if false   % not executed: needs 'mapt', which this script does not compute
s1=size(mapt,2);
s2=size(mapt,3);
s3=size(mapt,4);
parfor(i=1:s1,20)
    for j=1:s2
        for k=1:s3
                [Mh(i,j,k),ph(i,j,k)] = ttest(mapt(:,i,j,k));%
        end
    end
end
save(['CON_Density',num2str(tag),'final.mat'], 'Mh', "ph", '-v7.3')
end
end
function out=myselectVout(mask,Vout,raw,tag,sr)
    out=raw;
    if tag>0&&tag<=length(Vout)
        Tout=Vout{tag};
        out(Tout==0|mask==0)=nan;
        out(~isnan(out))=1;
    elseif tag==0
        out(mask==0)=nan;
    elseif tag<0
        Tout=Vout{-tag};
        out(Vout{min(-tag+5,25)}==1|Tout==0|mask==0)=nan;
        out(~isnan(out))=1;
    end
    out = mysizedownd(out,sr,mask);
    % out = mysizedown(out,sr);
end
function [Vdv,ALLi,C]=animalVdv(C,namelist,index,mask,tag,sr)
% Per-animal part of the epoch loop: loads the animal once and caches its block map in C.
if isempty(C)
data=[namelist(index).folder,'/',namelist(index).name];
load(data,'Vout');
if isempty(strfind(namelist(index).name,'38um'))
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_75um.nii_regraw+orig.BRIK']);
else
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_38um.nii_reg1+orig.BRIK']);
end
raw=raw(:,1:end-165,:);
C.ALL=sum(raw.*mask,"all")/sum(mask,"all");
C.Vdv=myselectVout(mask,Vout,raw,tag,sr);
C.name=namelist(index).name;
fprintf('   loaded %s\n', C.name);
end
if ~strcmp(C.name,namelist(index).name), error('file list changed'); end
Vdv=C.Vdv; ALLi=C.ALL;
end
