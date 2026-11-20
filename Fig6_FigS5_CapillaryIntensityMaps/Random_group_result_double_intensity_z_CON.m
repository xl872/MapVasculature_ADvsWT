clc;close all; clear all;
% Control permutations of the residual-signal maps (Fig S5, Methods): 5 WT + 5 AD
% against another 5 WT + 5 AD, 100 epochs, tag 1.
% One indexall=randperm(11) per epoch: group A = AD indexall(1:5) + WT
% indexall(6:10), group B = WT indexall(1:5) + AD indexall(6:10). Each animal's
% z-scored 450-um block map (top-50 % vessel voxels excluded) and its
% left-right mirror go into consecutive slots (20 maps per group); the slot
% array is cleared before each group. Per epoch: voxel-wise ttest2 (h and p)
% and the difference map A - B. The seed is fixed (rng(SEED,'twister')).
% Inputs (WORK_DIR): AD/ and WT/ as for Random_group_result_double_intensity_z.m.
% Output (WORK_DIR): CON_Density<tag>seg_rawznosr.mat (map, MVDV = ttest2 h, p,
%   draws, SEED). The Fig S5 null distributions use
%   mapt = imresize3(map,15).*imresize3(MVDV,15), with MVDV = 0 where it is NaN
%   (outside the brain).
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
%
SEED = 20261118;
PTHR=0.01;
ANIMAL = struct('AD',{cell(1,11)},'WT',{cell(1,11)});                  % per-animal cache
% tag=11;
for tag=[1]%0:raw, 1~10:0.8~0.98 percent brain voxels, 11 <0.8 voxels, -1~-9:from 0.8~0.96 to 0.98 percent
%50: 0:raw, 1~25:<0.5~0.98 percent brain voxels indensity
rng(SEED,'twister');
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask1=zeros(size(mask));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask=mask1.*mask;
mask=mask(:,1:end-165,:);
close all;
map=[];
epochs=100;                                                           % control epochs
fprintf('[%s] control intensity tag %d: rng(%d,''twister''), %d epochs\n', datestr(now), tag, SEED, epochs);
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

[hh,pVDV] = ttest2(AD0,WT0,'Dim',4);                                   % voxel-wise:
MVDV(epoch,:,:,:)=hh;                                                  %   [MVDV(epoch,i,j,k),pVDV(i,j,k)] = ttest2(...)

p(epoch,:,:,:)=pVDV;
map(epoch,:,:,:)=mean(WT0,4)-mean(AD0,4);
if mod(epoch,10)==0, fprintf('[%s] epoch %d\n', datestr(now), epoch); end

end
save(['CON_Density',num2str(tag),'seg_rawznosr.mat'], 'map', "MVDV","p",'sel','SEED', '-v7.3')  % maps, h, p, draws, seed
if false   % not executed: needs 'mapt', which this script does not compute
for i=1:size(mapt,2)
    for j=1:size(mapt,3)
        for k=1:size(mapt,4)
            [Mh(i,j,k),ph(i,j,k)] = ttest(mapt(:,i,j,k));%
        end
    end
end
save(['CON_Density',num2str(tag),'final_seg_rawz.mat'], 'Mh', "ph", '-v7.3')
end

end
%%
function out=myselectVout(mask,Vout,raw,tag,sr)
    out=raw;
    if tag>0&&tag<=length(Vout)
        Tout=Vout{tag};
        out(Tout==1|mask==0)=nan;
    elseif tag==0
        out(mask==0)=nan;
    elseif tag<0
        Tout=Vout{-tag};
        out(Vout{-tag+5}==1|Tout==0|mask==0)=nan;
    end
    out = mysizedown(out,sr);
end
function [Vdv,ALLi,C]=animalVdv(C,namelist,index,mask,tag,sr)
% Per-animal part of the epoch loop: loads and z-scores the animal once and
% caches its block map in C.
if isempty(C)
data=[namelist(index).folder,'/',namelist(index).name];
load(data,'Vout');
if isempty(strfind(namelist(index).name,'38um'))
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_75um.nii_regraw+orig.BRIK']);
else
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_38um.nii_reg1+orig.BRIK']);
end
raw=raw(:,1:end-165,:);
raw=(raw-mean(raw(mask==1),'all'))/std(raw(mask==1),0,'all');
C.ALL=sum(raw.*mask,"all")/sum(mask,"all");
C.Vdv=myselectVout(mask,Vout,raw,tag,sr);
C.name=namelist(index).name;
fprintf('   loaded %s\n', C.name);
end
if ~strcmp(C.name,namelist(index).name), error('file list changed'); end
Vdv=C.Vdv; ALLi=C.ALL;
end
