clc;close all; clear all;
% Residual-signal (capillary) group map, WT vs AD: Fig 6A (tag 1).
% Each animal's registered *_75um.nii_regraw volume is z-scored inside the brain
% mask; voxels of the top-50 % vessel mask are excluded and the rest is averaged
% in 450-um blocks. 100 permutation epochs: in each epoch, 10 of the 11 animals
% of each group are drawn at random (randperm(11), first 10); each drawn
% animal's block map and its left-right mirror go into consecutive slots (20
% maps per group); the slot array is cleared before each group. Per epoch:
% voxel-wise ttest2 (AD vs WT) and the difference map WT - AD, kept where
% p < 0.05. Over the epochs: one-sample t-test, Benjamini-Hochberg FDR inside
% the brain mask, smoothing of q (sigma = 1 voxel), q < PTHR, and
% |difference| >= 0.01. The seed is fixed (rng(SEED,'twister')).
% Inputs (WORK_DIR): AD/ and WT/ with *_resultinrange50.mat and the registered
%   *_75um.nii_regraw volumes; ambmcmask+orig, ambmc_template_30+orig.
% Outputs (WORK_DIR): Density<tag>WTAD_seg_rawz.mat (group maps of every epoch,
%   draws, SEED), Density<tag>seg_rawz_meanmapt.mat, Density<tag>final_seg_rawz.mat
%   (Mh, ph), Density<tag>seg_rawz_Ah_FDR.mat (significance mask), two preview
%   images and the per-animal block maps (ANIMAL).
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
%
SEED = 20261118;
SAVE_FULL_MAPS = false;                                               % true = also save map/mapt
PTHR=0.01;
ANIMAL = struct('AD',{cell(1,11)},'WT',{cell(1,11)});                  % per-animal cache
% tag=11;
for tag=[1]%0:raw, 1~10:0.8~0.98 percent brain voxels, 11 <0.8 voxels, -1~-9:from 0.8~0.96 to 0.98 percent
%50: 0:raw, 1~25:<0.5~0.98 percent brain voxels indensity
rng(SEED,'twister');
fprintf('[%s] intensity tag %d: rng(%d,''twister'')\n', datestr(now), tag, SEED);
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask1=zeros(size(mask));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask=mask1.*mask;
mask=mask(:,1:end-165,:);
close all;
map=[];
epochs=100;
selAD=zeros(epochs,10); selWT=zeros(epochs,10);                       % record the draws
for epoch=1:epochs
    indexall=randperm(11);
    selAD(epoch,:)=indexall(1:10);
namelist = dir('./AD/*resultinrange50.mat'); % or '.\AD\*.HEAD'
VD=[];VN=[];VDV=[];ALL=[];
slot=0;
sr=15;%size of average window: sr*sr*sr
for index=indexall(1:5)
[Vdv,ALL(index),ANIMAL.AD{index}]=animalVdv(ANIMAL.AD{index},namelist,index,mask,tag,sr); % z-scored block map of this animal (computed once, then cached)
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
namelist = dir('./AD/*resultinrange50.mat'); % or '.\AD\*.HEAD'
for index=indexall(6:10)
[Vdv,ALL(index),ANIMAL.AD{index}]=animalVdv(ANIMAL.AD{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
%

AD.VDV=VDV;AD.ALL=ALL;

indexall=randperm(11);
selWT(epoch,:)=indexall(1:10);
namelist = dir('./WT/*resultinrange50.mat'); % or '.\AD\*.HEAD'
VD=[];VN=[];VDV=[];ALL=[];
slot=0;
for index=indexall(1:5)
[Vdv,ALL(index),ANIMAL.WT{index}]=animalVdv(ANIMAL.WT{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
namelist = dir('./WT/*resultinrange50.mat'); % or '.\AD\*.HEAD'
for index=indexall(6:10)
[Vdv,ALL(index),ANIMAL.WT{index}]=animalVdv(ANIMAL.WT{index},namelist,index,mask,tag,sr);
slot=slot+1;
VDV(:,:,:,slot*2-1)=Vdv;
VDV(:,:,:,slot*2)=flip(Vdv,1);
end
WT.ALL=ALL;
WT.VDV=VDV;


%%
WT0=WT.VDV;%
AD0=AD.VDV;%
WT00{epoch}=WT.VDV;
AD00{epoch}=AD.VDV;
if size(WT0,4)~=20 || size(AD0,4)~=20 || any(all(WT0==0,[1 2 3])) || any(all(AD0==0,[1 2 3])) % guard: 20 non-empty maps per group
    error('epoch %d: expected 20 non-empty maps per group', epoch);
end

[MVDV,pVDV] = ttest2(AD0,WT0,'Dim',4);                                 % voxel-wise ttest2

% MVDVn=imresize3((mean(WT0,4)-mean(AD0,4)),sr);
MVDV=zeros(size(pVDV));
MVDV(pVDV<0.05)=1;
MVDVh=imresize3(MVDV,sr);
% MVDVWT=imresize3(mean(WT0,4),sr);
% MVDVAD=imresize3(mean(AD0,4),sr);
p(epoch,:,:,:)=pVDV;
% mapnotresize(epoch,:,:,:)=(mean(WT0,4)-mean(AD0,4)).*MVDV;
if epoch==1, map=zeros([epochs size(MVDVh)]); mapt=zeros([epochs size(MVDVh)]); end % preallocate
map(epoch,:,:,:)=imresize3((mean(WT0,4)-mean(AD0,4)),sr);
mapt(epoch,:,:,:)=imresize3((mean(WT0,4)-mean(AD0,4)),sr).*imresize3(MVDV,sr);
% maptr(epoch,:,:,:)=imresize3(imresize3((mean(WT0,4)-mean(AD0,4)),sr).*imresize3(MVDV,sr),1/sr);

if mod(epoch,10)==0, fprintf('[%s] epoch %d\n', datestr(now), epoch); end

end
save(['Density',num2str(tag),'WTAD_seg_rawz.mat'], 'WT00', "AD00", 'selWT', 'selAD', 'SEED', '-v7.3')  % group maps of every epoch, draws, seed
mapt_mean=squeeze(mean(mapt,1));
if SAVE_FULL_MAPS
save(['Density',num2str(tag),'seg_rawz.mat'], 'map', "mapt", '-v7.3')
end
save(['Density',num2str(tag),'seg_rawz_meanmapt.mat'], 'mapt_mean', 'p', 'SEED', 'selWT', 'selAD', '-v7.3')
%%
s1=size(mapt,2); s2=size(mapt,3); s3=size(mapt,4);                     % one-sample t-test over the epochs, in slabs of k:
Mh=zeros(s1,s2,s3); ph=zeros(s1,s2,s3);                                % [Mh(i,j,k),ph(i,j,k)] = ttest(mapt(:,i,j,k))
for k0=1:20:s3
    kk=k0:min(k0+19,s3);
    [h_,p_]=ttest(mapt(:,:,:,kk),0,'Dim',1);
    Mh(:,:,kk)=reshape(h_,s1,s2,numel(kk)); ph(:,:,kk)=reshape(p_,s1,s2,numel(kk));
end
clear h_ p_
%%
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask=imresize3(mask,1/sr);
mask=imresize3(mask,sr);
mask=permute(mask,[3,2,1]);
[U] = BrikLoad ('ambmc_template_30+orig.BRIK');
U=permute(U,[3,2,1]);U=imresize3(U,size(mask));
mask=mask(:,1:end-165,:);
shiftslide=20;
bn=3;
lb=0;
climB = [min(U(:)), max(U(:))];clear M;
climF = [-0.5,0.5];
xn=6;yall=18;ly=1;
yn=yall/xn;
step=10;

save(['Density',num2str(tag),'final_seg_rawz.mat'], 'Mh', "ph", '-v7.3')

% Ah=permute(imresize3(Mh,sr),[3,2,1]);
h = figure('Position',[100,100,1050,400]);

Ah=permute(Mh,[3,2,1]);
% significance mask: Benjamini-Hochberg FDR of the p map inside the brain mask,
% q smoothed (sigma = 1 voxel), q < PTHR
ph_c = permute(ph,[3,2,1]);
inb  = mask>0.5 & ~isnan(ph_c);                    % tested voxels
[ps,ix] = sort(ph_c(inb));                         % Benjamini-Hochberg (= pval_adjust 'fdr')
nP = numel(ps);
qs = min(1, cummin(ps.*nP./(1:nP)','reverse'));
qv = zeros(nP,1); qv(ix) = qs;
q_c = ph_c; q_c(inb) = qv;                         % q inside the brain
ph_cg = imgaussfilt3(q_c,1);                       % then smooth (sigma = 1 voxel)
Ah = zeros(size(ph_c));
pthr = PTHR;
Ah(ph_cg<pthr & mask>0.5) = 1;                     % then threshold
clear inb ps ix nP qs qv q_c

dF=squeeze(mean(mapt,1));
dF=permute(dF,[3,2,1]);
Ah(abs(dF)<0.01)=0;
save(['Density',num2str(tag),'seg_rawz_Ah_FDR.mat'], 'Ah', '-v7.3')
i=1;
 tiledlayout(yn*ly+bn*lb,xn,'TileSpacing','none','Padding','tight');

    for j=1:yall
        j=j+shiftslide;
        nexttile([ly 1]);B=squeeze(U(:,j*step,:));
        B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);
        hB = imagesc(B);axis xy
        axis image off;
        hold on;
        F=squeeze(dF(:,j*step,:));
        hF = imagesc(F,climF);axis xy
        colormap('jet')
        alphadata =0.6*squeeze(Ah(:,j*step,:)).*squeeze(mask(:,j*step,:));%0.6*squeeze(mask(:,j*5,:));%
        set(hF,'AlphaData',alphadata);
    end
    %
saveas(gcf,['new',num2str(tag),'WT-AD',num2str(PTHR),'seg_rawz.jpg'],'jpg')
%
h = figure('Position',[100,100,1050,400]);

% dF=MVDVn;%MVDVn MVDVWT MVDVAD
dF=squeeze(mean(mapt,1));
% dF=squeeze(mapt(1,:,:,:));
dF=permute(dF,[3,2,1]);
i=1;
 tiledlayout(yn*ly+bn*lb,xn,'TileSpacing','none','Padding','tight');

    for j=1:yall
        j=j+shiftslide;
        nexttile([ly 1]);B=squeeze(U(:,j*step,:));
        B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);
        hB = imagesc(B);axis xy
        axis image off;
        hold on;
        F=squeeze(dF(:,j*step,:));
        hF = imagesc(F,climF);axis xy
        colormap('jet')
        alphadata =0.6*squeeze(mask(:,j*step,:));%0.6*squeeze(Ah(:,j*5,:));%
        set(hF,'AlphaData',alphadata);
    end
    saveas(gcf,['new',num2str(tag),'WT-ADseg_rawz.jpg'],'jpg')
%%

end
save('rerun_animal_cache_INT.mat','ANIMAL','SEED','-v7.3')             % per-animal block maps (dir order)
%%
function out=myselectVout(mask,Vout,raw,tag,sr)
    out=raw;
    if tag>0&&tag<=length(Vout)
        Tout=Vout{tag};
        out(Tout==1|mask==0)=nan;
%         out(~isnan(out))=1;
%     elseif tag>length(Vout)
%         Tout=Vout{length(Vout)};
%         out(Tout==1|mask==0)=nan;
    elseif tag==0
        out(mask==0)=nan;
    elseif tag<0
        Tout=Vout{-tag};
        out(Vout{-tag+5}==1|Tout==0|mask==0)=nan;
%         out(~isnan(out))=1;
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
