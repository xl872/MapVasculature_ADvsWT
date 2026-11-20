clc; clear all; close all;
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
ROIall={};ROIname={};
namelist = dir('./INTROI*.mat');
for i=1:length(namelist)
    load(namelist(i).name);
    ROIall{i}=ROI;
    ROIname{i}=namelist(i).name(7:end-4);
end
% load ROICg.mat
% 
% load ROISeptal.mat
% ROIall{2}=ROI;

tag=1;
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask1=zeros(size(mask));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask=mask1.*mask;
mask=mask(:,1:end-165,:);
sr=15;%size of average window: sr*sr*sr
%%
namelist = dir('./AD/*resultinrange50.mat'); % or '.\AD\*.HEAD'
VD=[];VN=[];vdv=[];ALL=[];


for index=1:11
data=[namelist(index).folder,'/',namelist(index).name];
load(data);

if isempty(strfind(namelist(index).name,'38um'))
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_75um.nii_regraw+orig.BRIK']);
else
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_38um.nii_reg1+orig.BRIK']);
end
raw=raw(:,1:end-165,:);
raw=(raw-mean(raw(mask==1),'all'))/std(raw(mask==1),0,'all');
VV=raw;
ALL(index)=sum(raw.*mask,"all")/sum(mask,"all");

Vdv=myselectVout(mask,Vout,raw,tag,sr);
VDV(:,:,:,index*2-1)=Vdv;
VDV(:,:,:,index*2)=flip(Vdv,1);
end
AD.VDV=VDV;AD.ALL=ALL;

namelist = dir('./WT/*resultinrange50.mat'); % or '.\AD\*.HEAD'
for index=1:11
data=[namelist(index).folder,'/',namelist(index).name];%The data for segmentation (after reg by AFNI)
load(data);

if isempty(strfind(namelist(index).name,'38um'))
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_75um.nii_regraw+orig.BRIK']);
else
    raw=BrikLoad([namelist(index).folder,'/',namelist(index).name(1:19),'_38um.nii_reg1+orig.BRIK']);
end
raw=raw(:,1:end-165,:);
raw=(raw-mean(raw(mask==1),'all'))/std(raw(mask==1),0,'all');
VV=raw;

Vdv=myselectVout(mask,Vout,raw,tag,sr);
VDV(:,:,:,index*2-1)=Vdv;
VDV(:,:,:,index*2)=flip(Vdv,1);

ALL(index)=sum(raw.*mask,"all")/sum(mask,"all");
end
%

WT.VDV=VDV;WT.ALL=ALL;
save forROIbarplotINT WT AD
%%
% load forROIbarplotINT.mat
MdenWT=[];MdenAD=[];
for roi=1:length(ROIall)
    ROI=ROIall{roi};
    denWT=[];denAD=[];
    for i=1:22
        tVDV=imresize3(WT.VDV(:,:,:,i),sr);
        denWT(i,:)=tVDV(ROI==1);
        tVDV=imresize3(AD.VDV(:,:,:,i),sr);
        denAD(i,:)=tVDV(ROI==1);
    end
    MdenWT(:,roi)=mean(denWT,2);
    MdenAD(:,roi)=mean(denAD,2);
end
%%
for i=1:length(ROIall)
figure(Position=[100 100 150 250]);
h = boxplot([MdenWT(:,i),MdenAD(:,i)],'Labels',{'WT','AD'}, ...
    'colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980]);
set(h,'linewidth',1)
xt = get(gca, 'XTick');
yl=ylim;
yl=yl(1)+0.95*(yl(2)-yl(1));
[h,p]=ttest2(MdenWT(:,i),MdenAD(:,i));
hold on
if p<0.001
    text(mean(xt([1 2])), yl, '***',HorizontalAlignment='center')
elseif p<0.01
    text(mean(xt([1 2])), yl, '**',HorizontalAlignment='center')
elseif p<0.05
    text(mean(xt([1 2])), yl, '*',HorizontalAlignment='center')
end
ylabel('Normalized image intensity (a.u.)')
% title(ROIname{i})

fontsize(10,"points")
print(gcf,['ROIbar',ROIname{i},'.jpg'],'-djpeg','-r300');
print(gcf,['ROIbar',ROIname{i},'.eps'],'-deps','-r300');
end
%%
figure(Position=[100 100 1000 250]);
h = boxplotGroup({MdenWT,MdenAD}, 'PrimaryLabels', {'WT','AD'}, ...
    'SecondaryLabels', ROIname, ...
    'groupLabelType', 'vertical',  ...
    'Colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980], 'GroupType', 'betweenGroups');
ylabel('Vessel density')
print(gcf,['ROIbar.jpg'],'-djpeg','-r300');
print(gcf,['ROIbar.eps'],'-deps','-r300');
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