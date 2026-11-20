clc;clear all; close all;
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
%%
tag=1;

sr=15;
mask= BrikLoad ('ambmcmask+orig.BRIK');
mask=imresize3(mask,1/sr);
mask=imresize3(mask,sr);
mask=permute(mask,[3,2,1]);
[U] = BrikLoad ('Allen_label_AMBMC30+orig.BRIK');
U=permute(U,[3,2,1]);U=imresize3(U,size(mask),"Method","nearest");
U=U(:,1:end-165,:);

mask=mask(:,1:end-165,:);
shiftslide=20;
bn=3;
lb=0;
climB = [min(U(:)), max(U(:))];clear M;
climF = [-0.1, 0.1];
xn=6;yall=18;ly=1;
yn=yall/xn;
step=10;
%%
load(['Density',num2str(tag),'final_seg_rawz.mat'])
load(['CON_Density',num2str(tag),'seg_rawz.mat'])
mapC=map;maptC=mapt;%pC=p;
% load(['CON_Density',num2str(tag),'seg_rawznosr.mat'])
% mapC=map;
load(['Density',num2str(tag),'seg_rawz.mat'])
mapR=map;maptR=mapt;%pR=p;

%%
PTHR=0.01;

Ah=permute(Mh,[3,2,1]);
ph_c=permute(ph,[3,2,1]);
ph_cg=imgaussfilt3(ph_c,1);
Ah=zeros(size(ph_c));
pthr=PTHR;
Ah(ph_cg<pthr)=1;
% Ah = imgaussfilt3(Ah,1);

% dF=imresize3(squeeze(mean(mapnotresize,1)),sr);
dF=squeeze(mean(maptR,1));
% dF=squeeze(mapt(1,:,:,:));
dF=permute(dF,[3,2,1]);
Ah(abs(dF)<0.005)=0;

%%
j=20+12;

%%
figure;
B=squeeze(U(:,j*step,:));
hB = imagesc(B);axis xy
axis image off;

%%
figure;
B=squeeze(U(:,j*step,:));
BU=B;
B(:)=0;
% j=20+1;
hip=342:361;
vis=147:152;
Rchose=vis;
ROIname='VISp';
for roi=Rchose
    B(BU==roi)=1;
end
% B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);
hB = imagesc(B);axis xy
axis image off;
hold on;
F=squeeze(dF(:,j*step,:));
hF = imagesc(F,climF);axis xy
colormap('jet')
alphadata =0.6*squeeze(Ah(:,j*step,:)).*squeeze(mask(:,j*step,:));%0.6*squeeze(mask(:,j*5,:));% 
set(hF,'AlphaData',alphadata);

%%
figure;
R=B.*squeeze(Ah(:,j*step,:));
R(F<0.01)=0;
% R(F>-0.01)=0;

hB = imagesc(R);axis xy
%%
ROI=zeros(size(U));
ROI(:,j*step,:)=R;
ROI=permute(ROI,[3,2,1]);
save(['INTROI',ROIname,'.mat'],'ROI')

%%
close all
ROInamelist=dir(fullfile(WORK_DIR,'INTROIHi*.mat'));
for roi=1:length(ROInamelist)
    ROIname=ROInamelist(roi).name(7:end-4);
%ROIname='Amygdalar';
load(['INTROI',ROIname,'.mat'])

disR=[];disC=[];
for i=1:100
    tmap=squeeze(maptR(i,:,:,:));
    disR(i,:)=tmap(ROI==1);
end

for i=1:100
    % tmap=imresize3(squeeze(mapC(i,:,:,:)),sr).*imresize3(squeeze(MVDV(i,:,:,:)),sr);
    tmap=squeeze(maptC(i,:,:,:));
    disC(i,:)=tmap(ROI==1);
end
%
figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disR(disR~=0),20);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=1);
hold on;
[N,edges] = histcounts(disC(disC~=0),20);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=1);
legend({'WT-AD','Control'},'Box','off')
xlabel('Vessel density differential')
ylabel('Distribution')
title(ROIname)
print(gcf,['Compare_distribution',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['Compare_distribution',ROIname,'.epsc'],'-depsc','-r300');

%
figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disC(disC~=0),50);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
yl=ylim;

scatter(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all'),yl(2)/25,250,'filled',"pentagram")
f=fit((edges(1:end-1)+edges(2:end))'/2,N','gauss1');
plot(f)
xlabel('Vessel density differential')

if sum(permute(dF,[3,2,1]).*ROI,'all')>0
    xp=f.b1+f.c1/sqrt(2)*1.96;
else
    xp=f.b1-f.c1/sqrt(2)*1.96;
end
yl=ylim;
line([xp,xp],yl,'Color','r','LineWidth',1);
ylim(yl)
xl=ceil(max(abs(xp),abs(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all')))/0.01)*0.01;
xlim([-xl xl])
ylabel('Distribution of control')
title(ROIname)
fontsize(13,'points')
legend('off')
print(gcf,['distribution',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['distribution',ROIname,'.eps'],'-depsc','-r300');
end
%%
%%
close all
ROInamelist=dir(fullfile(WORK_DIR,'INTROI*.mat'));
for roi=1:length(ROInamelist)
    ROIname=ROInamelist(roi).name(7:end-4);
%ROIname='Amygdalar';
load(['INTROI',ROIname,'.mat'])

disR=[];disC=[];


for i=1:100
    % tmap=imresize3(squeeze(mapC(i,:,:,:)),sr).*imresize3(squeeze(MVDV(i,:,:,:)),sr);
    tmap=squeeze(maptC(i,:,:,:));
    disC(i,:)=mean(tmap(ROI==1));
end


%
figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disC,10);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
yl=ylim;

scatter(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all'),yl(2)/25,250,'filled',"pentagram")
f=fit((edges(1:end-1)+edges(2:end))'/2,N','gauss1');
plot(f)
xlabel('Vessel density differential')

if sum(permute(dF,[3,2,1]).*ROI,'all')>0
    xp=f.b1+f.c1/sqrt(2)*1.96;
else
    xp=f.b1-f.c1/sqrt(2)*1.96;
end
yl=ylim;
line([xp,xp],yl,'Color','r','LineWidth',1);
ylim(yl)
xl=ceil(max(abs(xp),abs(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all')))/0.01)*0.01;
xlim([-xl xl])
ylabel('Distribution of control')
title(ROIname)
fontsize(13,'points')
legend('off')
print(gcf,['mean_INTdistribution',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['mean_INTdistribution',ROIname,'.eps'],'-depsc','-r300');
end