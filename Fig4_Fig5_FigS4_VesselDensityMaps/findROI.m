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
load(['Density',num2str(tag),'final.mat'])
% load(['CON_Density',num2str(tag),'.mat'])
% maptC=mapt;

load(['CON_Density',num2str(tag),'nosr.mat'])
mapC=map;
load(['Density',num2str(tag),'.mat'])
maptR=mapt;

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
j=20+6;

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
cg=[19 20 169 170 171 172 173 174 175 176 177 178];
septal=[472 443 432];
pir=[311 312 313 437 438 439 ];
SS=[74:80];
MO=[18];

% j=20+6;
mediodorsal=[537 536 535];%pt=0.01
cm=[546];%pt=0.01
VISC=88:94;

%j=20+9;
dg=[359 361];
amygdalar=457;
PF=549;
PAG=[556];

%j=20+12;
mgn=[513 514 515];
vta=[654 655 686 ];
len=367:377;
SC=658:665;
PH=[626];
MM=602;

%j=20+17
men=367:380;
IC=646:648;
ENTl=366:377;
SC5=660:664;
PAG5=666;


Rchose=MO;
ROIname='MO';
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
% R(F<0.01)=0;
R(F>-0.01)=0;
% R(F>-0.05)=0;%SC

hB = imagesc(R);axis xy
%%
ROI=zeros(size(U));
ROI(:,j*step,:)=R;
ROI=permute(ROI,[3,2,1]);
save(['ROI',ROIname,'.mat'],'ROI')
%%


close all
ROInamelist=dir(fullfile(WORK_DIR,'ROISeptal.mat'));
for roi=1:length(ROInamelist)
    ROIname=ROInamelist(roi).name(4:end-4);
%ROIname='Amygdalar';
load(['ROI',ROIname,'.mat'])
disR=[];disC=[];disCs=[];
for i=1:100
    tmap=squeeze(maptR(i,:,:,:));
    disR(i,:)=tmap(ROI==1);    
end
%
for i=1:1000
    tmap=imresize3(squeeze(mapC(i,:,:,:)),sr);%.*imresize3(squeeze(MVDV(i,:,:,:)),sr);
    disC(i,:)=tmap(ROI==1);
    tmap=imresize3(squeeze(mapC(i,:,:,:)),sr).*imresize3(squeeze(MVDV(i,:,:,:)),sr);
    disCs(i,:)=tmap(ROI==1);
end
%
% figure(Position=[100 100 300 250]);
% [N,edges] = histcounts(disR(disR~=0),20);
% N=N/sum(N);
% plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=1);
% hold on;
% [N,edges] = histcounts(disC(disC~=0),20);
% N=N/sum(N);
% plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=1);
% legend({'WT-AD','Control'},'Box','off')
% xlabel('Vessel density differential')
% ylabel('Distribution')
% title(ROIname)
% print(gcf,['Compare_distribution',ROIname,'.jpg'],'-djpeg','-r300');
% print(gcf,['Compare_distribution',ROIname,'.epsc'],'-depsc','-r300');

%
figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disC,250);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
yl=ylim;
xlm=xlim;

scatter(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all'),yl(2)/25,250,'filled',"pentagram")
f=fit((edges(1:end-1)+edges(2:end))'/2,N','gauss1');
% plot(f)
xlabel('Vessel density differential')
% ledend('off')
if sum(permute(dF,[3,2,1]).*ROI,'all')>0
    xp=f.b1+f.c1/sqrt(2)*1.96;
else
    xp=f.b1-f.c1/sqrt(2)*1.96;
end
line([xp,xp],yl,'Color','r','LineWidth',1);
ylim(yl)
xl=ceil(max(abs(xp),abs(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all')))/0.01)*0.01;
xlim(xlm)
ylabel('Distribution of control')
title(ROIname)
fontsize(13,'points')
print(gcf,['distribution',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['distribution',ROIname,'.eps'],'-depsc','-r300');


figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disCs,250);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
yl=ylim;

scatter(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all'),yl(2)/25,250,'filled',"pentagram")
f=fit((edges(1:end-1)+edges(2:end))'/2,N','gauss1');
% plot(f)
xlabel('Vessel density differential')
% ledend('off')
if sum(permute(dF,[3,2,1]).*ROI,'all')>0
    xp=f.b1+f.c1/sqrt(2)*1.96;
else
    xp=f.b1-f.c1/sqrt(2)*1.96;
end
line([xp,xp],yl,'Color','r','LineWidth',1);
ylim(yl)
xl=ceil(max(abs(xp),abs(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all')))/0.01)*0.01;
xlim([-xl xl])
ylabel('Distribution of control')
title(ROIname)
fontsize(13,'points')
print(gcf,['distribution_stat',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['distribution_stat',ROIname,'.eps'],'-depsc','-r300');


figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disCs(disCs~=0),50);  % 50 bins
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
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
print(gcf,['distribution_stat_nonzero',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['distribution_stat_nonzero',ROIname,'.eps'],'-depsc','-r300');
end

%%


close all
ROInamelist=dir(fullfile(WORK_DIR,'ROI*.mat'));
for roi=1:length(ROInamelist)
    ROIname=ROInamelist(roi).name(4:end-4);
%ROIname='Amygdalar';
load(['ROI',ROIname,'.mat'])

disR=[];disC=[];disCs=[];

%
for i=1:1000
    tmap=imresize3(squeeze(mapC(i,:,:,:)),sr);%.*imresize3(squeeze(MVDV(i,:,:,:)),sr);
    disC(i,:)=mean(tmap(ROI==1));
    tmap=imresize3(squeeze(mapC(i,:,:,:)),sr).*imresize3(squeeze(MVDV(i,:,:,:)),sr);
    disCs(i,:)=mean(tmap(ROI==1));
end
%

figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disC,250);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
yl=ylim;
xlm=xlim;

scatter(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all'),yl(2)/25,250,'filled',"pentagram")
f=fit((edges(1:end-1)+edges(2:end))'/2,N','gauss1');
plot(f)
xlabel('Vessel density differential')
ledend('off')
if sum(permute(dF,[3,2,1]).*ROI,'all')>0
    xp=f.b1+f.c1/sqrt(2)*1.96;
else
    xp=f.b1-f.c1/sqrt(2)*1.96;
end
line([xp,xp],yl,'Color','r','LineWidth',1);
ylim(yl)
xl=ceil(max(abs(xp),abs(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all')))/0.01)*0.01;
xlim(xlm)
ylabel('Distribution of control')
title(ROIname)
fontsize(13,'points')
print(gcf,['0702test/distribution',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['0702test/distribution',ROIname,'.eps'],'-depsc','-r300');


figure(Position=[100 100 300 250]);
[N,edges] = histcounts(disCs,250);
N=N/sum(N);
plot((edges(1:end-1)+edges(2:end))/2,N,LineWidth=2,Color='k');
hold on
yl=ylim;

scatter(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all'),yl(2)/25,250,'filled',"pentagram")
f=fit((edges(1:end-1)+edges(2:end))'/2,N','gauss1');
plot(f)
xlabel('Vessel density differential')
ledend('off')
if sum(permute(dF,[3,2,1]).*ROI,'all')>0
    xp=f.b1+f.c1/sqrt(2)*1.96;
else
    xp=f.b1-f.c1/sqrt(2)*1.96;
end
line([xp,xp],yl,'Color','r','LineWidth',1);
ylim(yl)
xl=ceil(max(abs(xp),abs(sum(permute(dF,[3,2,1]).*ROI,'all')/sum(ROI,'all')))/0.01)*0.01;
xlim([-xl xl])
ylabel('Distribution of control')
title(ROIname)
fontsize(13,'points')
print(gcf,['0702test/distribution_stat',ROIname,'.jpg'],'-djpeg','-r300');
print(gcf,['0702test/distribution_stat',ROIname,'.eps'],'-depsc','-r300');

end