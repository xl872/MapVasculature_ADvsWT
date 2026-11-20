clc;close all; clear all;
% Fig S1: regional signal before vs after super-resolution, and per-animal linear fit (a, b).
% DATA_DIR must contain AD/ and WT/ with the raw 75 um *.nii and the matching *.sr.nii.gz
% (SR output of ../SuperResolution). data/ holds the values this script computes.
DATA_DIR = './SR_input';
for rs=0.01
mask = BrikLoad ('ambmcmask+orig.BRIK');
anat = BrikLoad ('Allen_label_AMBMC_30+orig.BRIK');
mask1=zeros(size(mask));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask1=mask1.*mask;
mask=mask(:,1:end-165,:);
anat = anat(:,1:end-165,:).*mask;
namelist = dir(fullfile(DATA_DIR,'AD','*.nii')); % or '.\AD\*.HEAD'
bf=[];af=[];bfs=[];afs=[];
for index=1:length(namelist) 
data=[namelist(index).folder,'/',namelist(index).name];%The data for segmentation (after reg by AFNI)
V0=niftiread(data);
V1=niftiread([data(1:end-4),'.sr.nii.gz']);
if size(V0,1)*2~=size(V1,1)
    V0=cat(1,V0,double(mode(V0,[1 2 3]))*ones(size(V1,1)/2-size(V0,1),size(V0,2),size(V0,3)));
end
if size(V0,3)*2~=size(V1,3)
    V0=cat(3,V0,double(mode(V0,[1 2 3]))*ones(size(V0,1),size(V0,2),size(V1,3)/2-size(V0,3)));
end
mask=zeros(size(V0));
mask(V0>1000)=1;
V0=double(V0).*mask;
v1=double(V1).*imresize3(mask,2);
bf(index)=mean(V0,'all');
af(index)=mean(V1,'all');

Vt=imresize3(V0, rs)./imresize3(mask, rs);
bfs(index,:)=Vt(:);
Vt=imresize3(V1, rs/2)./imresize3(mask, rs);
afs(index,:)=Vt(:);
end
AD.bf=bf;AD.af=af;
AD.bfs=bfs;AD.afs=afs;
%
linecolors = jet;
A = linecolors(round(linspace(1,256,size(AD.bfs,1))),:); 
myColors = reshape(repmat(A,size(afs,2),1),size(A,1)*size(afs,2),size(A,2));

BFs=single([AD.bfs(:)]);
AFs=[AD.afs(:)];
figure('Position',[100 100 350 300]);scatter(BFs,AFs,40,myColors,'filled');xlim([0 10000]);ylim([0 1600])
xlabel('Image intensity before super-resolution')
ylabel('Image intensity after super-resolution')
title('AD')
fontsize(12,"points")
exportgraphics(gcf,[num2str(rs),'ADbefore2after.png'],'Resolution',300)
%%
namelist = dir(fullfile(DATA_DIR,'WT','*.nii')); % or '.\AD\*.HEAD'
bf=[];af=[];bfs=[];afs=[];
for index=1:length(namelist) 
data=[namelist(index).folder,'/',namelist(index).name];%The data for segmentation (after reg by AFNI)
V0=niftiread(data);
V1=niftiread([data(1:end-4),'.sr.nii.gz']);
if size(V0,1)*2~=size(V1,1)
    V0=cat(1,V0,double(mode(V0,[1 2 3]))*ones(size(V1,1)/2-size(V0,1),size(V0,2),size(V0,3)));
end
mask=zeros(size(V0));
mask(V0>1000)=1;
V0=double(V0).*mask;
v1=double(V1).*imresize3(mask,2);
bf(index)=mean(V0,'all');
af(index)=mean(V1,'all');
Vt=imresize3(V0, rs)./imresize3(mask, rs);
bfs(index,:)=Vt(:);
Vt=imresize3(V1, rs/2)./imresize3(mask, rs);
afs(index,:)=Vt(:);
end
WT.bf=bf;WT.af=af;
WT.bfs=bfs;WT.afs=afs;

%
linecolors = jet;
A = linecolors(round(linspace(1,256,size(WT.bfs,1))),:); 
myColors = reshape(repmat(A,size(afs,2),1),size(A,1)*size(afs,2),size(A,2));

BFs=single([WT.bfs(:)]);
AFs=[WT.afs(:)];
figure('Position',[100 100 350 300]);scatter(BFs,AFs,40,myColors,'filled');xlim([0 10000]);ylim([0 1600])
xlabel('Image intensity before super-resolution')
ylabel('Image intensity after super-resolution')
title('WT')
fontsize(12,"points")
exportgraphics(gcf,[num2str(rs),'WTbefore2after.png'],'Resolution',300)
end
%%
% BF=zscore([AD.bf,WT.bf]);
% AF=zscore([AD.af,WT.af]);
% figure;scatter(BF,AF,'filled')
%%
a=[];b=[];
for i=1:size(WT.bfs,1)
f = fit(double([WT.bfs(i,:)])', [WT.afs(i,:)]',  'poly1');
a(i)=f.p1;
b(i)=f.p2;
end
WT.a=a;WT.b=b;
a=[];b=[];
for i=1:size(AD.bfs,1)
f = fit(double([AD.bfs(i,:)])', [AD.afs(i,:)]',  'poly1');
a(i)=f.p1;
b(i)=f.p2;
end
AD.a=a;AD.b=b;


E=[std(WT.a)/sqrt(length(WT.a)),std(WT.b)/sqrt(length(WT.a));std(AD.a)/sqrt(length(AD.a)),std(AD.b)/sqrt(length(AD.a))];
Y=[mean(WT.a),mean(WT.b);mean(AD.a),mean(AD.b)];

%%
% figure;barwitherr(E(:,1),Y(:,1))
% legend({'WT','AD'},'Location','northwest')
% ylabel('Correlation coefficient')
% set(gca,'XTickLabel',{'a'})
 rng('default')  % For reproducibility
    x1 = AD.b(:);
    x2 = WT.b(:);
    x = [x1; x2];
    g1 = repmat({'AD'},length(x1),1);
    g2 = repmat({'WT'},length(x2),1);
    g = [g1; g2];
    figure(position=[100 100 300 250]);CategoricalScatterplot(x,g,'MarkerSize',5)
    % ylim([0 0.3])
    %boxplot(x,g)
    [h,p]= ttest2(x1,x2)
    ylabel('b')
    fontsize(12,"points")
    exportgraphics(gcf,[num2str(rs),'before2after_b.png'],'Resolution',300)
%%
[h,p]=ttest2(WT.a,AD.a)
[h,p]=ttest2(WT.b,AD.b)
%%

% figure;scatter(WT.bf,WT.af,'filled')
% figure;scatter(AD.bf,AD.af,'filled')