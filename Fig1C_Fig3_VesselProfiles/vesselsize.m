clc;clear all;close all;
% data/vsall<WT|AD>.mat (Fig 3) are made by fig3_vsall.m.
folder='AD';
sn=410;%413
VC=[];VW=[];
vesscurve_all=[];
namelist = dir(['./',folder,'/*resultinrange50.mat']);
mask = BrikLoad ('ambmcmask+orig.BRIK');
anat = BrikLoad ('Allen_label_AMBMC_30+orig.BRIK');
mask1=zeros(size(mask));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask1=mask1.*mask;
anat = anat.*mask;

for tag=16%1:25
    for index=1:7%length(namelist)
        close all;
        data=[namelist(index).folder,'/',namelist(index).name];
        load(data);
        %0:raw, 1~25:0.5~0.98 percent brain voxels
        dF=Vout{tag};
        dF=dF(:,91:end,:);
        %%
        volRot = imrotate3(dF,-50,[1 1 0],"nearest","loose");
        slice = squeeze(volRot(:,131:end,sn));
        figure;imagesc(slice)
        %%
        CC=bwconncomp(slice);
        numPixels{tag,index} = cellfun(@numel,CC.PixelIdxList);
        % figure;histogram(numPixels{index})
        
        [biggest,idx] = find(numPixels{tag,index}>=100);
        Vout=slice;
        for i=1:length(idx)
        Vout(CC.PixelIdxList{idx(i)}) = 0;
        end
        raw = BrikLoad([namelist(index).folder,'/',namelist(index).name(1:end-20),'.BRIK']);%reg
        af_ROI(index,1)=mean(raw(anat>=7&anat<=276));%cortex
        af_ROI(index,2)=mean(raw(anat>=428&anat<=493));%cerebral nueclei
        af_ROI(index,3)=mean(raw(anat>=494&anat<=694));%Subcortical
        af_ROI(index,4)=mean(raw(anat>=695&anat<=806));%Brain stem
        % 
        namenew = dir(['./',folder,'/*',namelist(index).name(16:19),'*um.*nii*+orig.BRIK']);
        raw = BrikLoad ([namenew(1).folder,'/',namenew(1).name]);%not reged
        
        bf_ROI(index,1)=mean(raw(anat>=7&anat<=276));%cortex
        bf_ROI(index,2)=mean(raw(anat>=428&anat<=493));%cerebral nueclei
        bf_ROI(index,3)=mean(raw(anat>=494&anat<=694));%Subcortical
        bf_ROI(index,4)=mean(raw(anat>=695&anat<=806));%Brain stem
        
        raw=raw(:,91:end-165,:);
        rawRot=imrotate3(raw,-50,[1 1 0],"nearest","loose");
        sliceraw=squeeze(rawRot(:,131:end,sn));
        figure;imagesc(repmat(rawRot(:,:,sn),1,1,3)/max(rawRot(:,:,sn),[],'all'))
        bg=sliceraw(221:280,30:65);
        ug=Vout(221:280,30:65);
        figure;imagesc(repmat(bg,1,1,3)/max(bg,[],'all'));hold on; imagesc(ug,'AlphaData',0.3*ug);

        vesscurve=zeros(max(length(CC.PixelIdxList),1),15);
        for vi=1:length(CC.PixelIdxList)
            im=zeros(size(slice));
            im(CC.PixelIdxList{vi})=1;
            out = regionprops(im, 'Centroid');
            vesscurve(vi,:)=(sliceraw(round(out.Centroid(2)),round(out.Centroid(1))-7:round(out.Centroid(1))+7)+ ...
                sliceraw(round(out.Centroid(2))-7:round(out.Centroid(2))+7,round(out.Centroid(1)))')/2;
            [~,~,VW(tag,index,vi),~]=findpeaks(-vesscurve(vi,:),"SortStr","descend",'NPeaks',1);
        end
       VC(tag,index,:)=mean(vesscurve,1);
       vesscurve_all=[vesscurve_all;vesscurve];

       

    end
end



%%
load("data/before.mat")
before=vesscurve;
load('data/after.mat')
after=vesscurve;

% load("data/beforeall.mat")
% before=vesscurve_all;
% load('data/afterall.mat')
% after=vesscurve_all;
% vesscurve=vesscurve_all;

% figure;plot(vesscurve');
before=(before-repmat(mean(before,2),[1,15]))./repmat(mean(before,2),[1,15]);
after=(after-repmat(mean(after,2),[1,15]))./repmat(mean(after,2),[1,15]);

% before=mapminmax(before,-1,1);
% after=mapminmax(after,-1,1);
figure('Position',[100 100 300 250]);
time_all=1:15;
time_all=(time_all-8)*37.5;
[ymean,setting] = mapminmax(mean(before,1));
before=mapminmax(before,setting);
ystd= std(before,1)/sqrt(size(before,1));
plot(time_all,ymean,'r-','linewidth',2)
hold on
h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
set(h,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none');
xlim("tight")
hold on

time_all=1:15;
time_all=(time_all-8)*37.5;
[ymean,setting] = mapminmax(mean(after,1));
after=mapminmax(after,setting);
ystd= std(after,1)/sqrt(size(after,1));
plot(time_all,ymean,'b-','linewidth',2)
hold on
h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
set(h,'FaceColor','b','FaceAlpha',0.2,'EdgeColor','none');
xlim("tight")
leg=legend('75um','','37.5um');set(leg,'Box','off');leg.Location='southeast';
fontsize(10,'points')
hold off
xlabel('Radial distance (um)')
ylabel('Normalized signal intensity (a.u.)')
print(gcf,['vessExp.jpg'],'-djpeg','-r300');
print(gcf,['vessExp.eps'],'-deps','-r300');
%%
gene='WT';
load(['data/vsall',gene,'.mat'])
time_all=1:15;time_all=(time_all-8)*37.5;
ymean=squeeze(mean(VC,2));
for i=1:size(VC,1)
    [ymean(i,:),setting]=mapminmax(ymean(i,:),-1,1);
    VC(i,:,:)=mapminmax(VC(i,:,:),setting);
end

figure('Position',[100 100 350 250]);plot(time_all,ymean([1 11 16 21 25],:)','linewidth',2);
ystd=squeeze(std(VC,0,2))/sqrt(size(VC,2));
C=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
index=[1 11 16 21 25];
hold on
for i=1:5
h=fill([time_all,time_all(end:-1:1)],[ymean(index(i),:)-ystd(index(i),:),ymean(index(i),end:-1:1)+ystd(index(i),end:-1:1)],C(i,:));
set(h,'FaceColor',C(i,:),'FaceAlpha',0.2,'EdgeColor','none');
end
xlim("tight")
xlabel('Radial distance (um)')
ylabel('Normalized signal intensity (a.u.)')
title(gene)
fontsize(11,'points')

leg=legend('50%','30%','20%','10%','2%');set(leg,'Box','off');leg.Location='southeast';
print(gcf,['Vesselpercent',gene,'normalized.jpg'],'-djpeg','-r300');
saveas(gcf,['Vesselpercent',gene,'normalized.eps'],'epsc');

%%
gene='WT';
load(['data/vsall',gene,'.mat'])
time_all=1:15;time_all=(time_all-8)*37.5;
ymean=squeeze(mean(VC,2));
figure('Position',[100 100 350 250]);plot(time_all,ymean([1 11 16 21 25],:)','linewidth',2);
ystd=squeeze(std(VC,0,2))/sqrt(size(VC,2));
C=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
index=[1 11 16 21 25];
hold on
for i=1:5
h=fill([time_all,time_all(end:-1:1)],[ymean(index(i),:)-ystd(index(i),:),ymean(index(i),end:-1:1)+ystd(index(i),end:-1:1)],C(i,:));
set(h,'FaceColor',C(i,:),'FaceAlpha',0.2,'EdgeColor','none');
end
xlim("tight")
xlabel('Radial distance (um)')
ylabel('Signal intensity (a.u.)')
fontsize(11,'points')
leg=legend('50%','30%','20%','10%','2%');set(leg,'Box','off');leg.Location='southeast';
print(gcf,['Vesselpercent',gene,'.jpg'],'-djpeg','-r300');
print(gcf,['Vesselpercent',gene,'.eps'],'-depsc');
%%
VW(VW==0)=nan;
vsm=mean(VW,[2 3],"omitnan");

figure('Position',[100 100 300 200]);plot(50:-2:2,vsm*0.038*0.038,'LineWidth',2)
xlabel('Vessel volume ratio (%)')
ylabel('Mean vessel size (mm^2)')
fontsize(10,"points")
xlim("tight")
% figure;histogram(vs(vs<100));
%

figure;
distributionPlot(reshape(VW,size(VW,1),size(VW,2)*size(VW,3))','histOpt',2); % histOpt=2 works better for uniform distributions than the default
set(gca,'XTick',1:25)
set(gca,'XTickLabel',split(num2str(50:-2:2))')
% ylim([0 0.1])
ylabel('Vessel size (mm^2)','FontSize',12)
xlabel('Vessel volume ratio (%)')
fontsize(10,"points")