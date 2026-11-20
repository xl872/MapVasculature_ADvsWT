clc;clear all;close all;
folder='WT';
sn=410;%413
VC=[];VW=[];
vesscurvesu_all=[];
vesscurveraw_all=[];
flag=0;
namelist = dir(['./',folder,'/*resultinrange50.mat']);
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
        
        su = BrikLoad([namelist(index).folder,'/',namelist(index).name(1:end-20),'.BRIK']);%reg
        su=su(:,91:end-165,:);
        suRot=imrotate3(su,-50,[1 1 0],"nearest","loose");
        namenew = dir(['./',folder,'/*',namelist(index).name(16:19),'*um.*nii*+orig.BRIK']);
        raw = BrikLoad ([namenew(1).folder,'/',namenew(1).name]);%not reged
        raw=raw(:,91:end-165,:);
        rawRot=imrotate3(raw,-50,[1 1 0],"nearest","loose");
        
        figure('Position',[100,100,100+size(rawRot,2),100+size(rawRot,1)]);
        axes('Position',[0.1 0.1 0.1+0.0015*size(rawRot,2) 0.1+0.0015*size(rawRot,1)],'Box','on');
        % imagesc(repmat(rawRot(:,:,sn),1,1,3)/max(rawRot(:,:,sn),[],'all'));
        imagesc(volRot(:,:,sn));
        if flag == 0
            [x,y]=ginput(2);
            flag = 1;
        end
        figure('Position',[100,100,100+size(rawRot,2),100+size(rawRot,1)]);
        axes('Position',[0.1 0.1 0.1+0.0015*size(rawRot,2) 0.1+0.0015*size(rawRot,1)],'Box','on');
        imagesc(repmat(rawRot(:,:,sn),1,1,3)/max(rawRot(:,:,sn),[],'all'));hold on
        rectangle('Position',[min(x),min(y),max(x)-min(x),max(y)-min(y)],'EdgeColor','b','LineWidth',2);hold off
        figure('Position',[100,100,100+size(rawRot,2),100+size(rawRot,1)]);
        axes('Position',[0.1 0.1 0.1+0.0015*size(rawRot,2) 0.1+0.0015*size(rawRot,1)],'Box','on');
        imagesc(repmat(suRot(:,:,sn),1,1,3)/max(suRot(:,:,sn),[],'all'));hold on
        rectangle('Position',[min(x),min(y),max(x)-min(x),max(y)-min(y)],'EdgeColor','b','LineWidth',2);hold off

        sliceraw=squeeze(rawRot(min(y):max(y),min(x):max(x),sn));
        slicerawbig=squeeze(rawRot(min(y)-7:max(y)+7,min(x)-7:max(x)+7,sn));
        slicerawbig=(slicerawbig-mean(slicerawbig,'all'))/(max(slicerawbig,[],'all')-mean(slicerawbig,'all'));
        slicesu=squeeze(suRot(min(y):max(y),min(x):max(x),sn));
        slicesubig=squeeze(suRot(min(y)-7:max(y)+7,min(x)-7:max(x)+7,sn));
        slicesubig=(slicesubig-mean(slicesubig,'all'))/(max(slicesubig,[],'all')-mean(slicesubig,'all'));
        slicemask=volRot(min(y):max(y),min(x):max(x),sn);
        figure('Position',[100,100,100+3*size(sliceraw,2),100+3*size(sliceraw,1)]);
        axes('Position',[0.1 0.1 0.1+0.01*size(sliceraw,2) 0.1+0.01*size(sliceraw,1)],'Box','on');
        imagesc(repmat(sliceraw,1,1,3)/max(sliceraw,[],'all'));hold on; 
        imagesc(slicemask,'AlphaData',slicemask*0.5);colormap('autumn')

        %%
        CC=bwconncomp(slicemask);
        numPixels{tag,index} = cellfun(@numel,CC.PixelIdxList);
        % figure;histogram(numPixels{index})

        [biggest,idx] = find(numPixels{tag,index}>=100);
        Vout=slicemask;
        for i=1:length(idx)
        Vout(CC.PixelIdxList{idx(i)}) = 0;
        end
        vesscurveraw=zeros(max(length(CC.PixelIdxList),1),15);
        vesscurvesu=zeros(max(length(CC.PixelIdxList),1),15);
        for vi=1:length(CC.PixelIdxList)
            im=zeros(size(slicemask));
            im(CC.PixelIdxList{vi})=1;
            out = regionprops(im, 'Centroid');
            out.Centroid=out.Centroid+7;
            vesscurveraw(vi,:)=(slicerawbig(round(out.Centroid(2)),round(out.Centroid(1))-7:round(out.Centroid(1))+7)+ ...
                slicerawbig(round(out.Centroid(2))-7:round(out.Centroid(2))+7,round(out.Centroid(1)))')/2;
            vesscurvesu(vi,:)=(slicesubig(round(out.Centroid(2)),round(out.Centroid(1))-7:round(out.Centroid(1))+7)+ ...
                slicesubig(round(out.Centroid(2))-7:round(out.Centroid(2))+7,round(out.Centroid(1)))')/2;
            % [~,~,VW(tag,index,vi),~]=findpeaks(-vesscurve(vi,:),"SortStr","descend",'NPeaks',1);
        end
       VCraw(tag,index,:)=mean(vesscurveraw,1);
       vesscurveraw_all=[vesscurveraw_all;vesscurveraw];
       VCsu(tag,index,:)=mean(vesscurvesu,1);
       vesscurvesu_all=[vesscurvesu_all;vesscurvesu];
       index
    end
end

%%
%load("data/before.mat")
before=vesscurveraw_all;
%load('data/after.mat')
after=vesscurvesu_all;

% figure;plot(vesscurve');
%
figure('Position',[100 100 350 250]);
time_all=1:15;
time_all=(time_all-8)*37.5;
ymean = mean(before,1);
ystd= std(before,1)/sqrt(size(before,1));
plot(time_all,ymean,'r-','linewidth',2)
hold on
h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
set(h,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none');
xlim("tight")
hold on

ymean = mean(after,1);
ystd= std(after,1)/sqrt(size(after,1));
plot(time_all,ymean,'b-','linewidth',2)
hold on
h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
set(h,'FaceColor','b','FaceAlpha',0.2,'EdgeColor','none');
xlim("tight")
leg=legend('75um','','37.5um');set(leg,'Box','off');leg.Location='southeast';
ylabel('Normalized signal intensity (a.u.)')
xlabel('Radial distance(um)')
fontsize(10,'points')

%%
load("data/before.mat")
before=vesscurve;
load('data/after.mat')
after=vesscurve;

figure;plot(vesscurve');

figure;
time_all=1:15;
ymean = mean(before,1);
ystd= std(before,1)/sqrt(size(before,1));
plot(time_all,ymean,'r-','linewidth',2)
hold on
h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
set(h,'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none');
xlim("tight")
hold on

time_all=1:15;
ymean = mean(after,1);
ystd= std(after,1)/sqrt(size(after,1));
plot(time_all,ymean,'b-','linewidth',2)
hold on
h=fill([time_all,time_all(end:-1:1)],[ymean-ystd,ymean(end:-1:1)+ystd(end:-1:1)],'b');
set(h,'FaceColor','b','FaceAlpha',0.2,'EdgeColor','none');
xlim("tight")
leg=legend('75um','','37.5um');set(leg,'Box','off');leg.Location='southeast';
fontsize(10,'points')
hold on
%%

load('data/vsallAD.mat')
time_all=1:15;
ymean=squeeze(mean(VC,2));
figure('Position',[100 100 300 200]);plot(time_all,ymean([1 11 16 21 25],:)','linewidth',2);
ystd=squeeze(std(VC,0,2))/sqrt(size(VC,2));
C=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
index=[1 11 16 21 25];
hold on
for i=1:5
h=fill([time_all,time_all(end:-1:1)],[ymean(index(i),:)-ystd(index(i),:),ymean(index(i),end:-1:1)+ystd(index(i),end:-1:1)],C(i,:));
set(h,'FaceColor',C(i,:),'FaceAlpha',0.2,'EdgeColor','none');
end
xlim("tight")
fontsize(10,'points')
leg=legend('50%','30%','20%','10%','2%');set(leg,'Box','off');leg.Location='southeast';
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