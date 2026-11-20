clc;clear all;close all;
load data/AD_ROI.mat
AD_AF=af_ROI;
AD_BF=bf_ROI;
load data/WT_ROI.mat
WT_AF=af_ROI;
WT_BF=bf_ROI;
%
figure(Position=[100 100 400 200]);
y = [mean(WT_BF,1);mean(AD_BF,1)]';         % random y values (3 groups of 4 parameters) 
errY = [std(WT_BF,1);std(AD_BF,1)]'/49;% 10% error
% h = barwitherr(errY, y);% Plot with errorbars
% set(gca,'XTickLabel',{'cortex','cerebral nueclei','Subcortical','Brain stem'})
% lgd=legend('WT-raw','AD-raw');
% lgd.NumColumns = 2;
% legend('boxoff')
h = boxplotGroup({WT_BF, AD_BF}, 'PrimaryLabels', {'WT','AD'}, ...
    'SecondaryLabels', {'cortex','cerebral nueclei','Subcortical','Brain stem'}, ...
    'groupLabelType', 'vertical',  ...
    'Colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980], 'GroupType', 'betweenGroups');

ylabel('Averaged signal intensity')
% set(h(1),'FaceColor',"#0072BD");
% set(h(2),'FaceColor',"#4DBEEE");
%set(h(3),'FaceColor',"#A2142F");
%set(h(4),'FaceColor',"#D95319");
fontsize(10,"points")
print(gcf,'ROIbox.jpg','-djpeg','-r300');
print(gcf,'ROIbox.eps','-depsc','-r300');
%%
load data/vsallAD.mat
VW(VW==0)=nan;
ADVW=reshape(VW,[size(VW,1),size(VW,2)*size(VW,3)]);
ADvsm=mean(VW,[2 3],"omitnan");
ADvsstd=std(VW,0,[2 3],"omitnan")/sqrt((size(VW,2)*size(VW,3)));
load data/vsallWT.mat
index=[1 11 16 21 25];
VW(VW==0)=nan;
WTVW=reshape(VW,[size(VW,1),size(VW,2)*size(VW,3)]);
WTvsm=mean(VW,[2 3],"omitnan");
WTvsstd=std(VW,0,[2 3],"omitnan")/sqrt((size(VW,2)*size(VW,3)));

figure(Position=[100 100 400 200]);
% y = [WTvsm(index),ADvsm(index)];         % random y values (3 groups of 4 parameters) 
% errY = [WTvsstd(index),ADvsstd(index)];% 10% error
% h = barwitherr(errY, y);% Plot with errorbars
% set(gca,'XTickLabel',{'50%','30%','20%','10%','2%'})
% lgd=legend('WT','AD');
% lgd.NumColumns = 2;
% lgd.Location = "northwest";
% legend('boxoff')
% set(h(1),'FaceColor',"#0072BD");
% set(h(2),'FaceColor',"#A2142F");

h = boxplotGroup({WTVW(index,:)', ADVW(index,:)'}, 'PrimaryLabels', {'',''}, ...
    'SecondaryLabels', {'50%','30%','20%','10%','2%'}, ...
    'groupLabelType', 'horizontal', 'groupLines', true, ...
    'Colors', [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980], 'GroupType', 'betweenGroups');
xl=xlim;
ylabel('Vessel size (mm^2)')
xlabel('Vessel volume ratio',Position=[mean(xl),-1])

fontsize(gca,10,"points")
print(gcf,'Percentbox.jpg','-djpeg','-r300');
print(gcf,'PercentROIbox.eps','-depsc','-r300');
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

% figure('Position',[100 100 300 200]);plot(50:-2:2,vsm*0.038*0.038,'LineWidth',2)
% xlabel('Vessel volume ratio (%)')
% ylabel('Mean vessel size (mm^2)')
% fontsize(10,"points")
% xlim("tight")
% figure;histogram(vs(vs<100));
%

figure(Position=[100 100 600 300]);
distributionPlot(reshape(VW,size(VW,1),size(VW,2)*size(VW,3))','histOpt',2); % histOpt=2 works better for uniform distributions than the default
set(gca,'XTick',1:25)
set(gca,'XTickLabel',split(num2str(50:-2:2))')
% ylim([0 0.1])
ylabel('Vessel size (mm^2)','FontSize',12)
xlabel('Vessel volume ratio (%)')
fontsize(10,"points")
print(gcf,'ADallpercent.jpg','-djpeg','-r300');
print(gcf,'ADallpercent.eps','-depsc','-r300');

%%
V=BrikLoad('AZ_Vasculature_WT06_75um.nii_regraw+orig.BRIK');
mask=BrikLoad ('ambmcmask+orig'); 
raw=V(:,1:end-165,:);
mask1=mask(:,1:end-165,:);
[N,edges] = histcounts(raw(mask1==1));
N1=cumsum(N);
brainsize=sum(mask1(:),'all');
[~,i]=find(N1>0.9*brainsize,1);
V3=raw;
thr=edges(i)
x=[849 927 975 1035 1167];
figure;histogram(raw(mask1==1),'FaceColor', 'k');set(gca,'YScale','log');hold on;
xline(x, '--r', linewidth=2);hold off