clc;clear all; close all;
% WORK_DIR = folder holding AD/, WT/ (registered volumes + segmentation .mat),
% the AMBMC/Allen BRIKs (HarvardDataverse/AMBMC) and the intermediate results.
WORK_DIR = '/path/to/vess';   % <-- set this
cd(WORK_DIR);                 % BrikLoad needs the HEAD/BRIK pair in the cwd
%%
sr=15;
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
climF = [-0.1, 0.1];
xn=6;yall=18;ly=1;
yn=yall/xn;
step=10;
epochs=100;
for tag=[1]
%% Ah=permute(imresize3(Mh,sr),[3,2,1]);
load(['Density',num2str(tag),'final.mat'])
load(['Density',num2str(tag),'.mat'])
PTHR=0.01;
h = figure('Position',[100,100,1800,700]);

Ah=permute(Mh,[3,2,1]);
% Significance mask: Benjamini-Hochberg FDR of ph inside the brain (mask > 0.5),
% then Gaussian smoothing of q (sigma = 1 voxel), then q < PTHR (0.01).
% The |difference| threshold follows.
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
% Ah = imgaussfilt3(Ah,1);



% dF=imresize3(squeeze(mean(mapnotresize,1)),sr);
dF=squeeze(mean(mapt,1));
% dF=squeeze(mapt(1,:,:,:));
dF=permute(dF,[3,2,1]);
Ah(abs(dF)<0.005)=0;
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
print(gcf,['new',num2str(tag),'WT-AD',num2str(PTHR),'segDendity-brainmask_cfdr',num2str(epochs),'.jpg'],'-djpeg','-r300')
% %
% h = figure('Position',[100,100,1050,400]);
% 
% % dF=MVDVn;%MVDVn MVDVWT MVDVAD
% % dF=squeeze(mean(mapt,1));
% % dF=squeeze(mapt(1,:,:,:));
% % dF=permute(dF,[3,2,1]);
% i=1;
%  tiledlayout(yn*ly+bn*lb,xn,'TileSpacing','none','Padding','tight');
% 
%     for j=1:yall
%         j=j+shiftslide;
%         nexttile([ly 1]);B=squeeze(U(:,j*step,:));
%         B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);
%         hB = imagesc(B);axis xy
%         axis image off;
%         hold on;
%         F=squeeze(dF(:,j*step,:));
%         hF = imagesc(F,climF);axis xy
%         colormap('jet')
%         alphadata =0.6*squeeze(mask(:,j*step,:));%0.6*squeeze(Ah(:,j*5,:));% 
%         set(hF,'AlphaData',alphadata);
%     end
%     saveas(gcf,['new',num2str(tag),'WT-ADsegDensity-brainmask',num2str(epochs),'.jpg'],'jpg')
%%

end
function out=myselectVout(mask,Vout,raw,tag,sr)
    out=raw;
    if tag>0&&tag<=length(Vout)
        Tout=Vout{tag};
        out(Tout==0|mask==0)=nan;
        out(~isnan(out))=1;
%     elseif tag>length(Vout)
%         Tout=Vout{length(Vout)};
%         out(Tout==1|mask==0)=nan;
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