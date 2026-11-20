function seg_one(dataHead, outDir, mode, sigC)
%SEG_ONE  Segmentation of one registered volume, with the computation of
%  Fig2_VesselSegmentation/myseg_tissue_top2to50.m                (mode 'top2to50')
%  FigS6_S7_S8_AtlasComparison/upstream/myseg_tissue_less50.m     (mode 'less50')
% Input: one volume (full path to the +orig.HEAD). Output:
% outDir/<name>_resultinrange50.mat (or _resultinrangeless50.mat).
%   top2to50: Sthr 0.5, Vthr 0.98, step 0.02 (25 levels), sigmoid constant 1000
%             for both groups
%   less50  : Sthr 0.1, Vthr 0.5, step 0.1 (5 levels), sigmoid constant 500 for AD,
%             1000 for WT (as in myseg_tissue_less50.m)
% sigC must be passed explicitly so the log records it.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','_third_party','FastNonLocalMeans3D')));
cd('/path/to/vess');                      % <-- set this: WORK_DIR, where ambmcmask+orig lives (read-only)
t0 = tic;
fprintf('[%s] seg_one %s mode %s sigC %g\n', datestr(now), dataHead, mode, sigC);
mask = BrikLoad('ambmcmask+orig');
switch mode
    case 'top2to50', Vthr = 0.98; Sthr = 0.5; stp = 0.02; suffix = '_resultinrange50.mat';
    case 'less50',   Vthr = 0.5;  Sthr = 0.1; stp = 0.1;  suffix = '_resultinrangeless50.mat';
    otherwise, error('mode');
end
data = dataHead;                      % '<...>+orig.HEAD'
V = BrikLoad(data(1:end-5));

mask1=zeros(size(V));
mask1(1:339,:,1:243)=mask(1:339,:,7:249);
mask1=mask1.*mask;

Vout=maxfilt(V,3);                    % returns V unchanged (as in myseg_tissue_top2to50.m)
B = imgaussfilt3(Vout,10);
V1=exp(-(V-B)/sigC+1).*mask1;

V2=FastNonLocalMeans3D(V1,0.2);
raw=V2;
raw=raw(:,1:end-165,:);

mask1=mask1(:,1:end-165,:);
[N,edges] = histcounts(raw(mask1==1));
N1=cumsum(N);
brainsize=sum(mask1(:),'all');

Vout={};
percent=Sthr:stp:Vthr;
thrs = nan(1,numel(percent));
for perct= 1:length(percent)
[~,i]=find(N1>percent(perct)*brainsize,1);
V3=raw;
thr=edges(i);
thrs(perct)=thr;
V3(V3<thr)=0;
V3(V3~=0)=1;
Vout{perct}=V3;
end
stem = data(1:end-5);                       % '.../<name>+orig'
if ~endsWith(stem, '+orig'), error('expected +orig input'); end
nm = stem(find(stem=='/',1,'last')+1:end);  % not fileparts: names contain '.sr.nii.gz'
outFile = fullfile(outDir, [nm suffix]);   % same name rule as myseg_tissue_top2to50.m: data(1:end-5) + suffix
save(outFile,'Vout','raw',"percent",'-v7.3');
fprintf('[%s] wrote %s  nbins %d  brainsize %d  thr(1) %.6g thr(end) %.6g  (%.1f min)\n', ...
    datestr(now), outFile, numel(N), brainsize, thrs(1), thrs(end), toc(t0)/60);
end

function Vout=maxfilt(V,r)
Vout=V;
for i =1+r:size(V,1)-r
    for j=1+r:size(V,2)-r
        for k=1+r:size(V,3)-r
            V(i,j,k)=max(V(i-r:i+r,j-r:j+r,k-r:k+r),[],'all');
        end
    end
end
end
