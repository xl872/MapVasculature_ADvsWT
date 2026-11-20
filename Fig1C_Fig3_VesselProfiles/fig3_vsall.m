function fig3_vsall(workDir, folder, outFile, sn)
%FIG3_VSALL  VC, VW, numPixels of data/vsall<WT|AD>.mat (Fig 3): vessel profiles,
%widths and cross-section sizes of every animal of one group (folder 'WT' or 'AD')
%at the 25 segmentation thresholds (tag = 1:25), on coronal slice sn, with the
%profile computation of the first section of vesselsize.m. Each animal's Vout is
%loaded once and all 25 tags are processed.
%   fig3_vsall(workDir, 'WT', 'vsallWT.mat', 413)
cd(workDir);
% sn: coronal slice index (data/vsall*.mat: sn = 413)
VC=[];VW=[];
namelist = dir(['./',folder,'/*resultinrange50.mat']);
fprintf('[%s] %s: %d files in %s\n', datestr(now), folder, numel(namelist), workDir);
numPixels = cell(25, numel(namelist));
rawUsed = strings(numel(namelist),1);
for index=1:length(namelist)
    data=[namelist(index).folder,'/',namelist(index).name];
    S = load(data,'Vout');
    % Profile volume: the segmented reg1 volume itself, [name(1:end-20) '.BRIK'] (= 'su' in
    % beforeafter_vesselsize.m).
    rf = [namelist(index).folder,'/',namelist(index).name(1:end-20),'.BRIK'];
    raw = BrikLoad (rf);
    rawUsed(index) = string(namelist(index).name(1:end-20));
    namenew = struct('name', char(rawUsed(index)));
    raw=raw(:,91:end-165,:);
    rawRot=imrotate3(raw,-50,[1 1 0],"nearest","loose");
    sliceraw=squeeze(rawRot(:,131:end,sn));
    for tag=1:25
        %0:raw, 1~25:0.5~0.98 percent brain voxels
        dF=S.Vout{tag};
        dF=dF(:,91:end,:);
        volRot = imrotate3(dF,-50,[1 1 0],"nearest","loose");
        slice = squeeze(volRot(:,131:end,sn));
        CC=bwconncomp(slice);
        numPixels{tag,index} = cellfun(@numel,CC.PixelIdxList);
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
    end
    fprintf('[%s]   %2d %s  (raw for profiles: %s)  n vessels tag1 %d tag25 %d\n', datestr(now), index, ...
        namelist(index).name, namenew(1).name, numel(numPixels{1,index}), numel(numPixels{25,index}));
    clear S
end
files = string({namelist.name})';
save(outFile,'VC','VW','numPixels','files','rawUsed','sn','-v7');
end
