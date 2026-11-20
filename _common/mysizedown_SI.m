function [output,validFraction] = mysizedown_SI(input,sr,minValidFraction)
%MYSIZEDOWN_SI Average finite signal values within sr-by-sr-by-sr blocks.
%   Blocks with insufficient retained voxels are returned as NaN rather than
%   zero, so missing coverage is not interpreted as a measured zero signal.

    if nargin < 3 || isempty(minValidFraction)
        minValidFraction = 0.1;
    end
    validateattributes(minValidFraction,{'numeric'}, ...
        {'scalar','real','finite','>=',0,'<=',1},mfilename,'minValidFraction');
    
    input1=zeros(ceil(size(input,1)/sr)*sr,ceil(size(input,2)/sr)*sr,ceil(size(input,3)/sr)*sr);
    mask=input1;
    input1(1:size(input,1),1:size(input,2),1:size(input,3))=input;
    mask0=input;mask0(~isnan(input))=1;mask0(isnan(input))=0;
    mask(1:size(input,1),1:size(input,2),1:size(input,3))=mask0;
    
    Vdv=squeeze(sum(reshape(input1,sr,size(input1,1)/sr,size(input1,2),size(input1,3)),1,'omitnan'));
    Vdv=squeeze(sum(reshape(Vdv,size(input1,1)/sr,sr,size(input1,2)/sr,size(input1,3)),2,'omitnan'));
    Vdv=squeeze(sum(reshape(Vdv,size(input1,1)/sr,size(input1,2)/sr,sr,size(input1,3)/sr),3,'omitnan'));
    mask=squeeze(sum(reshape(mask,sr,size(input1,1)/sr,size(input1,2),size(input1,3)),1,'omitnan'));
    mask=squeeze(sum(reshape(mask,size(input1,1)/sr,sr,size(input1,2)/sr,size(input1,3)),2,'omitnan'));
    mask=squeeze(sum(reshape(mask,size(input1,1)/sr,size(input1,2)/sr,sr,size(input1,3)/sr),3,'omitnan'));

    validFraction = mask ./ sr^3;
    output = Vdv ./ mask;
    output(validFraction < minValidFraction) = NaN;
end
