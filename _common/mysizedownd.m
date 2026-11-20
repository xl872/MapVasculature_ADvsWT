function output = mysizedownd(input,sr,mask0)
    
    input1=zeros(ceil(size(input,1)/sr)*sr,ceil(size(input,2)/sr)*sr,ceil(size(input,3)/sr)*sr);
    mask=input1;
    input1(1:size(input,1),1:size(input,2),1:size(input,3))=input;
    % mask0=input;mask0(~isnan(input))=1;mask0(isnan(input))=0;
    mask(1:size(input,1),1:size(input,2),1:size(input,3))=mask0;
    
    Vdv=squeeze(sum(reshape(input1,sr,size(input1,1)/sr,size(input1,2),size(input1,3)),1,'omitnan'));
    Vdv=squeeze(sum(reshape(Vdv,size(input1,1)/sr,sr,size(input1,2)/sr,size(input1,3)),2,'omitnan'));
    Vdv=squeeze(sum(reshape(Vdv,size(input1,1)/sr,size(input1,2)/sr,sr,size(input1,3)/sr),3,'omitnan'));
    mask=squeeze(sum(reshape(mask,sr,size(input1,1)/sr,size(input1,2),size(input1,3)),1,'omitnan'));
    mask=squeeze(sum(reshape(mask,size(input1,1)/sr,sr,size(input1,2)/sr,size(input1,3)),2,'omitnan'));
    mask=squeeze(sum(reshape(mask,size(input1,1)/sr,size(input1,2)/sr,sr,size(input1,3)/sr),3,'omitnan'));

    output=Vdv./mask;
    output(isnan(output))=0;
    output(mask<sr^3*0.1)=0;
end