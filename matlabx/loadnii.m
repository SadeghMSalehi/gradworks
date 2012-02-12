function [ imgOut ] = loadnii( fname, flip )
%LOADNII 
    if (nargin == 1) 
        flip = 0;
    end
    niiIn = load_nii(fname);
    imgOut = rescaleimg(im2double(niiIn.img));
    if (flip ~= 0)
        imgOut = slicefun(@(x)(rot90(x,flip)), imgOut);
    end
    niiIn = [];
end