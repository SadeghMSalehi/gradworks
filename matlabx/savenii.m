function [ output_args ] = savenii( img, fname )
%savenii (img, fname)
    fname = sprintf('../Data/%s', fname);
    save_nii(make_nii(slicefun(@(x)(rot90(x,-1)), img)),fname);
end

