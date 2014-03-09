function [ out ] = slicefun( func, img )
%SLICEFUN execute func iteratively on each slice and gather it
%   func: function handle ex) @(x)(rot90(x,1)
%   img: 3 dimensional image

if size(img,3) < 1
    out = [];
    return;
end

% there is a case the result of func is not compatible with the original
% img
first = func(img(:,:,1));
out = zeros([size(first) size(img,3)]);
for it = 2:size(img,3)
    out(:,:,it) = func(img(:,:,it));
end

