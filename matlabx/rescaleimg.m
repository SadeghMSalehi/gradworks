function [ imgOut ] = rescaleimg( img )
%RESCALEIMG Summary of this function goes here
%   Detailed explanation goes here
    minI = min(img(:));
    maxI = max(img(:));
    imgOut = (img-minI)/(maxI-minI);
end

