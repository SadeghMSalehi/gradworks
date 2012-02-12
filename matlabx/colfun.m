function [ outM ] = colfun( func, M )
%COLFUN Summary of this function goes here
%   Detailed explanation goes here
    [r c] = size(M);
    outM = zeros([r c]);
    for j = 1:c
        outM(:,j) = func(M(:,j));
    end
end

