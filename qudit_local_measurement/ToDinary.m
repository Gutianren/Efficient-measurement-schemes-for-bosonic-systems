function [ vector ] = ToDinary(d, j, Nq )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    
    for k = 1 : Nq
        vector(Nq - k + 1) = mod(j,d);
        j = floor(j/d);
    end

end
