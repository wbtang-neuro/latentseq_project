function [xdiff] = circ_diff(x, pad)
%CIRC_DIFF Computes circular distance between adjacent elements of vector x
xdiff = circ_dist(x(2:end), x(1:end-1));
if nargin ==2
    xdiff(end+1)=nan;
end
end

