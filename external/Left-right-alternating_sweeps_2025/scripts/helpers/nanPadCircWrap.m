function [phiPad, ijump, indsOut, indsIn] = nanPadCircWrap(phi, thresh)
% Pad circular wrapping events in data with single NaN values for plotting.

% Assumes that circular axes are in columns. If phi is a matrix, it will be
% interpreted as a torus, and padding will be added if a jump occurs on
% any axis.

if nargin < 2 || isempty(thresh), thresh = pi; end

% assert(isvector(phi), "Input 'phi' must be a vector");

[npoints, naxes] = size(phi);

ijump = findCircularJumps(phi, thresh);

% repeat indices for remaining columns
% ijumpRep = ijump + (0 : naxes-1)*npoints; % no need to repeat
njump = numel(ijump);
[phiPad, indsIn, indsOut] = insertValues(phi, ijump, nan);

nout = npoints+njump;
phiPad = reshape(phiPad, nout, naxes);

end