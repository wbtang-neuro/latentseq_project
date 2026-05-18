function [distVec, dist] = hexDistance(posRef,pos,gridAxes)
%HEXDISTANCE Compute distance vector between pairs of points on a grid tile
%(the vector from pos to posRef)
%
% Calulate offset between ref and target pos
posDiff = posRef-pos;

% Tile offset 
[xDiff, yDiff] = tilePointsHex(posDiff(:, 1), posDiff(:, 2), gridAxes);

% Compute the distance from zero for each tiling
dists = sqrt(xDiff.^2 + yDiff.^2); 

% Find the smallest one
[minDist, idx] = min(dists, [], 2);
[inds] = sub2ind(size(xDiff), (1:numel(minDist))', idx);

% Find the distance vector corresponding to the smallest distance
distVec = [xDiff(inds), yDiff(inds)];
dist = minDist;
end

