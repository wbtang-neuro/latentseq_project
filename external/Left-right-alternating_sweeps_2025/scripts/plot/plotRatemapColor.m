function hmap = plotRatemapColor(rmap, color, plotcontour)
%PLOTRATEMAPCOLOR Summary of this function goes here
%   Detailed explanation goes here

c(1,1,:) = color;
c = repmat(c, [size(rmap), 1]);

hmap = imshow(c);

rmap = rmap - prctile(rmap(:), 10); rmap(rmap<0)=0;
hmap.AlphaData = rmap./prctile(rmap(:), 99);
% rmap = rmap - .05; rmap(rmap<0)=0;
% hmap.AlphaData = rmap./.5;
k = prctile(rmap(:), 90);
if plotcontour
    [m, c] = contour(rmap, 'LineColor', color);
    c.LevelList = k;
end

end

