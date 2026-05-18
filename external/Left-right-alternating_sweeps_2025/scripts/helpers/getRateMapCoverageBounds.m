function [Pbounds, mask] = getRateMapCoverageBounds(x, y, mask, nclose, upsample)
if nargin < 4 || isempty(nclose), nclose = 0; end
if nargin < 5 || isempty(upsample), upsample = 1; end


% mask = rmap.validBin;
% x = rmap.x;
% y = rmap.y;

if upsample>1
    mask = imresize(mask, upsample);
    x = linspace(x(1), x(end), numel(x)*upsample);
    y = linspace(y(1), y(end), numel(y)*upsample);
end

x = x(:);
y = y(:);

if nclose
    se = strel('disk', nclose*upsample, 0);
    mask = imclose(mask, se);
end

b = bwboundaries(mask);
Pbounds = [];

for i = 1:numel(b)
    p = b{i};
    px = x(p(:, 1));
    py = y(p(:, 2));
    Pbounds = [Pbounds; [px, py]; nan, nan];
end

end