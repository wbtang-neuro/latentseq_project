function [grid, gg, gv] = getGrid(X, ng, boundPercentile)
if nargin < 3 || isempty(boundPercentile), boundPercentile = 0; end
    
if boundPercentile
    p = boundPercentile;
    bounds = prctile(X, [p, 100-p], 1)';
else
    bounds = [min(X)', max(X)'];
end

[grid, gg, gv] = generateGridLocal(bounds, ng, X);

end

function [ggAll, gg, gv] = generateGridLocal(lims, ng, X)

nc = size(lims, 1);

for ii=1:nc
    gv{ii} = linspace(lims(ii,1), lims(ii,2), ng);
    gv{ii} = cast(gv{ii}, "like", X);
end
[gg{1:nc}] = ndgrid(gv{:});

ggAll = [];
for ii=1:nc
    ggAll = [ggAll, gg{ii}(:)];
end

end