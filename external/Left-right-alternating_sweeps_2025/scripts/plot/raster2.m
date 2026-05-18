function handl = raster2(units, kwargs)
%RASTER Summary of this function goes here
%   Detailed explanation goes here

arguments
    units
    kwargs.yVar = [units.shankPos];
    kwargs.timeRange = [0,100]+min(units(1).spikeTimes);
    kwargs.color = [0,0,0]
    kwargs.yOffset = 0;
    kwargs.yOffsetSpk = [];
    kwargs.vunits = [];
    kwargs.ax = gca;
    kwargs.yHeight = [];
    kwargs.tickWidth = 1;
    kwargs.nspikes = [];
end
p = kwargs;

yrange = [min(p.yVar), max(p.yVar)];
if isempty(p.yHeight)
   p.yHeight = .03*diff(yrange);
end

nUnits = numel(units);
xAll = nan;
yAll = nan;
for u = 1:nUnits
    unit = units(u);
    [~,t] = restrictq(unit.spikeTimes, p.timeRange);
    x = t';
    
    dv = p.yVar(u) - p.yOffset;
    y = repmat(dv, 1, numel(t));
    if~isempty(p.yOffsetSpk)
        offset = interp1(p.yOffsetSpk(:, 1), p.yOffsetSpk(:, 2), t, 'nearest');
        y = wrapToPi(y-offset');
    end

    xAll = [xAll,x];
    yAll = [yAll,y];
end

if ~isempty(p.nspikes)
    nspk = numel(xAll);
    if nspk>p.nspikes
        idx = randsample(1:nspk,p.nspikes);
        xAll = xAll(idx); yAll = yAll(idx);
    end
end
xAll = repmat(xAll, 3, 1);
xAll(3, :) = nan;
xAll = xAll(:);

yAll = repmat(yAll, 3, 1) + [.5;-.5; 0]*p.yHeight;
yAll(3, :) = nan;
yAll = yAll(:);

color = [p.color];

handl = plot(p.ax, xAll, yAll, 'Color',color, 'LineWidth', p.tickWidth);

end


