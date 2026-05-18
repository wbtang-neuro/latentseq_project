function [outputArg1,outputArg2] = plotSwsExample(D,newdec)
%% Plot decoded trajectory with ID and velocity vectors

m = 2;
nsm = 0;
res = newdec.gridDecoding;


posSlow = res(2).peakPos;
posFast = res(1).peakPos;
pos = posFast;
data = D.data;
% Get the offset between the fast and slow decoded position,
% and wrap these offsets onto hex tile
info = res(1).info; % doesn't matter which res we use here
rot = res(1).info.acorrRotation;
timeRange = [-.2,3] + 8725;
vt = restrictq(D.t, timeRange);
data.idprob = gsmooth(D.idspk, 1);
%%
clf; 
tiledlayout(6, 3, 'TileSpacing','tight');
timeRange = [9153, 9156];
vt = restrictq(data.t, timeRange);
axs(1) = nexttile([1, 3]);
ax = gca;
bar(data.t(vt), data.idprob(vt), 'FaceColor',[.2,.2,.2])
ax.FontSize = 12;
ax.YTick = [0, .2]; ax.XTick = [];
ax.XLim = timeRange;
axs(2) = nexttile([2, 3]);
ax = gca;
units = D.units.mec;
idmvl = getNestedField(units, 'rmf.id.mvl');
idmu = getNestedField(units, 'rmf.id.centerOfMass');
thetamvl = getNestedField(units,'rmf.theta.mvl');
thetamu = getNestedField(units, 'rmf.theta.centerOfMass');
isId = idmvl>.4 & thetamvl>.3;
idx = isId;
data.units = units(isId); data.idmu=(idmu(isId));
raster2(data.units, ax=ax, yVar=data.idmu, timeRange=timeRange, tickWidth=.5);
raster2(data.units, ax=ax, yVar=data.idmu+2*pi, timeRange=timeRange, tickWidth=.5);
plot(D.t(vt), gsmooth(D.hd(vt), 2), 'b')
ax.XTickLabel = [];
ax.YLim = [-pi, pi+.5];
ax.YTick = [-pi, pi]+.25; ax.YTickLabel = ["0", "360"];
ax.FontSize = 12;
% 1 s scale bar
ax = axs(1);
x0 = ax.XLim(1);
y0 = ax.YLim(2);
% ax.Clipping = "off";
ax.XLim = ax.XLim;
ax.YLim = ax.YLim;
ax.XTick = [];
scale_bar = plot(ax, x0+[0, .5], y0 + 0 + [0, 0], "k", "lineWidth", .5);

% timeRanges = [9154.22, 9154.3; 9154.695, 9154.99; 9155.021, 9155.12];
timeRanges = [9154.695, 9154.86; 9154.86, 9154.97; 9155.021, 9155.14];

timeRanges = [9154.695, 9155.14];
for t = 1:size(timeRanges, 1)
    timeRange = timeRanges(t, :);
    vt = restrictq(data.t, timeRange);
    
    id = data.id;
    inds = find(vt);
    
    [x, y] = unwrapPosHex(pos(inds, 1), pos(inds, 2), info.idealAxes);
    x(~isnan(x)) = gsmooth(x(~isnan(x)), 1);
    y(~isnan(y)) = gsmooth(y(~isnan(y)), 1);
    inan = isnan(x);

    nexttile([3,1])
    plotGridTiles(gca, info.idealAxes);
    plot(x, y, 'k', 'LineWidth',.25);
    colormap cool
    scatter(x, y, 8, inds, 'MarkerFaceAlpha',.8, 'MarkerEdgeColor','k')
    
    
    idx = round(numel(inds)./2);
     peakt = D.t(D.thetaChunks.icen);
    [u,v]=pol2cart(data.id(inds), data.idprob(inds)*1*0+.35);

    quiver(x(idx), y(idx), u(idx), v(idx), 'AutoScale','off', 'Color','r', 'LineWidth',.5, 'MaxHeadSize',0);
    scatter(x(idx), y(idx), 8, 'r')


    axis square
    
    idx = find(~isnan(x));
    set(gca, 'XLim', [-.75,.55]+.2, 'YLim', [-.75,.55]+.33, 'XTick', [], 'YTick', []);

    timeRanges(t, :)-timeRanges(1)
end
linkaxes(axs, 'x');
xline(axs(1), [timeRanges(1), timeRanges(end)], 'r')
xline(axs(2), [timeRanges(1), timeRanges(end)], 'r')
xline(axs(2), [timeRanges(:)], 'r')
peaks=data.idprob(D.thetaChunks.icen);
idx = peakt>timeRanges(1) & peakt<timeRanges(end);% & peaks>.15;
scatter(axs(1), peakt(idx), peaks(idx), 10, 'r')
end

