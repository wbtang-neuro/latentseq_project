function [outputArg1,outputArg2] = plotRemExample(D,newdec)
%PLOTREMFIXBAYES Summary of this function goes here
%   Detailed explanation goes here
%% Plot decoded trajectory with ID and velocity vectors
m = 2;
nsm = 0;
res = newdec.gridDecoding;


posSlow = res(2).peakPos;
posFast = res(1).peakPos;

data = D.data;
% Get the offset between the fast and slow decoded position,
% and wrap these offsets onto hex tile
info = res(1).info; % doesn't matter which res we use here
rot = res(1).info.acorrRotation;
%         info = dec.gridDecoding{1}(1).info
posOffsets = wrapPosToGridTile(posFast-posSlow, info.idealAxes);

% Unwrap the slow-moving path, and add the wrapped offsets
[xu, yu] = unwrapPosHex(posSlow(:, 1), posSlow(:, 2), info.idealAxes);
posFastU =  gsmooth([xu, yu], 2) + posOffsets;
%     posFastU = rotate2d(fliplr(posFastU), -rot);
%     posSlowU = rotate2d([yu, xu], -rot);
posSlowU = rotate2d([xu, yu], -rot);
posSlowU = [xu, yu];
dec.mods{m}.posSlow = posSlowU;
dec.mods{m}.posFast = posFastU;
dec.mods{m}.posOffset = posOffsets;

pos = dec.mods{m}.posFast;
possm = dec.mods{m}.posSlow;
possm = movmean(possm, 4);
possmAll = possm;
timeRange = [-.2,3] + 8725;
vt = restrictq(D.t, timeRange);
possm = possm(vt, :);
% possm = rotate2d(fliplr(possm), -rot);

% pos(res.peakProb<.03, :)=nan;
% pos(data.gprob<.02, :)=nan;

%%
clf; 
nswp = 4;
tiledlayout(6, nswp, 'TileSpacing','tight');
timeRange = [+.1,2.5] + 8725;
vt = restrictq(D.t, timeRange);
axs(1) = nexttile([1, nswp]);
ax = gca;

data.idprob = gsmooth(D.idspk, 1);
bar(D.t(vt), data.idprob(vt), 'FaceColor',[.2,.2,.2])
ax.FontSize = 12;

ax.YTick = [0, .2]; ax.XTick = [];
ax.XLim = timeRange;
axs(2) = nexttile([2, nswp]);
ax = gca;

units = D.units.mec;
idmvl = getNestedField(units, 'rmf.id.mvl');
idmu = getNestedField(units, 'rmf.id.centerOfMass');
thetamvl = getNestedField(units,'rmf.theta.mvl');
thetamu = getNestedField(units, 'rmf.theta.centerOfMass');
isId = idmvl>.4 & thetamvl>.3;
idx = isId;
data.units = units(isId); data.idmu=(idmu(isId));
% raster2(data.units, ax=ax, yVar=data.idmu, timeRange=timeRange, tickWidth=.5, yHeight=7);
raster2(data.units, ax=ax, yVar=data.idmu, timeRange=timeRange, tickWidth=.5);
raster2(data.units, ax=ax, yVar=data.idmu+2*pi, timeRange=timeRange, tickWidth=.5);
plot(D.t(vt), gsmooth(D.hd(vt), 2), 'b')
% plot()
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
% ax.Clipping = "off";
% timeRanges = [9154.22, 9154.3; 9154.695, 9154.99; 9155.021, 9155.12];
% timeRanges = [9154.695, 9154.86; 9154.86, 9154.97; 9155.021, 9155.14];
% 
% m = 3;
sweeps = newdec.sweeps;
vswp = restrictq([sweeps.tStart]', timeRange);
sweeps = sweeps(vswp);
for s = 1:numel(sweeps)
    timeRanges(s, :) = [sweeps(s).tStart, sweeps(s).tStop];
end
i = 4+2+4;
i = 9;
timeRanges = timeRanges(i:(i+3), :);
nsweeps = size(timeRanges, 1);

swps = sweeps(i:(i+3))
% pos0 = possm(sweeps(1).iStart, :);
pos0 = possm(1, :)+[.25,-.45];
pos = pos-pos0;
possm = possm-pos0;
possmAll = possmAll-pos0;
% possm = rotate2d(fliplr(possm), -rot);
for t = 1:nsweeps
    timeRange = timeRanges(t, :);
    vt = restrictq(D.t, timeRange);
%     vt = vts{t};
    id = data.id;
    inds = find(vt);
%     inds = swps(t).iStart:swps(t).iStop2;
%     inds = swps(t).iStart:swps(t).iStop;
%     inds = swps(t).iSweep;
    vt(:)=false; vt(inds)=true;
    
    [x, y] = unwrapPosHex(pos(inds, 1), pos(inds, 2), info.idealAxes);
    x = pos(inds, 1); y = pos(inds, 2)
    x(~isnan(x)) = gsmooth(x(~isnan(x)), 1);
    y(~isnan(y)) = gsmooth(y(~isnan(y)), 1);
    inan = isnan(x);
    x(inan) = interp1(1:sum(~inan), x(~inan), find(inan), 'linear', 'extrap');
    y(inan) = interp1(1:sum(~inan), y(~inan), find(inan), 'linear', 'extrap');
    

    nexttile([3,1])
    plotGridTiles(gca, info.idealAxes);
    x = gsmooth(x, nsm); y = gsmooth(y, nsm);
%     x = swps(t).posHpf(:, 1); y= swps(t).posHpf(:, 2);
    plot(possm(:, 1), possm(:, 2), 'LineWidth',.5, 'Color', [1,1,1]*.7)
    plot(x, y, 'LineWidth',.25,'Color', 'k');
    colormap cool
    scatter(x, y, 8, inds, 'MarkerFaceAlpha',.8, 'MarkerEdgeColor','k', 'LineWidth',.25)
    
    
    idx = round(numel(inds)./2);
%     peakt = data.t(data.idpeaks.loc); peaks=data.idpeaks.peak;
    peakt = D.t(D.thetaChunks.icen);
    idx = peakt>timeRange(1) & peakt<timeRange(end);% & peaks>.15;
    idx = interp1(D.t(vt), 1:numel(inds), peakt(idx), 'nearest');
%     idx = (1);
    [u,v]=pol2cart(data.id(inds), data.idprob(inds)*1*0+.25);


    x = possmAll(inds, 1); y = possmAll(inds, 2);
%     x = possmAll(:, 1); y = possmAll(:, 2);
%     idx = sweeps(t+11).iStart
    quiver(x(idx), y(idx), u(idx), v(idx), 'AutoScale','off', 'Color','r', 'LineWidth',.5, 'MaxHeadSize',0);
    scatter(x(idx), y(idx), 8, 'r')
%     quiver(x(1), y(1), u(1), v(1), 'AutoScale','off', 'Color','r', 'LineWidth',1.5, 'MaxHeadSize',.9);


    axis square
    
    idx = find(~isnan(x));
    set(gca, 'XLim', [-.75,.4]+.2, 'YLim', [-.75,.35]+.3, 'XTick', [], 'YTick', []);
%   axis off
%     for i = 1 % Iterate through samples and keep the one thats closest to the previous
% 
%         
%         posTiled = [x(:), y(:)];
%         posTiled = rotate2d(fliplr(posTiled), -rot);
%         x = posTiled(:, 1); y = posTiled(:, 2);
%         nexttile
%         scatter(x, y, 15, repmat(inds, 7, 1))
%         scatter(x(1), y(1), 30, 'r')
%     
%         [u,v]=pol2cart(repmat(data.id(inds), 7, 1), repmat(data.idprob(inds), 7, 1)*.3);
%         quiver(x, y, u, v);
%     
%     % set(gca, 'XLim', [-.5,.5], 'YLim', [-.5,.5]);
%         axis square
%         title(i)
%     
%     end
end
linkaxes(axs, 'x');
xline(axs(1), [timeRanges(1), timeRanges(end)], 'r')
xline(axs(2), [timeRanges(1), timeRanges(end)], 'r')
% xline(axs(2), [timeRanges(:)], 'r')
peaks=data.idprob(D.thetaChunks.icen);
idx = peakt>timeRanges(1) & peakt<timeRanges(end);% & peaks>.15;
scatter(axs(1), peakt(idx), peaks(idx), 10, 'r')
end

