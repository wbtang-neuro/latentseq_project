function [pax] = plotEgoHistPolar(p)
%PLOTEGOHISTPOLAR Summary of this function goes here
%   Detailed explanation goes here
arguments
    p.egodir = [];
    p.egohist = [];
    p.gve = linspace(-pi, pi, 50);
    p.tl = [];
    p.ax = [];
    p.color = [0,1,0]*.5;
    p.modes = [nan, nan];
    p.alpha = .5;
end
p.gv = edg2cen(p.gve);

if ~isempty(p.ax)
    pax = p.ax;
elseif ~isempty(p.tl)
    ax = nexttile(p.tl);
    tnum = tilenum(ax);
    ax.delete;
    pax = polaraxes(p.tl);
    pax.Layout.Tile = tnum; %fix the tile indexing
else
    p.tl = tiledlayout('flow');
    pax = polaraxes(p.tl);
end
  
pax.ThetaZeroLocation = "top"; 
pax.ThetaTick = [0, 90, 180, 270]; 
pax.ThetaTickLabel = ["Front", "Left", "Back", "Right"]; 
pax.RTick = []; pax.FontSize = 10;

if ~isempty(p.egohist)
    polarhistogram('BinEdges', p.gve, 'BinCounts', p.egohist,...
        'FaceColor', p.color, 'FaceAlpha',p.alpha);
elseif ~isempty(p.egodir)
    polarhistogram(p.egodir,'BinEdges', p.gve, 'FaceColor', p.color, 'FaceAlpha',p.alpha);
end
polarscatter(p.modes, [1,1]*pax.RLim(2), 'k'); 
if~isnan(p.modes)
    text(pi/4, pax.RLim(2)*1.2, sprintf("%.1f", rad2deg(p.modes(1)))); 
    text(-pi/4, pax.RLim(2)*1.2, sprintf("%.1f", rad2deg(p.modes(2))));
end
%
end

