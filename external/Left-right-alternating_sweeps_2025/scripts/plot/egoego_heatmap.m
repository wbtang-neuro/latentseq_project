function [res, handles] = egoego_heatmap(p)
%EGOEGO_PLOT Bivariate histogram of two sets of angular data (usually ego)

arguments
    p.ax = gca;
    p.x = [];
    p.y = [];
    p.gve = linspace(-pi, pi, 101); % Bin edges
    p.gv = [];
    p.lims = 100; % Plotting lims in degrees
    p.computeStats = false; % Fraction aligned, correlation
    p.plotMarginals = false; %
    p.labels = ["Ego angle", "Ego angle"];
    p.nsmooth = 0;
    p.plot = true;

end
if isempty(p.gv)
    p.gv = rad2deg(edg2cen(p.gve));
end
    
if p.plotMarginals
    tl = tiledlayout(3,3, 'TileSpacing','tight');
    nexttile(tl, 1, [1,2])
    h = histcounts(p.x, p.gve);
    bar(p.gv, h, 'FaceColor',[.5,.5,.5])
    ax = gca; ax.XAxis.Visible = "off"; axis off; ax.XLim = [-p.lims, p.lims]; ax.YLim(2)=2*max(h);
    nexttile(tl, 3, [1,1]), axis off
    p.ax = nexttile(tl, [2,2]);
end

res.hist = histcounts2(p.y, p.x, p.gve, p.gve);%
if p.nsmooth>0
    res.hist = imgaussfilt(res.hist, p.nsmooth);
end

if p.computeStats
    valid = ~isnan(p.x) & ~isnan(p.y);
    p.x(~valid)=nan; p.y(~valid)=nan; 
    isright = p.x>0;
    isleft = p.x<0;
    
    plgivenr = sum(p.y(isright)<0)./sum(isright);
    prgivenl = sum(p.y(isleft)>0)./sum(isleft);
    res.popposite = mean([plgivenr, prgivenl]);

    plgivenr = sum(p.y(isright)>0)./sum(isright);
    prgivenl = sum(p.y(isleft)<0)./sum(isleft);
    res.psame = mean([plgivenr, prgivenl]);
    
    valid = ~isnan(p.x) & ~isnan(p.y);
    [res.corr, res.pval] = circ_corrcc(p.x(valid), p.y(valid)); 
end

if p.plot
    handles.img = imagesc(p.ax, res.hist, 'XData', p.gv, 'YData', p.gv);
    axis image; 
    ax = gca;
    ax.XLim = [-1, 1]*p.lims; ax.YLim = [-1, 1]*p.lims;
    ax.XTick = [-90, 90]; ax.YTick = [-90, 90];
    ax.XTick = [-1, 1]*p.lims; ax.YTick = [-1, 1]*p.lims;
    ax.CLim = [0, max(res.hist(:))];
    xlabel(p.labels(1));
    ylabel(p.labels(2));
    xline(0, 'w')
    yline(0, 'w')
    colormap(gca, "inferno");
    c = colorbar(ax);
    c.Label.String = "N sweeps";
    c.Ticks = [];
    handles.cbar = c;
%     c.Ticks = c.Limits;

    
%     text([-80, 80], [-80, 80], ["same", "same"],"Color",'w', 'HorizontalAlignment','center', 'FontSize',12)
%     text([-80, 80], [80, -80], ["opp.", "opp."],  "Color",'w', 'HorizontalAlignment','center', 'FontSize',12)
    
    if p.computeStats
       title(sprintf("Prct aligned: %.0f, corr: %.1f, p: %.6f", res.psame*100, res.corr, res.pval))
    end
end

end

