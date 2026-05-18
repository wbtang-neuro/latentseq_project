function h = plotExampleIdCells(units, p)
arguments
    units
    p.plotPos = false;
end
%
S = SweepsSettings;
nu = numel(units);
tl = tiledlayout(2, nu, TileSpacing="compact", TileIndexing="columnmajor");
if p.plotPos
    tl = tiledlayout(3, nu, TileIndexing="columnmajor");
end

for u = 1:nu
    % Plot hd
    h = plotRateMap(units(u).hdTuning, "ax", nexttile, "colors", S.col_id);
    title(sprintf("%su%s", units(u).recName, units(u).id));
    pax = h.axes;
    text(45, pax.RLim(2), sprintf("%.0fHz", units(u).hdTuning.maxRate));
    if u==1
        text(pax, pi, 1.5, 'HD tuning', 'Rotation', 90, 'HorizontalAlignment', 'center')
    end
    % Plot spatial
    if p.plotPos
        h = plotRateMap(units(u).posTuning, "ax", nexttile);
        text(h.axes.XLim(2), h.axes.YLim(2), sprintf("%.0fHz", units(u).posTuning.maxRate))
        if u==1
            text(h.axes.XLim(1)*1.5, h.axes.YLim(2)*.2,"Pos tuning",'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
        end
    end
    
    % Plot acorr
    ax = nexttile;
    bar(ax, units(u).tempAcorr.lags, units(u).tempAcorr.rate, 1, 'FaceColor', 'k');
    axis square;
    xlabel("Lag (s)")
    if u==1
        ylabel("N spikes")
    end
    ax.XLim = [-.5, .5];
    ax.XTick = [-.5,0,.5];
    ax.YTick = [0, ax.YLim(2)];
    if u>1
        ylabel([])
    end
end

end
