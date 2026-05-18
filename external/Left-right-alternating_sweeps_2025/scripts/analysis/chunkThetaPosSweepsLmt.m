function sweeps = chunkThetaPosSweepsLmt(D, posNeural, Y, plt)
% new version which doesn't dependend on raw LMT fit data. Assumes that all
% of the supplies units are tuned to the position, with no weighting of
% firing rates.

if nargin < 4 || isempty(plt), plt = false; end

S = SweepsSettings();

% Y = reconstructSpikeCounts({Us.spikeInds}, D.nt);

pTrue = [D.x, D.y];
p = posNeural;
p = gsmooth(p, S.tsm.pos);
pSlow = gsmooth(p, S.tsm.pos_slow);
pHpf = p-pSlow;

t = D.t;
Ysum = sum(Y, 2); % activity weighted by tuning selectivity
Yv = Y>0;

chk = D.thetaChunks;
ncyc = numel(chk.iStart);

for c = ncyc:-1:1 % avoid re-allocating by running backwards
    
    i0 = chk.iStart(c);
    if c==ncyc
        i1 = numel(t);
    else
        i1 = chk.iStart(c+1);
    end
    ic = i0 : i1;
    Xc = p(ic, :);
    
    Yc = Y(ic, :);
    Ysumc = Ysum(ic);
    Yvc = Yv(ic, :);
    
    hiAct = Ysumc > median(Ysumc);
    
    % Find most distal point in the sweep (relative to the start position)
    dp = Xc-Xc(1, :);
    d = hypot(dp(:, 1), dp(:, 2));
    dn = d;
    dn(~hiAct) = nan;
    [dmx, imx0] = max(dn);
    imx = imx0 + i0 - 1;
    isweep = i0 : imx;
    
    % "Line score"
    line = [0, 0, dp(imx0, :)]; % sweep best-fit line
    dSweep = d(2:imx0-1);
    dpSweep = dp(2:imx0-1, :);
    dists = distancePointLine(dpSweep, line); % disance of each point to line
%     distsScaled = dists ./ dmx;
%     distsScaled = dists ./ dmx;
    distsScaled = dists ./ dSweep;

    s = struct();
    s.lineScore = mean(distsScaled);
    
    s.iStart = i0;
    s.iStop = imx;
    s.tStart = t(s.iStart);
    s.tStop = t(s.iStop);
    s.posStart = p(s.iStart, :);
    s.posStop = p(s.iStop, :);
    
    s.pos = p(isweep, :);
    s.posHpf = pHpf(isweep, :);

    s.posAll = p(ic, :);
    s.posHpfAll = pHpf(ic, :);

    s.posTrue = pTrue(s.iStart, :);
    s.length = dmx;
    s.rateScores = Ysumc;
    s.rateScoreTotal = sum(Ysumc);
    s.rateScoreMean = mean(Ysumc);

    s.activeCellsTotal = sum(Yvc, 1) > 0;       % logical vector indicating cells that fired >=1 spike in the cycle
    s.nActiveCells = uint16(sum(Yvc, 2));        % number of cells firing >=1 spike in current time bin
    s.nActiveCellsTotal = uint16(sum(s.activeCellsTotal)); % number of cells firing >=1 spike in the cycle
    s.nSpikes = uint16(sum(Yc, 2));             % summed spike count across cells, for each sample
    s.nSpikesTotal = uint16(sum(s.nSpikes));    % total summed spike count
    s.thetaCycleInds = [i0, i1];
    
    % Displacement vector of sweep path
    dp = s.posStop - s.posStart;
    s.travelVector = dp;
    s.travelDirection = atan2(dp(2), dp(1));
    
    % Displacement vector from true to final sweep position
    dp = s.posStop - s.posTrue;
    s.truePosVector = dp;
    s.truePosDirection = atan2(dp(2), dp(1));
    
    % Displacement vector from true to final sweep position
    dp = pHpf(s.iStop, :);
    s.hpfPosVector = dp;
    s.hpfPosDirection = atan2(dp(2), dp(1));

    sweeps(c, 1) = s;
    
end

if plt
    
    if numel(plt) == 1
        trng = t(1) + [0, 10];
    else
        trng = plt;
    end

    vt = restrict(t, trng);
    vs = restrict([sweeps.tStart], trng);
    
    markerArgs = {"lineWidth", 1, "markerSize", 4};
    
    figure
    clear h
    for n = 1:2
        axs(n) = subplot(3, 1, n);
        hold on
        h.true_position = plot(t(vt), pTrue(vt, n), "color", S.col_pos_true, "lineWidth", 2);
        h.encoded_position = plot(t(vt), p(vt, n), "color", S.col_pos_sweep, "lineWidth", 1);
        pstart = cat(1, sweeps(vs).posStart);
        pstop = cat(1, sweeps(vs).posStop);
        h.sweep_start = plot([sweeps(vs).tStart], pstart(:, n), 'o', "color", S.col_start, markerArgs{:});
        h.sweep_end = plot([sweeps(vs).tStop], pstop(:, n), 'o', "color", S.col_stop, markerArgs{:});
        if n==2, xlabel("Time (s)"); end
        ylabel(char('x'+n-1)+"-position (m)");
        xlim(trng);
        axs(n).XAxis.Visible = "off";
        icyc = cat(1, sweeps(vs).thetaCycleInds);
        tcycstart = t(icyc(:, 1));
        plotShadedThetaCycles(axs(n), tcycstart);
    end
    legendFromStruct(h);

    axs(3) = subplot(3, 1, 3);
    clear h
    z = Y(vt, :);
    z = gsmooth(z, 1);
    ymax = size(spikeCounts, 2);
    imagesc(t(vt), 1:ymax, z');
    clim([0, prctile(z(:), 99)]);
    z = sum(z, 2);
    z = z./max(z)*ymax;
    plot(t(vt), z, "w", "lineWidth", 1);
    plotShadedThetaCycles(axs(3), tcycstart);
    
    linkaxes(axs, "x");
    axs(2).Clipping = "off";
    axs(2).YLim = axs(2).YLim;
    x = axs(2).XLim(1) + [0, 0.5];
    y = axs(2).YLim(1) - 0.1 + [0, 0];
    plot(axs(2), x, y, "k", "lineWidth", 3);
    
    
    if numel(plt) == 2
        % Plot 2-D latent trajectory and sweep start/end pos
        figure
        
        x = pTrue(vt, 1);
        y = pTrue(vt, 2);
        
        ax = subplot(1, 2, 1);
        h = sweepTrajectoryPlot(ax, x, y, sweeps(vs), "fullTrajectory");
        title("Full latent trajectory");
        
        ax = subplot(1, 2, 2);
        h = sweepTrajectoryPlot(ax, x, y, sweeps(vs), "sweepsOnly");
        title("Sweeps only")
        
        % 0.2 m scale bar
        ax.Clipping = "off";
        pAll = cat(1, sweeps(vs).posAll);
        pmin = min(pAll);
        x = pmin(1) - 0.1 + [0, 0.2];
        y = pmin(2) - 0.2 + [0, 0];
        axis(axis(ax));
        plot(ax, x, y, "k", "lineWidth", 3);
    end
    
end

end