function sweeps = chunkThetaPosSweeps(D, dec, Y, plt)
% new version which doesn't dependend on raw LMT fit data. Assumes that all
% of the supplies units are tuned to the position, with no weighting of
% firing rates.

if nargin < 4 || isempty(plt), plt = false; end

S = sweepsSettings();

% Y = reconstructSpikeCounts({Us.spikeInds}, D.nt);

pTrue = [D.x, D.y];
p = dec.decpos;
p = gsmooth(p, .8);
% p = gsmooth(p, 1);
% pSlow = gsmooth(dec.possm, S.tsm.pos);
pSlow = gsmooth(dec.possm,1);
pHpf = p-pSlow;
% prob = dec.maxprob;

t = D.t;
Ysum = sum(Y, 2); % activity weighted by tuning selectivity
Yv = Y>0;
Ysm = gsmooth(Ysum, 1);

chk = D.thetaChunks;
ncyc = numel(chk.iStart);

maxAngleDiff = pi/2;
maxDist = .2;
% Compute travel direction between points
pdiff = diff(p);
travelDir = atan2(pdiff(:, 2), pdiff(:, 1));
travelDirdiff = circ_diff(travelDir);
travelDirdiff = [0;abs(travelDirdiff);0];
% nextdirdiff = [abs(travelDirdiff);0;0];
baddirdiff = travelDirdiff>maxAngleDiff;

% Copmute distance between points
nextdist = hypot(pdiff(:, 1), pdiff(:, 2));
nextdist(end+1) = nan;
prevdist = circshift(nextdist, 1);
% maxDist = prctile(nextdist(1:end-1), 90);
% Scrap bins with too few active cells
if isfield(dec, "hc")
    p = dec.decpos;
    p = gsmooth(p, 1);
    p(dec.nactive<1, :) = nan;
    p(find(dec.err>.5), :) = nan;
    p(dec.maxprob(:)<0, :) = nan;
else 
    p(dec.nactive<5, :) = nan;
    p(find(dec.err>.15), :) = nan;
    p(dec.maxprob(:)<dec.shuff99, :) = nan;
end
disp(size(p))
pHpf = p-pSlow;
%%
for c = ncyc:-1:1 % avoid re-allocating by running backwards
    %

    i0 = chk.iStart(c);
    if c==ncyc
        i1 = numel(t);
    else
        i1 = chk.iStart(c+1);
    end
    ic = i0 : i1;
    pos = p(ic, :);
    
    Yc = Y(ic, :);
    Ysumc = Ysum(ic);
    Yvc = Yv(ic, :);
    Ysmc = Ysm(ic);

%     hiAct = Ysumc > median(Ysumc); % Throws away half of the samples? It says that the distal point need to be over the 50th perc
     
    % Find the point with highest prob
%     probc = prob(ic);
%     probc(~hiAct)=nan;
%     [~,imx] = max(probc);
    [~,imx] = max(Ysmc);
%     posmx = pos(imx, :);
%     dirdiffc = travelDirdiff(ic(1:end-1));
    
    % Find adjacent samples less than 10cm away
    badnxt = find(nextdist(ic)>maxDist|baddirdiff(ic)|isnan(pos(:, 1)));
    badprv = find(prevdist(ic)>maxDist|baddirdiff(ic)|isnan(pos(:, 1)));

%     badpos = find(dpos>.2 | dirdiffc>maxAngleDiff);
    istart = 1;
    iend = numel(ic);
    if ~isempty(badprv)
        istart = max([badprv(badprv<imx); 1]);
    end
    if ~isempty(badnxt)
        iend = min([badnxt(badnxt>=imx); iend]);
    end
    valid = istart:iend;
    invalid = true(size(ic));
    invalid(valid)= false;
    pos(invalid, :)=nan;
%     pos(valid, :) = gsmooth(pos(valid, :), .5);
    % Define end point as the point that maximize distance from point 1
    dpos = triu(squareform(pdist(pos)));
    if isempty(dpos)
        dpos = [nan];
    end
%     dpos(:, ~hiAct)=nan;
    [dmx, imxdist]=max(dpos(:));
    [istart,imx0] = ind2sub(size(dpos), imxdist);
    nvalidold = numel(istart:imx0);
    %% Try an alternative approach that maximizes distance from theta0
    hpfpos = pHpf(ic, :);
    hpfpos(invalid, :) = nan;
    hpfdist = hypot(hpfpos(:, 1), hpfpos(:, 2));
    [~, imx0] = max(hpfdist);
    [~, istart] = min(hpfdist);
    if istart>imx0
        istart = 1;
    end
%     dn = d;
%     dn(~hiAct) = nan;
%     [dmx, imx0] = max(dn); %update iend
    imx = imx0 + i0 - 1;
    iswp0 = i0 + istart - 1;
    isweep = iswp0 : imx;
    
    invalid = true(size(ic));
    invalid(istart:imx0)= false;
    pos(invalid, :)=nan;
    hpfpos(invalid, :)=nan;
    
    % define sweep axis
    
%     swpline = [0, 0, hpfpos(imx0, :)]; % sweep best-fit line
    ortend = rotate2d(hpfpos(imx0, :), pi/2); % Orthogonal axis
    ortline = [0,0,ortend];
%     dists = distancePointLine(hpfpos, swpline); % disance of each point to line
    hpfpos(imx0,:)=nan;
    distsort = distancePointLine(hpfpos, ortline);
%     swpvar = var(dists, 'omitnan');
    totalvar = sum(var(hpfpos, 'omitnan'), 'omitnan');
    ortvar = var(distsort, 'omitnan');
    varratio = ortvar./totalvar;

    dpos = pos-pos(istart, :);
    d = hypot(dpos(:, 1), dpos(:, 2));
    
    %

    % "Line score"
    distsScaled = nan;
    dists = nan;
    if numel(isweep)>2
        line = [0, 0, dpos(imx0, :)]; % sweep best-fit line
        dSweep = d(2:imx0-1);
        dpSweep = dpos(2:imx0-1, :);
        dists = distancePointLine(dpSweep, line); % disance of each point to line
    %     distsScaled = dists ./ dmx;
    %     distsScaled = dists ./ dmx;
%         distsScaled = dists ./ dSweep;
    end
    % Turuosity
    diffs = [0; vecnorm(diff(pos), 2,2)];
    pathlength = sum(diffs(istart:imx0), 'omitnan');

    s = struct();
    s.lineScore = mean(dists);
    s.pathlength = pathlength;
    s.netlength = dmx;
    s.turtuosity = dmx./pathlength;
%     s.pos = pos;%tst
%     s.intrainds = istart:imx0;%tst
    s.iStart = iswp0;
    s.iStop = imx;
    s.iSweep = isweep;
    s.iStop2 = iswp0+iend-1;
    s.tStart = t(s.iStart);
    s.tStop = t(s.iStop);
    s.posStart = p(s.iStart, :);
    s.posStop = p(s.iStop, :);
    s.nvalid = numel(iswp0:imx);
    
    s.pos = p(isweep, :);
    s.posHpf = pHpf(isweep, :);

    s.posAll = p(ic, :);
    s.posHpfAll = pHpf(ic, :);

    s.posTrue = pTrue(s.iStart, :);
    s.length = dmx;

    s.activeCellsTotal = sum(Yvc, 1) > 0;       % logical vector indicating cells that fired >=1 spike in the cycle
    s.nActiveCells = uint16(sum(Yvc, 2));        % number of cells firing >=1 spike in current time bin
    s.nActiveCellsTotal = uint16(sum(s.activeCellsTotal)); % number of cells firing >=1 spike in the cycle
    s.nSpikes = uint16(sum(Yc, 2));             % summed spike count across cells, for each sample
    s.nSpikesTotal = uint16(sum(s.nSpikes));    % total summed spike count
    s.thetaCycleInds = [i0, i1];
    
    % Displacement vector of sweep path
    dpos = s.posStop - s.posStart;
    s.travelVector = dpos;
    s.travelDirection = atan2(dpos(2), dpos(1));
    
    % Displacement vector from true to final sweep position
    dpos = s.posStop - s.posTrue;
    s.truePosVector = dpos;
    s.truePosDirection = atan2(dpos(2), dpos(1));
    
    % Displacement vector from true to final sweep position
    dpos = pHpf(s.iStop, :); 
%     inds = max(s.iStop-1, s.iStart):s.iStop;
%     dpos = mean(pHpf(inds, :), 1);
%     dpos = pSlow(s.iStart, :) - s.posTrue; % Should prob change to this
    s.hpfPosVector = dpos;
    s.hpfLength = hypot(dpos(1), dpos(2));
    s.hpfPosDirection = atan2(dpos(2), dpos(1));

    s.straight = varratio;
    s.nvalidold = nvalidold;

    sweeps(c, 1) = s;
    % Some measures
    % Turtuosity (path length/net length)
    % Straightness index
    
    % Line fit
    % var expl by pc2 over pc1
    % Corr between 
    
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
    

    
end

end
