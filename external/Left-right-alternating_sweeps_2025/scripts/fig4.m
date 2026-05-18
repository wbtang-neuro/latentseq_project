%% PV examples
sweepsSetup
% S.dataRoot = ;% Specify data root
recNames = ["26648_1", "25843_1"];
stypes = ["lt", "ww"];
timeRanges = [7078.5, 7081.5;17152.2, 17154.2];

%% Get data
datafld = fullfile(S.dataRoot, "navigation");
for r = 1:2
    fname = fullfile(datafld, stypes(r), recNames(r)+".mat");
    tmp = load(fname);
    D = tmp.Dsession;
    D.occupancy = histcounts(D.hd, S.gve.angular);
    units = D.units.mec;
    pb = ProgressBar;
    for u = 1:numel(units)
        units(u).sc = reconstructSpikeCounts(units(u).spikeInds, D.nt); 
        tc = units(u).rmf.hd;
        D.tuning(:, u) = vec(gsmooth(tc.z, 1));
        pb.update(u/numel(units))
    end
    D.sc = [units.sc];
    recs(r) = D;
end

%% Decode
for r = 1:2
    D = recs(r);
    [prob] = decodePv(tuning=D.tuning, spikeCounts=D.sc, smoothSpikes = 1);
    gv = S.gv.angular;
    idmuAll = circ_mean(gv, prob')';
    chk = D.thetaChunks;
    chk.icen = round(movmean(chk.iStart, [0,1]));
    chk.id_corr = idmuAll(chk.icen);
    chk.hd = D.hd(chk.icen);
    chk.speed = D.speed(chk.icen);
    D.thetaChunks = chk;
    recs(r)=D;
    disp(recNames(r))
end

%% Plot results
tl = tiledlayout(3,3, "TileSpacing","tight", "TileIndexing","rowmajor");
for r = 1:2
    D = recs(r);
    chk = D.thetaChunks;
    chk.id = chk.id_corr;

    vspd = D.speed>.2;
    pos = [D.x, D.y];
    if r==2
        vspeed = D.speed>.15;
        pos(~vspeed, :) = nan;
        pos = gsmooth(pos, 10);
    end
    timeRange = timeRanges(r, :);
    idx = 1:min(D.nt, 1e5);
    plot(nexttile, pos(idx, 1), pos(idx, 2), 'lineWidth',.1,'Color', [.85,.85,.85]);
    sweepIdTrajectoryPlotShaded(gca, D.t, gsmooth(D.x, 7), gsmooth(D.y, 7), timeRange, sweeps =[], id = chk, plotPath=true, plotScaleBar=true);

    nexttile
    xticks([0,360]);
    xlabel("Decoded direction (deg)")
    ax = gca;
    ax.FontSize = 12;
    
    egoid = circ_dist(chk.id, chk.hd);
    egoid = egoid-circ_mean(egoid);
    h  = histcounts(chk.id_corr(chk.speed>.2), S.gve.angular);
    bar(linspace(0, 360,60),  h./max(h), 'faceColor', [.2,.8,.2], 'FaceAlpha', .5)
    plot(linspace(0, 360,60), D.occupancy./max(D.occupancy), 'Color', [.5,.5,.5], 'LineWidth',1)
    ymax = ax.YLim(2); ylim([0, ymax]);yticks([0, (ymax)]);
    xlim([0,360])
    ylabel("Norm count")

    plotEgoHistPolar(tl=tl, egodir=egoid(chk.speed>.2))
end

%% Decode pos on ww based on of rate maps
recName = "25843_1";
fname = fullfile(S.dataRoot, "navigation", "ww", recName+".mat");
tmp = load(fname);
D = tmp.Dsession;
fname = fullfile(S.dataRoot, "navigation", "of", recName+".mat");
tmp = load(fname);
unitsww=D.units.mec;
unitsof=tmp.Dsession.units.mec;
uids = intersect([unitsww.id], [unitsof.id]);
unitsof = unitsof(ismember([unitsof.id], uids));
unitsww = unitsww(ismember([unitsww.id], uids));
%%
units = unitsww;
[units.rmf] = dealArr([unitsof.rmf]);
D.units.mec = units;
timeRange = [17150, 17155];
%% Bin spikes
disp('Counting spikes...')
units = D.units.mec;
sc = reconstructSpikeCounts({units.spikeInds}, D.nt);

% Get baseline tuning
tic
vbins = restrictq(S.gv.pos_of_fine_extended, minmax(S.gv.pos_of_fine')+[-.01, .01]);

clear tuningAll
nu = numel(units);
for u = 1:nu
   tc = units(u).rmf.pos_id_shift;
   tc = units(u).rmf.pos;
   rmap = tc.z'+eps;
   rmap(isnan(rmap))=eps;
   rmap = imgaussfilt(rmap,2);
   rmap = rmap(vbins, vbins);
   tuningAll(:, u) = rmap(:);
end

% Find valid bins (4 first bins of each cycle)
chk = D.thetaChunks;
nchk = numel(chk.iStart);
sc0 = zeros(nchk, nu);
nbins = 3;
for i = 1:nchk
    [~, inds] = restrictq(chk.iStart(i)+(0:nbins), [1, D.nt]);
    sc0(i, :) = sum(sc(inds, :));
end
D.mec.sc0 = sc0;

%% Decoding
nmods = numel(D.unitAcorrClus.grid);
sigmas = [1.3, 2];
gvTile = (-1.5:0.01:1.5)';
if nmods>0
    clear res
    for m = 2
        disp("Module "+m);
        mod = D.unitAcorrClus.grid(m);
        units = D.units.mec;
        mod.unitInds = ismember([units.id]', mod.unitIds);
        units = D.units.mec(mod.unitInds);
        unique([units.acorrClu])
        spikeCounts{1} = sc(:, mod.unitInds);
        spikeCounts{2} = D.mec.sc0(:, mod.unitInds);
        vts{1} = restrictq(D.t, timeRange);
        vts{2} = restrictq(chk.tStart, timeRange);
        mod.nspikes = sum(spikeCounts{1}, 2);
        % Get baseline tuning
        clear tuning
        nu = numel(units);
        pb = ProgressBar();
        disp('Getting rate maps')
        tuning = tuningAll(:, mod.unitInds);
        mod.tuning = tuning;
        % Do single module linear dec
        clear tmp
        for s = 1:2
             vt = vts{s};
             tmp(s).peakPos = zeros(numel(vt), 2);
             tmp(s).peakProb = zeros(numel(vt), 1);
             tmp(s).torusPhase = zeros(numel(vt), 2);
             [peakPos, peakProb, torusPhase, tmp(s).info] = decodePvXcorr(tuning, spikeCounts{s}(vt, :), sigmas(s), gvTile, "regionalmax");
             tmp(s).peakPos(vt, :) = peakPos;
             tmp(s).peakProb(vt, :) = peakProb;
             tmp(s).torusPhase(vt, :) = torusPhase;
        end
        %
        mod.gridDecoding = tmp;
        mod = unwrapDec(mod, D);

        dec.decpos = mod.gridDecoding(1).peakPos;
        pslow = interp1(movmean(chk.tStart, [0,1]), mod.gridDecoding(2).peakPos, D.t, "nearest");
        dec.possm = pslow;
        dec.gridAxes = mod.gridDecoding(1).info.idealAxes;
        dec.maxprob = mod.gridDecoding(1).peakProb;
        dec.nactive = mod.nspikes;

        sweeps = chunkThetaPosSweepsModules(D, dec, mod.nspikes);
        mod.sweeps = sweeps;
        res.mods(m) = mod;
    end
end

%% Trajectory with sweeps
trng = timeRanges(2, :);
vt = restrictq(D.t, trng);
inds = find(vt);
pos = [D.x, D.y];
vspeed = D.speed>.15;
pos(~vspeed, :) = nan;
pos = gsmooth(pos, 10);

for m = 2
    nexttile(4); cla
    sweeps = res.mods(m).sweeps;
    res.mods(m).posOffset = wrapPosToGridTile(res.mods(m).posFast-res.mods(m).posSlow, res.mods(m).gridDecoding(1).info.idealAxes);
    hpfdec = (res.mods(m).posOffset);
    hpfdec = gsmooth(hpfdec, .5);
    % hpfdec=rotate2d(hpfdec, deg2rad(-20)); % Approximate map rotation
    % betwee OF and WW
    vswp = restrictq([sweeps.tStart], trng);
    allsweeps = sweeps;
    sweeps = sweeps(vswp);
      
    plot(pos(:, 1), pos(:, 2), 'Color',[1,1,1]*.7, 'LineWidth',.1)
    plot(pos(inds, 1), pos(inds, 2), 'Color',[1,1,1]*.3, 'LineWidth',.5)
    lmt = gsmooth(D.lmt.mec.pos.XA, 1);
    decpos = [D.x, D.y]+hpfdec.*res.mods(m).gridSpacingMeters;
    vt2 = D.speed>.2;
    for s = 1:numel(sweeps)
        if sweeps(s).nvalid>2
            swpinds = sweeps(s).iStart:(sweeps(s).iStop-1);
            lmt = decpos(swpinds, :);
            lmt = gsmooth(lmt, .8);
            if rem(s, 2)==0
                col = S.col_cyc_even;
            else
                col = S.col_cyc_odd;
                col = col+[0,.1, 0];
            end
            plot(lmt(:, 1), lmt(:, 2), 'Color', col)
            plot(lmt(end, 1), lmt(end, 2), '.', 'markerSize', 5, 'Color', col)
        end
    end
    scatter(pos(inds(1), 1), pos(inds(1), 2), 10, 'r')
    title(sprintf("Module %d", m))
    xlim([-1,1]); ylim([-1,1])
    axis square off
end

%% LMT examples
sweepsSetup
% S.dataRoot = ;% Specify data root
recName = "25843_1";
fname = fullfile(S.dataRoot, "navigation", "ww", recName+".mat");
tmp = load(fname);
D = tmp.Dsession;
timeRange = [17152.2, 17154.2];

%% Chunk LMT sweeps
units = D.units.mec;
sc = reconstructSpikeCounts({units.spikeInds}, D.nt);
D.sweeps.mec = chunkThetaPosSweepsLmt(D, D.lmt.mec.pos.XA, sc);

%%
units = D.units.mec;
tedg = timeRange(1):.005:(timeRange(2)+2);
tcen = tedg(2:end)-.0025;
sc = zeros(numel(tcen), numel(units));
for u = 1:numel(units)
    sc(:, u) = histcounts(units(u).spikeTimes, tedg);
    rm = units(u).rmf.pos_lmt_mec.z;
    rm = imgaussfilt(rm,3);
    tuning(:, u) = rm(:);
end
sc = gsmooth(sc, 1.5);
[prob] = decodeBayes(sc,tuning,.005,0);
dec.t = tcen;
dec.prob = reshape(prob', [96,96, numel(tcen)]);
dec.grid.x=units(u).rmf.pos_lmt_mec.bins{1};
dec.grid.y=units(u).rmf.pos_lmt_mec.bins{1};

%% 
sweeps = D.sweeps.mec;
vswp = restrictq([sweeps.tStart]', [17153.3, 17154.3]);
sweeps = sweeps(vswp);
chkinds = find(vswp);
chk = D.thetaChunks;
length = .5;
tl=tiledlayout("flow", "TileSpacing","tight");
clear hdecs;
set(gcf, 'Color', 'k')
for s = 1:6
    nexttile;
    trange = [sweeps(s).tStart+.02, sweeps(s).tStop];
    vframes = restrictq(dec.t, trange);
    inds = find(vframes);
    prob = dec.prob(:, :, vframes);
    frame = prob(:, :, 1);
    nsteps = sum(vframes);
    for i = 1:nsteps
        color =[];
        color(1,1,:) = mapcolors(i, [1,nsteps], 'cool');
        color = repmat(color, [size(frame), 1]);
        hdec(i) = imshow(color, 'XData',dec.grid.x, 'YData',dec.grid.y);
        hdec(i).AlphaData = gather(frame./max(frame(:)));
    end

    frameidx = nsteps;
    for stepidx = 1:nsteps
        frame = prob(:, :, stepidx);
        frame = frame';
        frame(frame<prctile(frame(:), 99))=0;
        hdec(stepidx).AlphaData = gather(frame./max(frame(:)));
    end
    x = D.x(chk.iStart(chkinds(s)));
    y = D.y(chk.iStart(chkinds(s)));
    id = D.id(chk.iStart(chkinds(s)));
    
    id_x = [x, x+cos(id).*length];
    id_y = [y, y+sin(id).*length];
    plot(D.x, D.y, 'Color', [1,1,1,.2], 'LineWidth',.1);
    plot(id_x(1, :), id_y(1, :), 'g','LineWidth',2); 
    clear hdec;

end

%% Plot WW and OF rate maps for example cells
exampleIds = ["2_0642", "2_0730", "2_0812"];
units = D.units.mec;
units = units(ismember([units.id], exampleIds));
rm = units(1).rmf.pos;
gv = rm.bins{1};
[covBnds, mask] = getRateMapCoverageBounds(gv, gv, rm.validBin, 2);
tl=tiledlayout(numel(units), 3);
for u = 1:numel(units)
    plotRateMap(units(u).rmf.pos, axes=nexttile), xlim([-1,1]), ylim([-1,1]); text(-1,.9,units(u).id);
    plotRateMap(units(u).rmf.pos_lmt_mec, axes=nexttile), xlim([-1,1]), ylim([-1,1]); plot(covBnds(:, 1), covBnds(:, 2), 'w')
    plotRateMap(units(u).rmf.pos_id_shift, axes=nexttile), xlim([-1,1]), ylim([-1,1]); plot(covBnds(:, 1), covBnds(:, 2), 'w')
end
title(nexttile(1), "Tracked pos")
title(nexttile(2), "LMT")
title(nexttile(3), "GLM")

%% LMT wagon wheel and LT plots
sweepsSetup
stypes = ["ww", "lt"];
i=0;
for stype = stypes
    recNames = S.("recs_"+stype);
    for recName = recNames
        fname = fullfile(S.dataRoot, "navigation", stype, recName+".mat");
        tmp = load(fname);
        D = tmp.Dsession;
        D.stype = stype;
        D.recName = recName;
        units = D.units.mec;
        sc = reconstructSpikeCounts({units.spikeInds}, D.nt);
        D.sweeps.mec = chunkThetaPosSweepsLmt(D, D.lmt.mec.pos.XA, sc);
        i = i+1;
        recs(i)=D;
        disp(recName);
    end
end
%% Plot averaged sweeps
clf
tl = tiledlayout(2,2);
inds = [1, 5];
for i = inds
    D = recs(i);
    chk = D.thetaChunks;
    chk.hd = D.hd(chk.iStart);
    chk.speed = D.speed(chk.iStart);
    egoid = circ_dist(chk.id,chk.hd);
    egoid(chk.speed<.15)=nan;
    egoid=egoid-circ_mean(egoid(~isnan(egoid)));
    plotEgoHistPolar("egodir",egoid, "tl",tl);
    if D.stype=="ww"
        title("Wagon wheel")
    else
        title("Linear track")
    end
end
fds = ["right", "left"];
cols.left = [1,0,0, .4];
cols.right = [0,0,1, .4];
stypes = ["Wagon wheel", "Linear track"];
for a = 3:4
    axs(a) = nexttile(a);
    xline(0); yline(0);
    xlim([-.15,.15]);
    ylim([-.1,.2]);
    axis square off
    if a ==3
        plot([.15,.15, .05], [0.1,.2, .2], 'k','LineWidth',2)
    end
    title((stypes(a-2)))
end 

for r = 1:numel(recs)
    D = recs(r);
    sweeps = D.sweeps.mec;
    [sweeps.nvalid] = dealArr(5); [sweeps.straight]=dealArr(5);[sweeps.straight]=dealArr(5*[sweeps.length]);
    [sweeps.hpfPosDirection] = dealArr(D.thetaChunks.id);
    dec.sweeps = sweeps;
    dec.decpos = D.lmt.mec.pos.XA;
    dec.poshpf = gsmooth(dec.decpos, S.tsm.pos)-gsmooth(dec.decpos, S.tsm.pos_slow);
    chk = D.thetaChunks;
    chk.tCen = movmean(chk.tStart, [0,1]);
    chk.iCen = round(movmean(chk.iStart, [0,1]));
    chk.hd = D.hd(chk.iCen);
    chk.speed = D.speed(chk.iCen);
    dec.chk = chk;
    tmp = plotAvgSweeps(D, dec);
    
    a = 3;
    if strcmpi(D.stype, "lt")
        a=4;
    end
    for fd = fds               
        plot(axs(a), tmp.(fd)(:, 1), tmp.(fd)(:, 2), Color=cols.(fd))
        scatter(axs(a), tmp.(fd)(end, 1), tmp.(fd)(end, 2), 20,cols.(fd)(1:3), 'MarkerFaceAlpha',.7)
    end
end


%%
%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%
function dec = unwrapDec(dec, D)
    res = dec.gridDecoding;
    chk = D.thetaChunks;

    info = res(1).info; % doesn't matter which res we use here

    posSlow = res(2).peakPos; % Get the slow-moving decoded path
    posFast = res(1).peakPos;
    pslow = nan(size(posFast));
    for i = 1:numel(chk.iStart)-1

        inds = chk.iStart(i):(chk.iStart(i+1)-1);
        pslow(inds, :) = posSlow(i, :)+zeros(numel(inds), 1);
    end
    posSlow = pslow;
    
    dec.posSlow = posSlow;
    dec.posFast = posFast;
end