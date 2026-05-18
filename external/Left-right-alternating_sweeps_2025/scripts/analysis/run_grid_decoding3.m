% 27.10.24 extend rate maps and adjust smoothing 
sweepsSetup;

recs = findSweepsRecs("mec", "sleep_best");
nrecs = size(recs, 1);
%%
% Define a meshgrid of template coordinates for the unit phase tile in 2D.
% This spans the range +/- 1.5, so that when plotting we'll see all of the
% central tile, plus part of the 6 neighboring tiles.
gvGridTile = (-1.5:0.01:1.5)';

% Create idealized hexagonal grid axes
idealGridAngles = [0, pi/3, 2*pi/3];
idealGridAxes = zeros(3, 2);
for a = 1:3
    idealGridAxes(a, :) = rotate2d([1, 0], idealGridAngles(a));
end
%% Load each rec, and decode grid tile position for each module
usegpu=1;

force = 0;

sigmaId = 1;

clear allData

for r = nrecs:-1:1

    rec = recs(r, :);
    recName = getRecName(rec);
    fprintf("%s ...\n", recName);

    clear placeholder
    fnOutId = S.filepath("id_grid_decoding_v3", "v2", recName + ".mat");
    if exist(fnOutId, "file") && ~force
        fprintf("File '%s' already exists. Skipping.\n", fnOutId);
        continue;
    end

    placeholder = TempFile(fnOutId);

    try
        sType = findOfSessions(rec, true);
    catch
        warning("No OF session type found. Skipping.");
        continue;
    end

    D = loadBasicSweepsData(recName, sType, "units", "complete");
    Us = D.units.mec;

    ksd = parseKsDirname(rec.ksdir_1);
    fullTimeRange = ksd.timeRange;
    if isempty(fullTimeRange)
        fullTimeRange = "all";
    end
    
    Dfull = loadSweepsDataset(recName, "all", ...
        "validUnits", "all", ...
        "filterTime", false, ...
        "customTimeRange", fullTimeRange, ...
        "loadWaveforms", false);
    %%% decode ID
    % Recalculate mean firing rates
    for u = 1:numel(Dfull.units)
        Dfull.units(u).meanRate = numel(Dfull.units(u).spikeTimes) / (Dfull.nt * Dfull.dt);
    end

    % Load LMT
    fn = S.filepath("lmt_fits", "mec", "id+pos", sType, recName + ".mat");
    lmtMdl = S.load(fn, "mdl", true);
    mPos = lmtMdl.getModel("pos");
    mId = lmtMdl.getModel("id");

    % Create LMT copy based on full-recording data
    lmtArgs = getLmtArgs("id+pos");
    lmtArgs = struct2arglist(lmtArgs);
    UsId = findUnits(Dfull.units, lmtMdl.neuronIds);
    lmtMdlFull = createLMT(Dfull, UsId, lmtArgs{:}, ...
        "chunkLength", Dfull.nt, ...
        "useGpu", false);
    lmtMdlFull.copyTuning(lmtMdl, false);

    % Construct "thetaChunks" struct in the same way we do for the "basic"
    % datasets
    [cycdat, timedat] = decodeThetaChunkedId(lmtMdlFull, false);
    cycdat.tStart = Dfull.t(cycdat.iStart);
    cycdat.tStartInterp = interp1(Dfull.t, cycdat.iStartInterp);

    % Decode ID and save
    Lraw = decodeFromLmt(mId, UsId, Dfull.t, [], [], false);
    [~, Pid] = normalizeLL(Lraw, sigmaId);
    clear Lraw

    fprintf("Saving ID data...\n");
    data = struct( ...
        "t", Dfull.t, ...
        "id", single(Pid), ...
        "idBins", mId.Fgridv{1}, ...
        "thetaChunks", cycdat);
    clear Pid cycdat
    %%%%
    % Match up units
    units = D.units.mec;
    unitsfull = Dfull.units;
    for u = 1:numel(units)
        unit = unitsfull([unitsfull.id]==units(u).id & [unitsfull.location]=="mec");
        units(u).spikeTimes = unit.spikeTimes;
        units(u).spikeCounts = (unit.spikeCounts);
    end
%     mvl = @(units, vname) arrayfun(@(u) u.rmf.(vname).mvl, units);
%     isdir = mvl(units, "hd")>.2 & mvl(units, "theta")>.3;
%     idunits = units(isdir);
%     nu = numel(idunits);
%     tuning = nan(60, nu);
%     %
%     pb = ProgressBar();
%     for u = 1:nu
%        tc = angularTuningCurve(D.t, D.hd, idunits(u).spikeInds, "smooth", 2, "nbins", 60);
%        tuning(:, u) = vec(tc.z);
%        pb.update(u./nu)
%     end
%     
%     [prob] = linearDecBatch(tuning=tuning, spikeCounts=full([idunits.spikeCounts]), ...
%             chunkSize=3e5, smoothing=0, vt = "all", smoothSpikes = 1, normalize="meandivide");
% 
%     fprintf("Saving ID data...\n");
%     data = struct( ...
%         "t", Dfull.t, ...
%         "id", single(prob), ...
%         "sc_id", sum(full([idunits.spikeCounts]), 2),...
%         "idBins", S.gv.angular);
%     clear prob

    % Now do the single-module grid stuff
    gmods = D.unitAcorrClus.grid;
    fprintf("%-15s: n=%u gmods\n", recName, numel(gmods));
    %
    % GRID DECODING
    vbins = restrictq(S.gv.pos_of_fine_extended, [min(S.gv.pos_of_fine)-.2, max(S.gv.pos_of_fine)+.2]);
    clear tuningAll
    nu = numel(units);
    for u = 1:nu
       rmap = units(u).rmf.pos_lmt_mec;
       z = rmap.z(vbins, vbins);
       z = z';
       z(isnan(z)) = 0;
       z = imgaussfilt(z, 1);
       tuningAll(:, u) = vec(z);
       % pb.update(u./nu)
    end

    sigmaGrid = [1, 20];

    for g = 1:numel(gmods)

        gmod = D.unitAcorrClus.grid(g);
        gunits = units(gmod.unitInds);
        tuning = tuningAll(:, gmod.unitInds);
        res = struct;
        if numel(gunits)>40
            clear res
            for s = 1:numel(sigmaGrid)
                fprintf("\tgmod #%u, sigma=%.1f ...\n", g, sigmaGrid(s));
                fn = sprintf("grid_decoding_%s_g%u.mat", recName, g);
                fn = S.filepath(fn);
                if isfile(fn) && ~force
                    fprintf("File '%s' already exists: skipping\n", fnOutId);
                    continue
                end
    
                r = struct();
                
                [r.peakPos, r.peakProb, r.torusPhase, r.info] = decXcorrBayes(tuning, full([gunits.spikeCounts]), sigmaGrid(s), gvGridTile, "regionalmax");
 
                r.torusPhase = 0;
                r.peakPos = single(r.peakPos);
                r.peakProb = single(r.peakProb);
                r.sigma = sigmaGrid(s);
                r.spikeCounts = sum(full([gunits.spikeCounts]),2);
                r.nu = numel(gunits);
                res(s) = r;
                clear r
            end
            data.gridDecoding{g} = res;
        end

        
        clear res Dtmp
    end

    %

    placeholder.deactivate();
    S.save(fnOutId, data, true);
    clear Dfull D

end

%%
sweepsSetup
recs = findSweepsRecs("mec", "sleep");
idfld = S.filepath("id_grid_decoding_corr\id");
%%
for r = 7:height(recs)
    recName = getRecName(recs(r, :));
    [D, Dfull, Dtmp, dec] = loadFullRecSweeps(recName);
    res = rmfield(Dfull, "units");
    
    % sleep stuff
    sleepTimes = Dfull.sleepData.times;
    ofTimes = Dfull.sessions(strcmpi([Dfull.sessions.type], "open_field")).validTimeRanges;
    sleepTimes.of = ofTimes;
    
    states = struct();
    states.run = restrictq(Dfull.t, ofTimes) & D.speed>.1;
    states.rem = restrictq(Dfull.t, sleepTimes.rem);
    states.sws = restrictq(Dfull.t, sleepTimes.sws);
    res.states = states;

    % pop activity
    units = D.units.mec;
    ctype = classifyGridCellType(units);
    [units.cellType] = dealArr(ctype);
    isdir = [units.cellType]=="id";
    sc = [units(isdir).spikeCounts];
    res.idspk = (full(sum(sc, 2)));
    
    % ID lmt
    res.id_lmt = circ_mean(dec.idBins', dec.id')';
    res.idprob_lmt = max(dec.id, [], 2);
    
    % ID corr
    tmp = S.load(fullfile(idfld, recName+".mat"));
    res.id_corr = circ_mean(tmp.idBins, tmp.id')';
    res.idprob_corr = max(tmp.id')';

    % grid decoding
    res.thetaChunks = dec.thetaChunks;
    res.gridDecoding = dec.gridDecoding;
    res.mods = dec.mods;
    
    res.recName = recName;
    resAll(r) = res;
    clear D Dfull Dtmp dec tmp res
    disp(recName)
end


%% Save
for r = 1:numel(recs)
    res = recs(r);
    disp(res.recName)
    S.save(S.filepath("sleepData", "wsweeps", res.recName+".mat"), res)
    disp("done")
end