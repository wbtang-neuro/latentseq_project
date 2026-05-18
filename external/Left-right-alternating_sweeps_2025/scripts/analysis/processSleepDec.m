function [outputArg1,outputArg2] = processSleepDec()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
sweepsSetup
recs = findSweepsRecs("mec", "sleep_best");
idfld = S.filepath("id_grid_decoding_corr\id");
decfld = S.filepath("id_grid_decoding_v3\v2");
savefld = S.filepath("sleepData\withsweeps2");
%%
for r = 1:height(recs)
    %%
    recName = getRecName(recs(r, :));
    [D, Dfull, Dtmp, dec] = loadFullRecSweeps(recName, "loadDec",1);
    res = rmfield(Dfull, "units");
    res.recName = recName;
    
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
    sc = [units.spikeCounts];
    res.poprate = (full(sum(sc, 2)));
    res.idspk = (full(sum(sc(:,isdir), 2)));
    
    %% ID
    tmp = S.load(fullfile(decfld, recName+".mat"));
    res.id_lmt = gather(circ_mean(tmp.idBins', tmp.id')');
    res.idprob_lmt = max(tmp.id, [], 2)';
    
    % ID corr
    tmpid = S.load(fullfile(idfld, recName+".mat"));
    res.id_corr = circ_mean(tmpid.idBins, tmpid.id')';
    res.idprob_corr = max(tmpid.id')';
    %%
    %%% POP BURSTS CHUNKS
    sc = gsmooth(res.idspk, 2);
    % find chunks (centered around peaks and extending up to nearest trough)
    disp(res.recName)
    disp("chunking..")
    [pks, iloc] = findpeaks(sc);
    [~, ilow] = findpeaks(-sc);
    chk = struct;
    chk.icen = iloc;
    chk.iStart = chk.icen;
    chk.iStop = chk.icen;
    chk.vchk = false(size(iloc));
    pb = ProgressBar();
    for c = 5:(numel(iloc)-5)
        icen = chk.icen(c);
        chk.iStart(c) = max(icen-5, max(ilow(ilow<icen)));
        chk.iStop(c) = min(icen+5, min(ilow(ilow>icen)));
        chk.vchk(c) = true;
        pb.update(c/numel(iloc))
    end
    %
    res.thetaChunks = chk;
    

    %% Grid decoding
    griddec = tmp;
    nmod = numel(griddec.gridDecoding);
    for m = 1:nmod
        if ~isempty(griddec.gridDecoding{m})
            res.gridDecoding{m} = rmfield(griddec.gridDecoding{m}, "spikeCounts");
            
            dec.decpos = res.gridDecoding{m}(1).peakPos;
            dec.possm = res.gridDecoding{m}(2).peakPos;
            dec.maxprob = res.gridDecoding{m}(1).peakProb;
            dec.gridAxes = res.gridDecoding{m}(1).info.idealAxes;
            dec.acorrRotation = res.gridDecoding{m}(1).info.acorrRotation;
            sweeps = chunkThetaPosSweepsModSleep(res, dec,gsmooth(res.idspk, 2));
            
            res.mods{m}.sweeps = sweeps;
        end
    end

    %% Remove unnecessary stuff
    fdskeep = ["poprate", "idspk", "states", "id_lmt", "t", "id_corr", "recName", "gridDecoding", "mods", "hd", "x", "y", "sessions", "sleepData", "thetaChunks"];
    fds = fieldnamesstr(res);
    fdsrm = fds(~ismember(fds,fdskeep));
    res = rmfield(res, fdsrm);
    %% Save
    S.save(fullfile(savefld, recName+".mat"), res, 1)
    disp(recName)
end
end

