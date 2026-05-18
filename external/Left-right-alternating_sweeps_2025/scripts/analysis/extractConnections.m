function [outputArg1,outputArg2] = extractConnections(p)
%EXTRACTCONNECTIONS Extract candidate connections from xcorrdata

arguments
    p.sourceFolder = [];
    p.recNames = [];
    p.outFolder = [];
    p.fnames = [];
    p.min_spikes = 1000;
    p.min_udist =0;
    p.min_peakZ = 1;
    p.min_peakExcess = 20;
    p.lumpGrid = true;
    
end
S = SweepsSettings;
%%
if isempty(p.recNames)
    p.recNames = S.recs_of_mec(:);
end
if isempty(p.outFolder)
    p.outFolder = fullfile(p.sourceFolder, "filtered");
    if~isfolder(p.outFolder)
        mkdir(p.outFolder)
    end
end
datafld = fullfile(S.dataRoot, "navigation", "of");

nrecs = numel(p.recNames);
%% Save an info file with general settings

for r = 1:nrecs
    %%
    rec = p.recNames(r);
    tic
    
    fname = fullfile(p.sourceFolder, rec);
    fprintf("Loading conndat...")
    connDat = load(fname);
    fprintf("Done\n")

    % Select valid xcorrs
    fname = fullfile(datafld, rec+".mat");
    tmp = load(fname);
    D = tmp.Dsession;
    Us = D.units.mec;
    units = Us;
    nu = numel(Us);
    assert(isequal([Us.id]', connDat.unitIds));

    %% Load "lite" units which we can save easily
    if p.lumpGrid
        gridIdx = ismember([Us.cellType], ["bursty", "nonbursty", "prospective"]);
        [Us.cellTypeAll] = dealArr([Us.cellType]);
        [Us(gridIdx).cellType]= dealArr("grid");
    end
    
    for u = 1:nu
        Us(u).burstScore = burstScore(Us(u).spikeTimes);
    end
    Us = rmfield(Us, ["spikeTimes", "spikeInds"]);
    [Us.nspkFull] = dealArr(connDat.conns.nspk_i);
    %%
    C = connDat.conns;
    C.nspk_i = repmat(C.nspk_i, 1, nu);
    C.nspk_j = repmat(C.nspk_j, nu, 1);

    hasDistantUnits = isnan(C.uDist) | (C.uDist >= p.min_udist); % automatically excludes same-cell pairs
    hasEnoughEvents = C.nspk_i>p.min_spikes & C.nspk_i>p.min_spikes;
    isValidXcorr = hasDistantUnits & hasEnoughEvents;
    [iAll, jAll] = find(isValidXcorr);

    hasPeakMn = C.pkExcess >= p.min_peakExcess; % excess peak is at least X spikes
    hasPeakZ = C.pkvalZ >= p.min_peakZ; % excess peak is at least X times the standard deviation of the *excess* xcorr
    isValidPeak = hasPeakMn & hasPeakZ;

    isConnected = isValidXcorr & isValidPeak;
    ipair = find(isConnected);
    [i, j] = find(isConnected);
    % Add totalCntAll somewhere here
    %% Save connections
    res = structfun(@(x) x(ipair), C, "uni", 0);
    res.unitIds = connDat.unitIds([i, j]);
%     res.unitIds = connDat.unitIds([iAll, j]);
    res.units = Us;
    
    res.recName = p.recNames(r);
    res.recDuration = connDat.recinfo.recDuration;

    fnOut = fullfile(p.outFolder, "res_"+p.recNames(r));
    save(fnOut, "-struct", "res", "-v7", "-nocompression");
    
    %% Save valid xcorrs
    xc = (connDat.xcorrs);
    res.isValidXcorr = isValidXcorr;
    res.isConnected = isConnected;
    res.xcorrs = structfun(@(x) single(x(:, ipair)), xc, "uni", 0);
    
    %%
    fnOut = fullfile(p.outFolder, "res_heavy_"+p.recNames(r));
    save(fnOut, "-struct", "res", "-v7", "-nocompression");
    %% ------------------------------------------------------------------------------------------
    % Save the total number of outgoing pairs from each celltype (w and wo
    nPairs = struct();
    nPairs.all = sum(isValidXcorr(:));
    pretype = [Us(iAll).cellType];
    posttype = [Us(jAll).cellType];
    types = ["id", "conjunctive", "bursty", "nonbursty", "prospective"];
    
    % All combs
    if p.lumpGrid
    types = ["id", "conjunctive", "grid"];
    end
    ntypes = numel(types);
    combinds = combvec(1:ntypes, 1:ntypes)';
    combs = fliplr(types(combinds));
    for c = 1:size(combs, 1)
        combname = combs(c, 1)+"_"+combs(c, 2);
        nPairs.(combname) = sum(pretype(:)==combs(c, 1) & posttype(:)==combs(c, 2)); % n pairs from type a to type b
        nPairs.(combs(c, 1)+"_other") = sum(pretype(:)==combs(c, 1) & posttype(:)~=combs(c, 1)); %total outgoing excl recurrent
%         nPairs.(combs(c, 1)+"other_") = sum(pretype(:)==combs(c, 1) & posttype(:)~=combs(c, 1)); %total outgoing excl recurrent
    end
    
    fnOut = fullfile(p.outFolder, "nPairs_"+p.recNames(r)); 
    save(fnOut, "-struct", "nPairs", "-v7", "-nocompression");
    %------------------------------------------------------------------------------
    %%
    % Save tuning relationships with conn stats (only for somewhat connected pairs
%     P = p;
    tuning = struct();
    
    pretype = [Us(i).cellType];
    posttype = [Us(j).cellType];
    unitIds = connDat.unitIds([i, j]);
    linearInds = 1:numel(pretype);
    % Get prefdirs for directional pairings
    meandir = @(units, vname) arrayfun(@(u) circ_mean(D.gv.hd, u.rmf.(vname).z), units);
    predirs = meandir(units(i), "id");
    postdirs = meandir(units(j), "id");
    %
    dirtypes = ["id", "conjunctive"];
    ntypes = numel(dirtypes);
    combinds = combvec(1:ntypes, 1:ntypes)';
    combs = fliplr(dirtypes(combinds));
    tmp = struct();
    for c = 1:size(combs, 1)
        combname = combs(c, 1)+"_"+combs(c, 2);
        pairinds = find(pretype(:)==combs(c, 1) & posttype(:)==combs(c, 2)); % vlaid combo
        npairs = numel(pairinds);
        for pp = 1:npairs% Iterate through pairs and get tuning dirs
            
            tmp.(combname).predir(pp, 1) = predirs(pairinds(pp));
            tmp.(combname).postdir(pp, 1) = postdirs(pairinds(pp)); 
            tmp.(combname).unitIds(pp, :) = unitIds(pairinds(pp),:);
            tmp.(combname).linearInds(pp, 1) = linearInds(pairinds(pp));
        end
    end
    tuning.prefdir = tmp;
    %
%     clf
    vmod = ([units(i).acorrCluId] == [units(j).acorrCluId])';
    posttype = [Us(j).cellType];
    gridIdx = ismember(posttype, ["bursty", "nonbursty", "prospective"]);
    posttype(gridIdx) ="grid";
    gridtypes = ["conjunctive", "grid"];
    ntypes = numel(gridtypes);
    combinds = combvec(1:ntypes, 1:ntypes)';
    combs = fliplr(gridtypes(combinds));
    tmp = struct();
    for c = 1:size(combs, 1)
        combname = combs(c, 1)+"_"+combs(c, 2);
        disp(combname)
        vcombo = (pretype(:)==combs(c, 1) & posttype(:)==combs(c, 2)); % vlaid combo
        % Only consider same module pairings
        
        pairinds = find(vcombo & vmod);
        npairs = numel(pairinds);
        for pp = 1:npairs% Iterate through pairs and get tuning dirs
            tmp.(combname).predir(pp, 1) = predirs(pairinds(pp));
            unit1 = units(i(pairinds(pp))); unit2 = units(j(pairinds(pp)));
            [phaseoffset, shiftdir, dist] = getPhaseOffset(unit1, unit2);
            tmp.(combname).phaseoffset(pp, :) = phaseoffset;
            tmp.(combname).phaseoffsetdir(pp, 1) = shiftdir; 
            tmp.(combname).phasedist(pp, 1) = dist;  
            tmp.(combname).unitIds(pp, :) = unitIds(pairinds(pp),:);
            tmp.(combname).linearInds(pp, 1) = linearInds(pairinds(pp));
        end
    end
    tuning.gridphase = tmp;
    %%
    fnOut = fullfile(p.outFolder, "tuning_"+p.recNames(r));
    %
    save(fnOut, "-struct", "tuning", "-v7", "-nocompression");
end
fprintf("Done with %s", p.recNames(r))    
toc
clearvars -except r p S nrecs
end


function  [phaseoffset, shiftdir, shiftdist,spacing] = getPhaseOffset(unit1, unit2)
S = SweepsSettings;
vbins = restrictq(S.gv.pos_of_fine_extended, [min(S.gv.pos_of_fine)-.01, max(S.gv.pos_of_fine)+.01]);
% vbins(:)=true;
rms1 = unit1.rmf.pos_id_shift;
rms2 = unit2.rmf.pos_id_shift;

% FIX      
z1 = rms1.z(vbins, vbins); z1 = regionfill(z1, isnan(z1));%
z2 = rms2.z(vbins, vbins); z2 = regionfill(z2, isnan(z2));%
%         
% xc = xcorr2(z2', z1');
xc = normxcorr2_general(z2', z1');
bw = imregionalmax(xc);
s = regionprops(bw);
peaks = cat(1, s.Centroid);
origo = size(xc)./2;
pkrel = peaks-origo;
if ~isempty(peaks)
    [~,minIdx] = min(hypot(pkrel(:, 1), pkrel(:, 2)));
end
nearestPk = peaks(minIdx, :)-origo;

[shiftdir, shiftdist] = cart2pol(nearestPk(1), nearestPk(2));
phaseoffset = nearestPk;
end
