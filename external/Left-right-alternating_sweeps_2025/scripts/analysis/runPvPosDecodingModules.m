function [rez] = runPvPosDecodingModules(p)
%RUNLINEARDECPOSRECS Runs corr-based decoder on a set of recordings and
%saves results in specified folder
arguments
    p.recs = [];
    p.fld = [];
    p.process = 1;
    p.load = 0;
    p.save = 0;
    p.brainregion = "mec"
    p.sessionType = "of"
end

S = sweepsSettings;
if isempty(p.fld)
    fld = fullfile(S.dataRoot, "code", "results", "pv_dec", p.brainregion, "modules");
else
    fld = p.fld;
end

%%
if isempty(p.recs)
    recs = S.recs_of_mec';
else
    recs = p.recs;
end
nrecs = size(recs, 1);
useGpu = 1;
datafld = fullfile(S.dataRoot, "navigation", p.sessionType);
if strcmpi(p.sessionType, "sleep")
    datafld = fullfile(S.dataRoot, "sleep");
end

sigmas = [1.3, 2];
gvTile = (-1.5:0.01:1.5)';
%% Iterate through recs
if p.process
    for r = 2:nrecs
        %% Load rec
        rec = recs(r);
        disp('Loading '+rec+'...')
        fname = fullfile(datafld, rec+".mat");
        tmp = load(fname);
        D = tmp.Dsession;

        chk = D.thetaChunks;
        chk = rmfield(chk, ["L", "P"]);
        chk.tCen = movmean(chk.tStart, [0,1]);
        chk.iCen = round(movmean(chk.iStart, [0,1]));
        chk.hd = D.hd(chk.iCen);
        chk.speed = D.speed(chk.iCen);
        %
        units = D.units.(p.brainregion);
    
        %% Bin spikes
        disp('Counting spikes...')
        % for u = 1:numel(units)
        %     units(u).sc = gather(binSpikes(units(u).spikeTimes, [], D.t, useGpu));
        % end
        
        sc = reconstructSpikeCounts({units.spikeInds}, D.nt);
    %     %%

    
        % Get baseline tuning
        tic
        vbins = restrictq(S.gv.pos_of_fine_extended, minmax(S.gv.pos_of_fine')+[-.01, .01]);
    
        clear tuningAll
        nu = numel(units);
        for u = 1:nu
           tc = units(u).rmf.pos;
           rmap = tc.z'+eps;
           rmap(isnan(rmap))=eps;
           rmap = imgaussfilt(rmap,1);
           rmap = rmap(vbins, vbins);
           tuningAll(:, u) = rmap(:);
        end

        % Find valid bins (4 first bins of each cycle)
        nchk = numel(chk.iStart);
        sc0 = zeros(nchk, nu);
        nbins = 3;
        for i = 1:nchk
            [~, inds] = restrictq(chk.iStart(i)+(0:nbins), [1, D.nt]);
            sc0(i, :) = sum(sc(inds, :));
        end
        D.(p.brainregion).sc0 = sc0;
        %% Decoding
        nmods = numel(D.unitAcorrClus.grid);
        if nmods>0
            clear res
            for m = 1:nmods
                %%
                disp("Module "+m);
                mod = D.unitAcorrClus.grid(m);
                units = D.units.mec(mod.unitInds);
                spikeCounts{1} = sc(:, mod.unitInds);
                spikeCounts{2} = D.mec.sc0(:, mod.unitInds);
                mod.nspikes = sum(spikeCounts{1}, 2);
                % Get baseline tuning
                clear tuning
                nu = numel(units);
                pb = ProgressBar();
                disp('Getting rate maps')
                tuning = tuningAll(:, mod.unitInds);
                mod.tuning = tuning;
                % Do single module linear dec
                dbstop if error
                clear tmp
                for s = 1:2
                     [tmp(s).peakPos, tmp(s).peakProb, tmp(s).torusPhase, tmp(s).info] = decodePvXcorr(tuning, spikeCounts{s}, sigmas(s), gvTile, "regionalmax");
                end
                %%
                mod.gridDecoding = tmp;
                %
                % Unwrap
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
                % Store
            end
        
            res.recName = rec;
            res.t = D.t;
            res.x = D.x;
            res.y = D.y;
            res.speed = D.speed;
            res.hd = D.hd;
            res.id = D.id;
            res.chk = chk;
            res.theta = D.theta;
                
            % Save results (log likelihood, dec stuff (peak, prob, etc), ratemaps) 
            
            % store
            rez(r) = res;
            
            disp("Saving..")
            if p.save
                if ~isfolder(fld), mkdir(fld); end
                subfld = fullfile(fld, rec);
                save(subfld, "-struct", "res", "-v7", "-nocompression");
            end
            clear D dec res
        end
        disp("Done with"+rec)
    end
end


%% Load results
if p.load
    clear rez
    nrecs = size(recs, 1);
    for r = 1:nrecs
        rec = recs(r);
        disp('Loading '+rec+'...')
        fname = fullfile(fld, rec+".mat");
        rez(r) = load(fname);
    end
end

end
%%%%%%%%%%%%%%%% functions

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