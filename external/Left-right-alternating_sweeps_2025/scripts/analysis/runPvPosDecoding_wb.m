function [res] = runPvPosDecoding_wb(p)
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
    p.selectedUnits = [];
end

S = sweepsSettings;
if isempty(p.fld)
    fld = fullfile(S.dataRoot, "code", "results", "pv_dec", p.brainregion, "pos");
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
useGpu = 0;
% datafld = fullfile(S.dataRoot, "navigation", p.sessionType);
datafld = S.dataRoot;
if strcmpi(p.sessionType, "sleep")
    datafld = fullfile(S.dataRoot, "sleep");
end
%% Iterate through recs
if p.process
    for r = 1:nrecs
        % Load rec
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
        if ~isempty(p.selectedUnits)
            units = units(p.selectedUnits);
        end
    
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
    
        clear tuning
        nu = numel(units);
        for u = 1:nu
           tc = units(u).rmf.pos;
           rmap = tc.z'+eps;
           rmap(isnan(rmap))=eps;
           rmap = imgaussfilt(rmap,1);
           rmap = rmap(vbins, vbins);
           tuning(:, u) = rmap(:);
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
        %%
        tic
        disp("dec fast")
        % prob = decodePv(tuning=tuning, spikeCounts=sc, smoothSpikes=S.tsm.id, chunked=true);
        prob = decodeBayes(sc, tuning, 0.01,S.tsm.id);
        toc

        tic
        disp("dec slow")
        % probsm = decodePv(tuning=tuning, spikeCounts=sc0, smoothSpikes=1.7, chunked = true);
        probsm = decodeBayes(sc0, tuning, 0.04,1.7);
        toc
           %%
        dec = processDec(prob');
        decsm = processDec(probsm');

        %% Get some shuffled probs
        sc_shuff = sc;
        nt = size(sc_shuff, 1);
        nu = size(sc_shuff, 2);
        for u = 1:nu
            sc_shuff(:, u)  = circshift(sc_shuff(:, u), randi(nt));
        end
        tmax = min([1e4, nt]);
        probshuff = decodePv(spikeCounts=sc_shuff(1:tmax, :), tuning=tuning, smoothSpikes=1, chunked = true);
        dec.shuff99 = prctile(probshuff(:), 99);
        %%
        decsm.t = chk.tStart + .5*D.dt*(nbins+1);
        decsm.decpos = gsmooth(decsm.decpos, .5);
        dec.possm = interp1(decsm.t, decsm.decpos, D.t, "linear");
    
        dec.thetaChunks = chk;
        
        active = sc; active(sc>0) = 1;
        dec.nactive = sum(active, 2);
        dec.nspikes = sum(sc, 2);
        err = dec.possm-[D.x, D.y];
        dec.err = hypot(err(:, 1), err(:, 2));
        dec.(p.brainregion) = p.brainregion;
        sweeps = chunkThetaPosSweeps(D, dec, dec.nspikes, 0);
        %%
        dec.sweeps=sweeps;
        dec.chk = chk;
        % store
            
        % Save results (log likelihood, dec stuff (peak, prob, etc), ratemaps) 
        disp("Done with"+rec)
        %% store
        res(r) = dec;
        
        disp("Saving..")
        if p.save
            if ~isfolder(fld), mkdir(fld); end
            subfld = fullfile(fld, rec);
            save(subfld, "-struct", "dec", "-v7", "-nocompression");
        end
        clear D dec
        disp('Done')
    end
end


%% Load results
if p.load
    clear res
    nrecs = size(recs, 1);
    for r = 1:nrecs
        rec = recs(r);
        disp('Loading '+rec+'...')
        fname = fullfile(fld, rec+".mat");
        res(r) = load(fname);
    end
end
%%
% if p.save
%     S.save(fullfile(decfld, "resAll_"+p.brainregion), data, true);
% end
end
%%%%%%%%%%%%%%%% functions
function dec = processDec(prob)
fast = 0; % fix
calcSpread = 0;
S = sweepsSettings;
[xi, yi] = meshgrid(S.gv.pos_of_fine);
% vbins = restrictq(S.gv.pos_of_fine_extended, minmax(S.gv.pos_of_fine')+[-.03, .03]);
% [xi, yi] = meshgrid(S.gv.pos_of_fine_extended(vbins));
xx = xi(:); yy = yi(:); % This is where the transpose happens!!
gv = S.gv.pos_of_fine;
% gv = (S.gv.pos_of_fine_extended(vbins));
npos = numel(S.gv.pos_of_fine);
% npos = sum(vbins);
[maxprob, imx] = max(prob);
dec.maxprob = maxprob;
dec.posmax(:, 1) = xx(imx); dec.posmax(:, 2)= yy(imx);
% dec.prc80(:, 1) = prctile(prob, 80);
% dec.prc90(:, 1) = prctile(prob, 90);
% dec.prc99(:, 1) = prctile(prob, 99);
dec.pkmu(:, 1) = mean(prob);
dec.pkstd(:, 1) = std(prob);

if fast
    dec.decpos = dec.posmax;
else
    spreadthresh = (prctile(prob, 99));
    thresh = (prctile(prob, 99));
    maxdist = hypot(xx-dec.posmax(:, 1)', yy-dec.posmax(:, 2)');
    if calcSpread
        % weight dists by decprob
        posprob = prob; posprob(posprob<spreadthresh)=0;
        dec.spread = sqrt(sum(maxdist.*posprob)./sum(posprob));
    end
    isbad = prob<thresh & maxdist>.1;
    prob(isbad) = 0;
    sumweights = sum(prob, 'omitnan');
    nt = size(prob, 2);
    prob = reshape(prob, [npos, npos, nt]);
    xweights = squeeze(sum(prob, 1, 'omitnan'));
    yweights = squeeze(sum(prob, 2, 'omitnan'));
    
    dec.decpos(:, 1) = sum(gv.*xweights, 'omitnan')./sumweights;
    dec.decpos(:, 2) = sum(gv.*yweights, 'omitnan')./sumweights;
end

xy = dec.decpos-dec.posmax;
dec.maxoffset = hypot(xy(:, 1),xy(:, 2));
end
