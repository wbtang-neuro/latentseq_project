function [res] = runPosDecoding_wb_cheeseboard(p)
%runPosDecoding_wb_cheeseboard, runs Bayesian decoder on a set of recordings
%Modified from runPvPosDecoding.m (Vollan et al., 2024)
arguments
    p.rec = [];
    p.fld = [];
    p.save = 0;
    p.validbins = [];
    p.validbins_posxy = [];
    p.selectedUnits = [];
end

S = sweepsSettings;
if isempty(p.fld)
    fld = fullfile(S.dataRoot, "code", "results", "pv_dec", "pos");
else
    fld = p.fld;
end
%% loading file
datafld = S.dataRoot;
rec = p.rec;

disp('Loading '+rec+'...')
fname = fullfile(datafld, rec+".mat");
tmp = load(fname);
D = tmp.Dsession;
%% theta cycles
chk = D.thetaChunks;
chk.tCen = movmean(chk.tStart, [0,1]);
chk.iCen = round(movmean(chk.iStart, [0,1]));
chk.hd = D.hd(chk.iCen);
chk.speed = D.speed(chk.iCen);
%% select units and get spikes
units = D.units;
if ~isempty(p.selectedUnits)
    units = units(p.selectedUnits);
end
%% Bin spikes
disp('Counting spikes...')
sc = reconstructSpikeCounts({units.spikeInds}, D.nt);
for u = 1:numel(units)
    sc(:,u) =locsmooth(sc(:,u),1/0.01,0.05);
end
%% Get baseline tuning
tic
clear tuning
nu = numel(units);
vbins = p.validbins;
for u = 1:nu
    tc = units(u).rmf;
    rmap = tc.z+eps;
    rmap(isnan(rmap))=eps;
    rmap = imgaussfilt(rmap,1);
    tuning(:, u) = rmap(vbins);
end
%% Find valid bins (4 first bins of each cycle)
nchk = numel(chk.iStart);
sc0 = zeros(nchk, nu);
nbins = 3;
for i = 1:nchk
    [~, inds] = restrictq(chk.iStart(i)+(0:nbins), [1, D.nt]);
    sc0(i, :) = sum(sc(inds, :));
end
%% Bayesian decoding
tic
disp("dec fast")
prob = decodeBayes(sc, tuning, 0.01,1);
toc

tic
disp("dec slow")
probsm = decodeBayes(sc0, tuning, 0.04,1.7);
toc

% assign decoded postion based on posterior probability
dec = processDec(prob',p.validbins_posxy);
decsm = processDec(probsm',p.validbins_posxy);
%% Get some shuffled probs
sc_shuff = sc;
nt = size(sc_shuff, 1);
nu = size(sc_shuff, 2);
for u = 1:nu
    sc_shuff(:, u)  = circshift(sc_shuff(:, u), randi(nt));
end
tmax = min([1e4, nt]);
probshuff = decodeBayes(sc_shuff(1:tmax, :), tuning, 0.01,1);
dec.shuff99 = prctile(probshuff(:), 99);

%% Distance between decoded (slow) and actual positions
decsm.t = chk.tStart + .5*D.dt*(nbins+1);
decsm.decpos = gsmooth(decsm.decpos, .5);
dec.possm = interp1(decsm.t, decsm.decpos, D.t, "linear");

dec.thetaChunks = chk;

active = sc; active(sc>0) = 1;
dec.nactive = sum(active, 2);
dec.nspikes = sum(sc, 2);
err = dec.possm-[D.x, D.y]; 
dec.err = hypot(err(:, 1), err(:, 2)); % distance to actual position
%% get sweep info
sweeps = chunkThetaPosSweeps_wb(D, dec, dec.nspikes);

dec.sweeps=sweeps;
dec.chk = chk;
%% Save results (log likelihood, dec stuff (peak, prob, etc), ratemaps)
disp("Done with"+rec)
%% store
res = dec;

disp("Saving..")
if p.save
    if ~isfolder(fld), mkdir(fld); end
    subfld = fullfile(fld, rec);
    save(subfld, "-struct", "dec", "-v7", "-nocompression");
end
clear D dec
disp('Done')
end

%%%%%%%%%%%%%%%% functions
function dec = processDec(prob,bins_xy)
fast = 0; % fix
calcSpread = 0;
% gv = (S.gv.pos_of_fine_extended(vbins));
npos = numel(bins_xy(:,1));
[maxprob, imx] = max(prob);
dec.maxprob = maxprob;
dec.posmax(:, 1) = bins_xy(imx,1); dec.posmax(:, 2)= bins_xy(imx,2);
dec.pkmu(:, 1) = mean(prob);
dec.pkstd(:, 1) = std(prob);

if fast
    dec.decpos = dec.posmax;
else
    spreadthresh = (prctile(prob, 99));
    thresh = (prctile(prob, 99));
    maxdist = hypot(bins_xy(:,1)-dec.posmax(:, 1)', bins_xy(:,2)-dec.posmax(:, 2)');
    if calcSpread
        % weight dists by decprob
        posprob = prob; posprob(posprob<spreadthresh)=0;
        dec.spread = sqrt(sum(maxdist.*posprob)./sum(posprob));
    end
    isbad = prob<thresh & maxdist>.1;
    prob(isbad) = 0;
    sumweights = sum(prob, 'omitnan');
    
    dec.decpos(:, 1) = sum(bins_xy(:,1).*prob, 'omitnan')./sumweights;
    dec.decpos(:, 2) = sum(bins_xy(:,2).*prob, 'omitnan')./sumweights;
end

xy = dec.decpos-dec.posmax;
dec.maxoffset = hypot(xy(:, 1),xy(:, 2));
end
