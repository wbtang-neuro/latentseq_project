function [res] = runPvIdDecoding(p)
%RUNLINEARDECPOSRECS Runs corr-based decoder on a set of recordings and
%saves results in specified folder
arguments
    p.recs = [];
    p.fld = [];
    p.process = 1;
    p.load = 0;
    p.save = 0;
end

S = sweepsSettings;
if isempty(p.fld)
    fld = fullfile(S.dataRoot, "code", "results", "pv_dec", "mec", "id");
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
datafld = fullfile(S.dataRoot, "navigation", "of");
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
        units = D.units.mec;
    
        %% Bin spikes
        disp('Counting spikes...')
        % for u = 1:numel(units)
        %     units(u).sc = gather(binSpikes(units(u).spikeTimes, [], D.t, useGpu));
        % end
        
        sc = reconstructSpikeCounts({units.spikeInds}, D.nt);
    %     %%

        dectmpl = getDectmpl(D);
        % sc = [units.sc];
        active = sc; active(sc>0) = 1;
        dectmpl.nactive = sum(active, 2);
        dectmpl.nspikes = sum(sc, 2);
    
        % Get baseline tuning
        tic
        clear tuning
        nu = numel(units);
        for u = 1:nu
           tc = units(u).rmf.hd;
           tuning(:, u) = tc.z(:);
        end
        %%
        tic
        prob = decodePv(tuning=tuning, spikeCounts=sc, smoothSpikes=S.tsm.id);
        toc
           %%
        dec = dectmpl;
        [maxprob, imx] = max(prob, [], 2);
        dec.imx = imx;
        dec.idmuAll = circ_mean(dec.gv, prob')'; 
        dec.idmu = dec.idmuAll(chk.iCen);
        chk.id = dec.idmu;

        %% n cells and n dir cells
        chk.nu = numel(units);
        chk.ndir = sum(ismember([units.cellType], ["id", "conjunctive"]));
        chk.nconj = sum([units.cellType]=="conjunctive");
        chk.ngrid = sum([units.isGrid]);
        %% store
        res(r) = chk;
        
        if p.save
            disp("Saving..")
            if ~isfolder(fld), mkdir(fld); end
            subfld = fullfile(fld, rec);
            save(subfld, "-struct", "chk", "-v7", "-nocompression");
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
function dectmpl = getDectmpl(D)
S = sweepsSettings;
vt = true(size(D.t));

gve = linspace(0, 2*pi, 101)';
dec.gve = gve;
dec.gv = edg2cen(dec.gve);
dec.gv = S.gv.angular;
dec.gve = S.gve.angular;
dec.prob = [];
dec.t = D.t(vt);
dec.x = D.x(vt);
dec.y = D.y(vt);
dec.id = D.hd(vt);
if isfield(D, "id") && ~isempty(D.id)
    dec.id = D.id(vt);
end
dec.theta = D.theta;
dec.lmtpos = D.lmt.mec.pos.XA;
dec.speed = D.speed;
dec.hd = D.hd;          
dec.idprob = ones(size(dec.id(vt)))*.1;
dectmpl = dec;
end
