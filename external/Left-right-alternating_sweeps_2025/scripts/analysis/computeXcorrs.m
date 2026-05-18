function [] = computeXcorrs(recs, S, force)
%COMPUTEXCORRS Summary of this function goes here
%   Detailed explanation goes here
% sweepsSetup;
S = SweepsSettings;
% For monosynaptic connections
baseDirname = fullfile(S.dataRoot, "results", "xcorrs");
if ~isfolder(baseDirname)
    mkdir(baseDirname)
end
datafld = fullfile(S.dataRoot, "navigation", "of");
connLagRange = [.1, 8]*1e-3;
connLagRange = [.3, 5.5]*1e-3;

if nargin==2
    force = 1;
end

nrecs = size(recs, 1);
kernelWidth = 5;
binSize = 1e-3;

% Calculate raw xcorrs for all pairwise combs of units
for r = 1:nrecs
    rec = recs(r, :);
    fprintf("Rec '%s' ...\n ", rec);
    fnOut = fullfile(baseDirname, rec + ".mat");
    if isfile(fnOut) && ~force
        fprintf("Already exists. Skipping\n");
        continue;
    end
    fname = fullfile(datafld, rec+".mat");
    tmp = load(fname);
    D = tmp.Dsession;
    Us = D.units.mec;
    dat = struct();
    [dat.conns, dat.tlags, dat.xcorrs] = ...
        calc_connection_metrics(Us, [], ...
        kernelWidth=kernelWidth, ...
        connTimeRange=connLagRange,...
        binSize=binSize, ...
        useParfor=true);
    dat.unitIds = [Us.id]';
    dat.kernelWidth = kernelWidth;
    dat.binSize = binSize;
    dat.connLagRange = connLagRange;
%
    % Add info abaout rec
    info.recDuration = numel(D.t)*median(diff(D.t));
    dat.recinfo = info;
    fprintf("Rec duration: %.2fmin\n", info.recDuration./60)
    save(fnOut, "-struct", "dat", "-v7.3", "-nocompression");

end
end

