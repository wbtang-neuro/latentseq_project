function [conns, tlags, xcorrs] = calc_connection_metrics(Us1, Us2, P)
%CALC_CONNECTION_METRICS Summary of this function goes here
%{
connection yes/no (excitatory, inhibitory)
zscore peak height and lag
pairwise distance on shank
waveform similarity
%}

arguments
    Us1
    Us2
    P.binSize = 1/30e3 * 30
    P.nBins = 101
    P.connTimeRange = [.3, 5.5]*1e-3
    P.kernelWidth = 21
    P.useParfor = true
    P.alpha = .001;
    P.deconv = false;
    P.lag = -0.2*1e-3;
end

% conns = struct;
samePairing = nargin < 2 || isempty(Us2) || isequal(Us1, Us2);
if samePairing
    Us2 = Us1;
end

nu1 = numel(Us1);
nu2 = numel(Us2);

uDists = getDists(Us1, Us2);

collectXcs = false;
collectRes = false;
if nargout > 2
    c0 = zeros(P.nBins, nu1, nu2);
    xcR = uint32(c0);
    xcD = single(c0);
    xcP = single(c0);
    xcE = single(c0);
    xcUb = single(c0);
    xcPval = single(c0);
    collectXcs = true;
    if nargout == 5
        collectRes = true;
    end
end

[~,~,res] = runOnePair(Us1(1).spikeTimes, Us1(1).spikeTimes, P, samePairing);
tlags = res.tlags;
if P.useParfor
    parfor i = 1:nu1
        fprintf("Iter %u of %u\n", i, nu1);
        t1 = Us1(i).spikeTimes;
        for j = 1:nu2

            t2 = Us2(j).spikeTimes;
            samePairing = i==j;
            [C, xc, res] = runOnePair(t1, t2, P, samePairing);
            C.uDist  = uDists(i,j);
            conns0(i, j) = C;
            if collectXcs
                xcR(:, i, j) = xc.raw;
                xcP(:, i, j) = single(xc.pred);
                xcE(:, i, j) = single(xc.excess);
                xcPval(:, i, j) = single(xc.pvals);
            end

            if collectRes, resAll(i, j) = res; end
        end
    end

else
    for i = 1:nu1
        fprintf("Iter %u of %u\n", i, nu1);
        t1 = Us1(i).spikeTimes;
        for j = 1:nu2

            t2 = Us2(j).spikeTimes;
            samePairing = i==j;
            samePairing = false;
            [C, xc, res] = runOnePair(t1, t2, P, samePairing);
%             C.uDist  = uDists(i,j);
            conns0(i, j) = C;
            if collectXcs
                xcR(:, i, j) = xc.raw;
                xcP(:, i, j) = single(xc.pred);
                xcE(:, i, j) = single(xc.excess);
                xcPval(:, i, j) = single(xc.pvals);
            end

            if collectRes, resAll(i, j) = res; end
        end
    end
end

% Output as a scalar struct
conns = struct();
intFields = ["nspk_i", "nspk_j", "totalCnt"];
for fd = fieldnamesstr(conns0)
    val = reshape([conns0.(fd)], size(conns0));
    if any(fd == intFields)
        if all(val(:)<=intmax("uint32"))
            val = uint16(val);
        else
            val = uint32(val);
        end
    else
        val = single(val);
    end
    conns.(fd) = val;
end

conns.nspk_i = conns.nspk_i(:, 1);
conns.nspk_j = conns.nspk_j(1, :);

if collectXcs
    xcorrs = struct( ...
        "raw", xcR, ...
        "pred", xcP, ...
        "excess", xcE,...
        "pvals", xcPval);
end
end

function [C, xc, res] = runOnePair(t1, t2, P, samePairing)

C = struct();

[res, tlags] = xcorr_nodeconv2(t1, t2+P.lag, ...
nBins=P.nBins, ...
binSize=P.binSize, ...
kernelWidth=P.kernelWidth);
res.tlags = res.tlags-P.lag;
xc = res.xc;
% elseif
ipeak = res.ipeak;

tpeak = res.tlags(ipeak);

% Compute CCG noise
acausalrng = restrictq(res.tlags, [-5,.9]*1e-3);
C.maxacausal = max(xc.excess(acausalrng));

notpeak = true(P.nBins, 1); 
if tpeak>=P.connTimeRange(1) && tpeak<=P.connTimeRange(2)
    notpeak((-1:1)+ipeak) = false;
else
    notpeak(ipeak) = false;
end
C.maxOther = max(xc.excess(notpeak));
C.stdAll = std(xc.excess);
C.std = res.std;
isHigh = xc.excess>2*C.stdAll | xc.pvals<.01 | xc.excess>.5*xc.excess(ipeak);

if tpeak>=P.connTimeRange(1) && tpeak<=P.connTimeRange(2)
    C.pkvalM = res.pkvalM;
    C.pkvalZ = res.pkvalZ;
    C.pkExcess = xc.excess(ipeak);
    C.pkPval = xc.pvals(ipeak);
    C.ipeak = ipeak;
    C.tpeak = tpeak;
else
    C.pkvalM = nan;
    C.pkvalZ = nan;
    C.pkExcess = nan;
    C.ipeak = nan;
    C.tpeak = nan;
    C.pkPval = nan;
end


if ~samePairing && ~isnan(C.pkvalM) && isHigh(ipeak)
    % Get peak width, start and end
%     size(isHigh)
    lowInds = find(~isHigh)';
    C.pkStart = min(max(lowInds(lowInds<ipeak)+1), ipeak);
    C.pkEnd = max(min(lowInds(lowInds>ipeak)-1), ipeak);
    C.pkWidth = C.pkEnd - C.pkStart + 1;
    C.conn = true;
    C.zerolag = isHigh(52); %assuming 51 is centerbin
    C.acausal = any(isHigh(acausalrng));
    C.connectivity = sum(xc.excess(C.pkStart:C.pkEnd));
else
    C.pkStart = nan;
    C.pkEnd = nan;
    C.pkWidth = nan;
    C.conn = false;
    C.connectivity = nan;
    C.zerolag =nan;
    C.acausal=nan;
end

C.nspk_i = res.nspk_i;
C.nspk_j = res.nspk_j;
C.totalCnt = sum(xc.raw);

end

% functions
function uDists = getDists(units1, units2)
p1 = cat(1, units1.shankPos) + 99999 * cat(1, units1.probeId);
p2 = cat(1, units2.shankPos) + 99999 * cat(1, units2.probeId);
uDists = pdist2(p1, p2);
uDists(uDists>10000) = nan;
end



