function [res, P] = xcorr_nodeconv2(t1, t2, P)
% xcorr
arguments
    t1
    t2
    P.binSize = 2/30e3 * 10; % should be a whole number of Npx samples (avoid aliasing)
    P.nBins = 101; %change to maxlag in ms
    P.kernelType = "gaussian"
    P.kernelWidth = 7; %change to ms
    P.lag = [];
    %     P.normMethod = "mean"
    % opts.peakTimeRange = [-10, 10]*1e-3; % time range for detecting peak
    P.peakTimeRange = []; % time range for detecting peak
%     P.ignoreExact = false; % Don't ignore center bin for xcorrs by default
end

[xc.raw, tlags] = xcpp(t1, t2, P.binSize, P.nBins, 0);
[p, xc.pred, ~] = cch_conv(xc.raw, P.kernelWidth, P.kernelType);
xc.excess = xc.raw-xc.pred;
xc.excessM = xc.excess / mean(xc.raw);
xc.pvals = single(p);
% Check if we have a peak above the significance threshold
cpk = xc.excess;

[~, imx] = max(abs(cpk));
notpk = true(size(cpk));
notpk(imx) = false;

res.std = std(cpk(notpk));
xc.excessZ = xc.excess./res.std;
res.pkvalM = xc.excessM(imx);
res.pkvalZ = xc.excessZ(imx);
res.ipeak = imx;

res.totalCnt = sum(xc.raw);
res.tlags = tlags;
res.xc = xc;
res.nspk_i = numel(t1);
res.nspk_j = numel(t2);
end

