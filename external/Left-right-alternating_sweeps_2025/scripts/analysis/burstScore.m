function [scoreDiff, scoreRatio, acorr] = burstScore(spikeTimes)
binSize = 1e-3;
nBins = 101;
normalization = 0;

[c, b] = acpp(spikeTimes, binSize, nBins, normalization);
c(51)=0;
[scoreDiff, scoreRatio] = burstScoreAcorr(b, c);

if nargout == 3
    acorr = struct("lags", b, "counts", c);
end

end