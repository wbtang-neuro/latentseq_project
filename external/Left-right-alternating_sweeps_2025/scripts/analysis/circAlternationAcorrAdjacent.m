function [lags, rvals, altScore, h] = circAlternationAcorrAdjacent(angles, nLags, p)
% Calculate autocorrelogram of a circular signal for quantifying alternation
% This version has two modifications from circAlternationAcorr:
%   1. There are often gaps in the data bc of speed filtering etc. Invalid samples should be replaced with nans instead of removed.
%      Only considers adjacent samples (i.e. does not correlate angle
arguments
    angles
    nLags
    p.ax = [];
    p.vt = ~isnan(angles); % vt can be specified as a logical vector or implemented by replacing invalid angles with nans
    p.nocirc=false;
    p.t = [];
end
lags = -nLags : nLags;

a0 = angles;
v = true(size(a0));
v(1:nLags) = false;
v(end-nLags+1:end) = false;

rvals = zeros(nLags, 1);

for i = 1:numel(lags)
    lag = lags(i);
    anglesShifted = circshift(angles, lag);
    isgood = ~isnan(angles) & ~isnan(anglesShifted);
    if p.nocirc
        rvals(i) = corr(angles(v&isgood), anglesShifted(v&isgood));
    else
        rvals(i) = circ_corrcc(angles(v&isgood), anglesShifted(v&isgood));
    end
end


% Alternation score is the difference in explained variance between the
% second- and first-lagged angles
rsq1 = rvals(lags==1)^2;
rsq2 = rvals(lags==2)^2;
altScore = rsq2 - rsq1;

% PLOT
h = struct();
ax = p.ax;
if ~isempty(ax)
    h.bar = bar(ax, lags, rvals, 'faceColor', 'k');
    % h.xline = xline(ax, 0, 'r');
    xlabel(ax, "Lag (cycles)");
    ylabel("Correlation (r)");
    xlim(ax, nLags*[-1, 1]);
end

end