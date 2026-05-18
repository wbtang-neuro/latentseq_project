function [trigAvg,h] = plotTriggeredAverage(data,trigIdx,maxlag, method, ax, circ)
%PLOTTRIGGEREDAVERAGE Summary of this function goes here
%   Detailed explanation goes here

lags = -maxlag:maxlag;
trigAvg = zeros(numel(lags), 1);
trigStd = zeros(numel(lags), 1);

for l = 1:numel(lags)
    idx = trigIdx+lags(l);
    idx = min(idx, numel(data));
    idx = max(idx, 1);
    dat = data(idx);
    trigAvg(l) = mean(dat, 'omitnan');
    if nargin>3 && method=="median"
        trigAvg(l) = median(dat, 'omitnan');
    end
    trigStd(l) = std(dat, 'omitnan');
    if nargin ==5
        trigAvg(l) = circ_mean(dat);
    end
end

if nargin>4
    h = plot(lags, trigAvg, 'k');
    xline(0, 'r');
end
end

