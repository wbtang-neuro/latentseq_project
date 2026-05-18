function [h, validCycles] = plotShadedThetaCycles(ax, cycleTimes, timeRange)
if nargin < 3 || isempty(timeRange), timeRange = ax.XLim; end

t1 = cycleTimes(:, 1);
if size(cycleTimes, 2) == 1
    t2 = [t1(2:end); Inf];
else
    t2 = cycleTimes(:, 2);
end

v2 = restrictq(t2, timeRange);
v1 = restrictq(t1, timeRange);

truncstart = find(~v1 & v2);
truncend = find(v1 & ~v2);

if ~isempty(truncstart)
    t1(truncstart(1)) = timeRange(1);
end
if ~isempty(truncend)
    t1(truncend(1)+1) = timeRange(2);
end

validCycles = restrictq(t1, timeRange);
tnewcyc = t1(validCycles);

hold(ax, "on");

ncyc = numel(tnewcyc);
yl = ax.YLim;

if ncyc < 2
    h = [];
    return
end

for i = 1:ncyc-1
    x = tnewcyc([i, i, i+1, i+1, i]);
    y = yl([1, 2, 2, 1, 1])';
    if rem(i, 2) == 1
        alpha = 0.1;
    else
        alpha = 0;
    end
    h(i) = patch(ax, x, y, "k", 'faceAlpha', alpha, 'edgeColor', [0, 0, 0]+0.6);
end


end