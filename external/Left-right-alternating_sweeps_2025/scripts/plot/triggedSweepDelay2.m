function [swpt, swpx, iends] = triggedSweepDelay2(D,br,multi)
%TRIGGEDSWEEPDELAY2 Summary of this function goes here
if nargin==2
    multi = 0;
end
if ~multi
maxlag = 15;

poshpf = D.dec.(br).decpos-D.dec.(br).possm;
posr = rotate2d(poshpf, -D.hd);
decpos= posr(:, 1);
% decpos = hypot(posr(:, 1), posr(:, 2));

[trigavg] =plotTriggeredAverage(decpos, D.chk.iStart, maxlag);

[pks, loc] = findpeaks(trigavg);
lags = -maxlag:maxlag;
loct = lags(loc);
loct(loct<0) = nan;
[~, iloc] = min(loct);
iend = loc(iloc);

[pks, loc] = findpeaks(-trigavg);

loct = lags(loc);
[~, iloc] = min(abs(loct));
istart = loc(iloc);

swpt = 1000*D.dt*lags(istart:iend)';
swpx = trigavg(istart:iend);
iends = numel(swpx);
else
maxlag = 40;

poshpf = D.dec.(br).decpos-D.dec.(br).possm;
posr = rotate2d(poshpf, -D.hd);
decpos= posr(:, 1);
% decpos = hypot(posr(:, 1), posr(:, 2));

[trigavg] =plotTriggeredAverage(decpos, D.chk.iStart, maxlag);
lags = -maxlag:maxlag;
[pks, iends] = findpeaks(trigavg);
[pks, istarts] = findpeaks(-trigavg);

% Make sure we start with a trough
iends(iends<istarts(1)) = [];
istarts(istarts>iends(end)) = [];

swpt = nan(size(lags));
swpx = nan(size(lags));
for i = 1:numel(iends)
    swpt(istarts(i):iends(i)) = 1000*D.dt*lags(istarts(i):iends(i))';
    swpx(istarts(i):iends(i)) = trigavg(istarts(i):iends(i));
end

end

