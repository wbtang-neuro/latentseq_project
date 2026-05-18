function s = alternationScore4(x)
% Quantifies alternation between adjacent values in a circular timeseries
% signal. The score is calculated at each time point in the signal, using a 
% 3-sample sliding window. The Score ranges from 0 (a straight line) to 1 
% (a perfect zigzag). Random, uniformly distributed angles will yield a 
% score of around 0.4034. The score is independent of the magnitude of the 
% change between samples. If consecutive values are exactly the same, the 
% resultant score is NaN.

x = x(:);

% The input vector x may be of any length
iblk = [-1, 0, 1];
s = nan(size(x));

for b = 2 : (numel(x)-1)
    inds = b + iblk;
    xb = x(inds);
    s(b) = calcScore(xb);
end
end

function score = calcScore(p)
a = circ_dist(p(2), p(1));
b = circ_dist(p(3), p(2));
d = circ_dist(a, b);
score = abs(d) / (2 * max(abs(a), abs(b)));
end