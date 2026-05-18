function [scoreDiff, scoreRatio] = burstScoreAcorr(b, c)
% Calculate burst score from a temporal autocorrelogram

assert(all(c>=0), "Input autocorrelogram must be nonnegative"); % make sure it's a histogram (normalized is OK, though).

c = c/mean(c);

bsz = b(2)-b(1);
assert(bsz<=0.002, "bin size must not be larger then 2 ms");

rngB = [0.002, 0.010];
rngNb = [0.013, 0.050];

vB = b>=rngB(1) & b<=rngB(2);
vNb = b>=rngNb(1) & b<=rngNb(2);

scoreB = mean(c(vB));
scoreN = mean(c(vNb));

scoreDiff = scoreB - scoreN;
scoreRatio = scoreB ./ scoreN;

end