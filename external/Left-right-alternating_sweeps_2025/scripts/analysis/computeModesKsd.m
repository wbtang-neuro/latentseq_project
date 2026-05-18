function [modes,lefthist,righthist] = computeModesKsd(egosd)
%COMPUTEMODESKSD Summary of this function goes here
%   Detailed explanation goes here
egosd = egosd(:);
egosddiff = circshift([circ_diff(egosd); nan],1);
egosign = sign(circshift(egosddiff, 1));

gv = linspace(-pi,pi, 100);
righthist = ksdensity(egosd(egosign>0), gv);
lefthist = ksdensity(egosd(egosign<0), gv);
[~, imx(1)]=max(righthist);
[~, imx(2)]=max(lefthist);
modes = gv(imx);
end


