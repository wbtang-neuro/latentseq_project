function [percentflick, pval] = computeAlternationPercent(sd, nodiff)
%COMPUTEALTERNATIONPERCENT Summary of this function goes here
%   Detailed explanation goes here
if nargin==1
    sddiff = circ_diff(sd);
else
    sddiff = sd;
end
diffsign = sign(sddiff);
diffdiff = diff(diffsign);

isalternation = abs(diffdiff)==2;
npossible = sum(~isnan(diffdiff));
percentflick = sum(isalternation)./npossible;
pval = nan;
end

