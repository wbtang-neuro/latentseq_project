function [coalternation, pval] = computeCoAlternation(sd, usediff)
%COMPUTECOALTERNATION Summary of this function goes here
%   Detailed explanation goes here
sddiff = sd;
if nargin>1 && usediff
    sddiff(:, 1) = [nan;circ_diff(sd(:, 1))];
    sddiff(:, 2) = [nan;circ_diff(sd(:, 2))];
end
diffsign = sign(sddiff);
invalid = isnan(sum(diffsign, 2));
diffsign(invalid, :)=[];
nsameside = sum(diffsign(:, 1) == diffsign(:, 2));
npossible = sum(~invalid);
coalternation = nsameside./npossible;
pval=myBinomTest(nsameside,npossible,.5,'one');
end

