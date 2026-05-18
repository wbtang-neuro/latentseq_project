function [muoffset, rho, pval, muabsoffset, res] = computeDirAlignment(dir1,dir2, norot)
%COMPUTEDIRALIGNMENT Summary of this function goes here
%   Detailed explanation goes here
invalid = (isnan(dir1)|isnan(dir2));
dir1(invalid)=[]; dir2(invalid)=[];
dir1 = double(dir1);
dir2=double(dir2);
if nargin>2 && norot
    [rho, pval] = circ_corrcc_no_rotation(dir1, dir2);
else
    [rho, pval] = circ_corrcc(dir1, dir2);
end
offset = circ_dist(dir1, dir2);
muoffset = circ_mean(offset);
muabsoffset = median(abs(offset));
res.mvl = circ_r(offset);
res.rayleigh = circ_rtest(offset);
res.std = circ_std(offset);
res.sem = res.std./numel(dir1);
end

