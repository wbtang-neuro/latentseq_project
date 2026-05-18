function [legal] = isLegal(v)
vv = v(:);
legal = isreal(vv) & ~any(isnan(vv)) & ~any(isinf(vv));
% legal = sum(any(imag(vv)))==0 & sum(isnan(vv))==0 & sum(isinf(vv))==0;