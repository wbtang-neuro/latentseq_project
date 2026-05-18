function [v, xv] = restrictq(x, xranges)
% RESTRICTQ simplified version of RESTRICT using default behavior.

nranges = size(xranges, 1);
v = false(size(x));
x0 = x;

for r = 1:nranges
    rng = xranges(r, :);
    x = x0;
    vtmp = x >= rng(1) & x <= rng(2);
    v = v | vtmp;
end

if nargout==2
    xv = x(v);
end

end