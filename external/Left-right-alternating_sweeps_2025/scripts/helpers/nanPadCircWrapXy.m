function [xpad, ypad] = nanPadCircWrapXy(x, y)
% insert NaNs in X and Y vectors at locations where Y wraps circularly
[ypad, ~, indsOut] = nanPadCircWrap(y);
xpad = nan(size(ypad));
xpad(indsOut) = x;
end