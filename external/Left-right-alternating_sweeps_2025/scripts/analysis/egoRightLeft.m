function [isRight,isLeft, prevRight, prevLeft] = egoRightLeft(egodir)
%EGORIGHTLEGT Summary of this function goes here
%   Detailed explanation goes here
egodir = egodir(:);
egodirdiff = circshift([circ_diff(egodir); nan],1);
isRight = egodirdiff<0;
isLeft = egodirdiff>0;
prevRight = circshift(isRight, 1);
prevLeft = circshift(isLeft, 1);
end

