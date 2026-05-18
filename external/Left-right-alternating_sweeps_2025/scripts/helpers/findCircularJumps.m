function [ijump, hasJump] = findCircularJumps(phi, thresh)
% Find discontinuities in N-dimensional circular/toroidal data.

% Assumes that circular axes are in columns. If phi is a matrix, it will be
% interpreted as a torus, and a jump on any axis will be counted.

if nargin < 2 || isempty(thresh), thresh = pi; end

dp = diff(phi);                % Incremental phase variations
hasJump = abs(dp) >= thresh;
ijump = find(any(hasJump, 2)) + 1;

end