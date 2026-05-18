function varargout = centeredGridVector(varargin)
% Generate a vector of equally spaced points centered on zero
%
% X = CENTEREDGRIDVECTOR(N) generates a vector of length N, with values
% equal to ((-N/2)-1) : ((N/2)-1)
%
% X = CENTEREDGRIDVECTOR(N, STEP) functions similarly, but also multiplies 
% output X by the scalar STEP, such that STEP is the increment between
% points in X.

np = varargin{1};
nd = numel(np);

if nargin == 1
    step = 1;
elseif nargin == 2
    step = varargin{2};
end

if isscalar(step)
    step = repmat(step, 1, nd);
end

for d = 1:nd
    x = (1:np(d))*step(d);
    x = x-mean(x);
    varargout{d} = x;
end

end