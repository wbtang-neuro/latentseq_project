function [xout, indsIn, indsOut] = insertValues(x, inds, insertVals)
% Insert values into an array at specified indices
if nargin < 3 || isempty(insertVals), insertVals = nan; end

x0 = x;

x = gather(x);
inds = gather(inds);

% n = numel(x);
ncolsx = size(x, 2);
nrowsx = size(x, 1);
inds = inds(:); % gpu issorted() only accepts vector input
ninsert = numel(inds);
nRowsout = nrowsx+ninsert;

xout = zeros(nRowsout, ncolsx, "like", x);
indsIn = zeros(nRowsout, 1);
indsOut = zeros(nrowsx, 1);

if isscalar(insertVals)
    insertVals = repmat(insertVals, ninsert, 1);
end

if ~issorted(inds)
    [inds, isort] = sort(inds);
    insertVals = insertVals(isort);
end

c = 0;

for i = 1:nrowsx
    iout = i + c;
    % If padding is to be inserted at current position
    while c < ninsert && i==inds(c+1)
        c = c+1;
        xout(iout, :) = insertVals(c);
        iout = i + c;
    end
    indsIn(iout) = i;
    indsOut(i) = iout;
    xout(iout, :) = x(i, :);
end

xout = cast(xout, "like", x0);

end