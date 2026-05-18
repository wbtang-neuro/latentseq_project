function [Us, inds, wasFound] = findUnits(Us, ids, allowUnmatched)
% Returns units in Us indexed by ids 
if nargin < 3 || isempty(allowUnmatched), allowUnmatched = false; end
if isempty(ids)
    wasFound = [];
    inds = [];
else
    [wasFound, inds] = ismember(ids, [Us.id]);
end
if allowUnmatched
    inds = inds(wasFound);
else
    assert(all(wasFound), "Failed to find all requested units");
end
Us = Us(inds);
end