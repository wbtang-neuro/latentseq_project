function [XSh, prmInds, prmScales] = parseShiftTerm(D, shiftTerms, prmIndsAll, nt, cls)

if nargin < 5 || isempty(cls), cls = 'double'; end

shiftTermParts = cellfun(@(s) {strsplit(s, '*')}, shiftTerms);
nTerms = numel(shiftTerms);
idxPrm = calcNParams(prmIndsAll);

XSh = zeros(nt, 0, cls);
prmInds = struct();
prmScales = struct();

for t = 1:nTerms
    parts = shiftTermParts{t};
    nParts = numel(parts);
    if nParts == 1 && strcmpi(parts{1}, 'intercept')
        X = ones(nt, 1);
        subScales = 1;
    else
        clear cSub* combs
        for p = 1:nParts
            subterm = parts{p};
            cSubN(p) = size(D.(subterm).Y, 2);
            cSubGrid{p} = (1:cSubN(p))';
        end
        [combs{1:nParts}] = ndgrid(cSubGrid{:});
        combs = cellfun(@(x) {x(:)}, combs);
        combs = [combs{:}];
        X = ones(nt, size(combs, 1), cls);
        for p = 1:nParts
            subterm = parts{p};
            X = X .* D.(subterm).Y(:, combs(:, p));
        end
        subScales = std(X);
    end
    
    np = size(X, 2);
    cInds = 1:np;
    fd = strjoin(parts, '_');
    X = X ./ subScales;
    XSh = [XSh X];
    prmInds.(fd) = idxPrm + cInds;
    prmScales.(fd) = subScales;
    idxPrm = idxPrm + np;
end

end