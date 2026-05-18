function lmtDat = createLmtDataStruct(D, Us, mdl)

if mdl.enableXAlignment()
    mdl.calcXTform();
end
mdl.alignX();

nd = mdl.ndims;
% Interpolation grid of "fixed" tuning coordinates
gv0 = {D.gv.(mdl.name)};
gv0 = repmat(gv0, 1, nd);
ggi0 = cell(1, nd);
[ggi0{1:nd}] = ndgrid(gv0{:});
fszi = size(ggi0{1});
ggi0 = cellfun(@(x) {x(:)}, ggi0);
ggi0 = [ggi0{:}];
ggiT = mdl.alignX(ggi0, "reverse"); % Xinit -> X0

% N.B. ggi dims will be ordered as [X, (Y)], i.e. as a
% meshgrid rather than NDgrid

fsz = repmat(mdl.nf, 1, nd);
if isscalar(fsz), fsz(2) = 1; end

% For each unit in LMT model, find its location in the
% complete array of units
v = ismember(mdl.unitIds(:), [Us.id]');
assert(all(v));

lmtDat = struct( ...
    "XA",            mdl.XAligned, ...
    "F",             mdl.logTuningCurvesF, ...
    "gv",            {mdl.Fgridv}, ...
    "ggi",           ggiT, ...
    "unitIds",       mdl.unitIds(:), ...
    "gsz",           fsz, ...
    "gszi",          fszi, ...
    "hparams",       mdl.hparams);
lmtDat = arrayfunRecursive(lmtDat, @gather, 'gpuArray');

end