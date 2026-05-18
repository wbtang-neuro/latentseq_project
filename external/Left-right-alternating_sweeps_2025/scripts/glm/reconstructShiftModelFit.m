function [tuning, ggs, gvs] = reconstructShiftModelFit(fits, fspaces, sessionType, paramGroups)

if nargin < 4, paramGroups = []; end

warning("off", "MATLAB:structOnObject");
S = struct(SweepsSettings()); % some bugginess with parfor "threads"
warning("on", "MATLAB:structOnObject");

% Get binning grid for all covariates
gvNames = struct("hd", "angular", "id", "angular", "theta", "angular");
[gvNames.pos, gvNames.pos_coarse] = getPosGvField(sessionType);

fitTypes = ["alpha", "beta"];


nfits = numel(fits);
tuning = [];
gvs = [];
ggs = [];

for t = 1:numel(fitTypes)

    fitType = fitTypes(t);
    if isempty(paramGroups)
        pInds = fits(1).paramGroups.(fitType);
    else
        pInds = paramGroups.(fitType);
    end
    
    pGrpNames = string(fieldnames(pInds));
    hasIntercept = any(pGrpNames=="intercept");
    if hasIntercept
        pGrpNames(pGrpNames == "intercept") = [];
        for f = 1:nfits
            tuning(f, 1).(fitType).intercept = fits(f).params(pInds.intercept);
        end
    end

    for g = 1:numel(pGrpNames)

        % Extract parameters for this variable's contribution to the
        % fit
        varname = pGrpNames(g);
        iprm = pInds.(varname);
        fspace = fspaces.(varname);

        % find the bin grid for reconstructing this variable
        gvname = gvNames.(varname);
        gv = S.gv.(gvname);
        gg = S.gg.(gvname);
        gg1 = S.gg1.(gvname);
        sz = size(gg{1});

        ggs.(varname) = gg;
        gvs.(varname) = gv;

        z1 = fspace.evaluate(gg1, ones(1, numel(iprm)));

        if nfits>1
            % Now iterate through fits for different cells
            parfor f = 1:nfits
                tuning(f).(fitType).(varname) = calcOneTuning(fits(f).params, iprm, z1, sz);
%                 p = fits(f).params(iprm);
%                 z = sum(z1 .* p', 2);
%                 z = reshape(z, sz);
%                 tuning(f).(fitType).(varname) = z;
            end
        else
            tuning(1).(fitType).(varname) = calcOneTuning(fits(1).params, iprm, z1, sz);
        end
    end
end

end

function z = calcOneTuning(params, iprm, z1, sz)
p = params(iprm);
z = sum(z1 .* p', 2);
z = reshape(z, sz);
end