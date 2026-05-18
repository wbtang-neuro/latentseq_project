function [fitData, scores, P] = runShiftModelComparisonCv( ...
    spikeCounts, variablesD, interpD, fspaces, shiftData, cv, varargin)
% Run model-selection comparison for the three model variants (non-shift,
% constant-shift, theta-dependent-shift). The three models are compared
% using cross-validated log-likelihood, computed on held-out data.
%
% INPUTS
% spikeCounts: vector of 
%

inp = inputParser();
inp.addParameter("plot", false);
inp.addParameter("display", "off");
inp.addParameter("alphaPenalty", 0.01);
inp.addParameter("betaPenalty", 1e-4);
inp.addParameter("box", []);
inp.addParameter("useGpu", false);
inp.addParameter("floatClass", "single");
inp.addParameter("tol", 1e-4);
inp.addParameter("models", ["null", "const", "theta"]);

inp.parse(varargin{:});
P = inp.Results;

runOneFitLocal = @(mdlType, timeBinSelection) runOneFit( ...
    mdlType, variablesD, interpD, fspaces, shiftData, spikeCounts, ...
    "binSelection", timeBinSelection, ...
    "plot", P.plot, ...
    "display", P.display, ...
    "alphaPenalty", P.alphaPenalty, ...
    "betaPenalty", P.betaPenalty, ...
    "box", P.box, ...
    "useGpu", P.useGpu, ...
    "floatClass", P.floatClass, ...
    "tol", P.tol );

nf = cv.nFolds;

for mdlType = P.models

    if ~strcmpi(P.display, "off")
        fprintf("FITTING MODEL '%s' ...\n\n", mdlType);
    end

    F = struct();
    clear tdatCv
    for i = 1:nf
        if ~strcmpi(P.display, "off")
            fprintf("CROSSVALIDATION FOLD #%u ...\n\n'%s', ", i);
        end
        [F.cv(i), tdatCv(i)] = runOneFitLocal(mdlType, cv.train(i));
    end

    % Concatenate the LL across testing folds
    if nf
        cvL= cat(1, tdatCv.LHeldOut);
        F.cvLSum = sum(cvL);
    else
        F.cvLSum = [];
    end

    % Get whole-session fit
    if ~strcmpi(P.display, "off")
        fprintf("FITTING ON ALL DATA ...\n\n");
    end
    F.all = runOneFitLocal(mdlType, true(size(spikeCounts)) );
    fitData.(mdlType) = F;
end

% For each model type, calculate the total LLH across all held-out 
% crossvalidation folds.
L = structfun(@(s) s.cvLSum, fitData, "uni", 0);

n = numel(P.models);
scores = struct();
for m1 = 1:n
    mdl1 = P.models(m1);
    for m2 = (m1+1):n
        mdl2 = P.models(m2);
        scores.(mdl2+"_vs_"+mdl1) = L.(mdl1) ./ L.(mdl2);
    end
end

end

function [fdat, tdat] = runOneFit(mdlName, varargin)

[fdat, ~, tdat] = fitPosShiftGlm(varargin{:}, ...
    "alphaTerms",  getAlphaTerms(mdlName) );

fdat = rmfield(fdat, ["opts", "fitPosShiftGlmParams", "output"]);

end

function alphaTerms = getAlphaTerms(modelName)
alphaTerms = {};
if any(modelName==["const", "theta"])
    alphaTerms{end+1} = 'intercept';
end
if modelName=="theta"
    alphaTerms{end+1} = 'theta';
end
end