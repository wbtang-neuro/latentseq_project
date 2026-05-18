function [fitData, iterData, tData] = ...
    fitPosShiftGlm(D, interpD, fSpaces, shiftData, spikeCounts, varargin)

inp = inputParser();
inp.addParameter('plot', false);
inp.addParameter('display', 'iter');
inp.addParameter('useParallel', false);
inp.addParameter('binSelection', true(size(spikeCounts))); % use this for crossvalidation
inp.addParameter('alphaTerms', {'intercept', 'theta_mua'});
inp.addParameter('betaTerms', ['intercept'; fieldnames(D)]);
inp.addParameter('tauTerms', {});
inp.addParameter('box', []);
inp.addParameter('limitToBox', false);
inp.addParameter('fminuncOptions', struct());
inp.addParameter('checkLlh', false);
inp.addParameter('useGpu', true);
inp.addParameter('floatClass', 'single');
inp.addParameter('alphaPenalty', 1);
inp.addParameter("betaPenalty", 1);
inp.addParameter('tol', 1e-4);
inp.parse(varargin{:});
P = inp.Results;

assert(islogical(P.binSelection));

nt = numel(spikeCounts);
[XB, prmGroups.beta] = buildXB(P, D, nt);
[XA, prmGroups.alpha, prmScales.alpha] = parseShiftTerm(D, P.alphaTerms, prmGroups, nt, P.floatClass);
[opts, params0, iterData] = getOptimOptions(P, spikeCounts, prmGroups, fSpaces);

[~, prmInds.alpha] = extractParams(params0, prmGroups.alpha);
[~, prmInds.beta] = extractParams(params0, prmGroups.beta);
prmInds.betaPos = prmGroups.beta.pos;

if P.useGpu
    interpD.Y = gpuArray(cast(interpD.Y, P.floatClass)); % should already be gpuArray
    interpD.XGrid = cellfun(@(x) {gpuArray(cast(x, P.floatClass))}, interpD.XGrid);
end


    function fcn = getCostFcn(v)

        % Apply bin selection
        spikeCountsv = spikeCounts(v);
        XBv = XB(v, :); % GLM design matrix
        XAv = XA(v, :); % Shift covariates (alpha)
        shiftDatav = structfun(@(x) x(v, :), shiftData, "uni", 0);

        % Cast data
        spikeCountsv = cast(spikeCountsv, P.floatClass);
        if P.useGpu
            XBv = gpuArray(XBv);
            XAv = gpuArray(XAv);
            spikeCountsv = gpuArray(spikeCountsv);
            shiftDatav = structfun(@(x) gpuArray(x), shiftDatav, "uni", 0);
        end

        fcn = @(params) posShiftCostFcn( ...
            params, spikeCountsv, shiftDatav, interpD, ...
            prmInds, XAv, XBv, P.checkLlh, ...
            P.alphaPenalty, P.betaPenalty, P.floatClass);

        % Clear persistent vars in cost function
        posShiftCostFcn();
    end

costFcn = getCostFcn(P.binSelection);
[params, L, exitFlag, output] = fminunc(costFcn, params0, opts);

% Evaluate the cost function with the optimized parameters, to retrieve LLH
% values for all time bins
[~, tData.L, tData.pos] = costFcn(params);

fitData = struct( ...
    'params', params, ...
    'paramGroups', prmGroups, ...
    'paramGroups2', prmInds, ...
    'llh', L, ...
    'exitFlag', exitFlag, ...
    'output', output, ...
    'paramScales', prmScales, ...
    'opts', opts, ...
    'date', datetime('now'), ...
    'fitPosShiftGlmParams', P);

% Calculate llh on held-out data
if all(P.binSelection)
    tData.posHeldOut = zeros(2, 0);
    tData.LHeldOut = [];
else
    costFcn = getCostFcn(~P.binSelection);
    [~, tData.LHeldOut, tData.posHeldOut] = costFcn(params);
end

tData = structfun(@(x) gather(x), tData, "uni", 0);

end

function pltFcns = makePlotFcns(fSpaces, prmInds, P)

fields = fieldnames(fSpaces);
pltFcn = @(b, v, s, fd, pltData, varargin) plotOptim(b, v, s, prmInds.(fd), pltData, fd, varargin{:});
pltFcns = {};

for f = 1:numel(fields)
    clear coords
    fd = fields{f};
    if ~isfield(prmInds, fd); continue, end
    F = fSpaces.(fd);
    nDims = F.nDims;
    nDimsFcn = F.fcnNDimsIn;
    gridX = F.defaultPlotX;
    gridAx = F.defaultPlotAxesGrid;
    if isempty(gridAx), gridAx = gridX; end
    gridSz = cellfun(@numel, gridX);
    [coords{1:nDimsFcn}] = ndgrid(gridX{:});
    coords = cellfun(@(x) {x(:)}, coords);
    coords = cat(2, coords{:});
    if P.useGpu
        coords = gpuArray(coords);
    end

    fdPltArgs = {'grid', gridAx};
    if strcmp(fd, 'pos')
        fdPltArgs = [fdPltArgs, {'box', P.box}];
    end

    clear pltDataTmp
    zAll =  gather(F.evaluate(coords, 1));
    for i = 1:nDims
        z = zAll(:, i);
        if nDimsFcn > 1
            z = reshape(z, gridSz);
        end
        pltDataTmp{i} = z;
    end
    pltData = cat(nDimsFcn+1, pltDataTmp{:});
    pltFcns{end+1} = @(b, v, s) pltFcn(b, v, s, fd, pltData, fdPltArgs{:});
end

end

function stop = plotParams(b, v, s)
bar(b);
set(gca, 'xlim', [0, numel(b)+1]);
title("model parameters");
stop = false;
end

function [XB, prmInds] = buildXB(P, D, nt)
betaTerms = P.betaTerms;
XB = zeros(nt, 0, P.floatClass);

for f = 1:numel(betaTerms)
    fd = betaTerms{f};
    if strcmp(fd, 'intercept')
        XTmp = ones(nt, 1, P.floatClass);
    else
        XTmp = D.(fd).Y;
    end
    nc = size(XTmp, 2);
    prmInds.(fd) = size(XB, 2) + (1:nc)';
    XB = [XB XTmp];
end
end

function [opts, params0, iterData] = getOptimOptions(P, spikeCounts, prmInds, fSpaces)
% Intialize parameters and fminunc "TypicalX" values

nParamsTot = calcNParams(prmInds);
params0 = 0.0001 * randn(nParamsTot, 1);

meanSpikeCount = gather(mean(spikeCounts));
params0(1) = log(meanSpikeCount);

stepsz = 5 * double(eps(P.floatClass)) .^ (1/2); % why are we setting this manually?

% if isempty(opts)
opts = optimoptions( 'fminunc', ...
    'Hessian', 'off', ...
    'Display', P.display, ...
    'OptimalityTolerance', P.tol, ...
    'StepTolerance', P.tol, ...
    'FunctionTolerance', P.tol, ...
    'SpecifyObjectiveGradient', false, ...
    'Algorithm', 'quasi-newton', ...
    'UseParallel', P.useParallel, ...
    'MaxFunctionEvaluations', Inf, ...
    'OutputFcn', {@outputFcn}, ...
    'HessUpdate', 'bfgs', ...
    'FiniteDifferenceType', 'forward', ...
    'FiniteDifferenceStepSize', stepsz, ...
    'DiffMinChange', 0);
% else
%     opts.OutputFcn = {@outputFcn};
% end

opts.TypicalX = ones(nParamsTot, 1) * 0.01;
opts.TypicalX(1) = abs(params0(1));

if P.plot
    plotPrmInds = prmInds.beta;
    plotFSpaces = fSpaces;

    fds = string(fieldnames(prmInds.alpha));
    fds = fds(fds ~= "intercept");
    for f = 1:numel(fds)
        fd = fds(f);
        pfd = sprintf("alpha_%s", fds(f));
        plotPrmInds.(pfd) = prmInds.alpha.(fd);
        plotFSpaces.(pfd) = fSpaces.(fd);
    end
    opts.PlotFcn = [{@plotParams}, makePlotFcns(plotFSpaces, plotPrmInds, P)];
end

iterData = [];
    function [stop, optnew, changed] = outputFcn(params, optValues, state)
        if strcmpi(state, 'iter')
            tmp = optValues;
            tmp.params = params;
            iterData = [iterData; tmp];
        end
        stop = false;
        optnew = false;
        changed = false;
    end

end