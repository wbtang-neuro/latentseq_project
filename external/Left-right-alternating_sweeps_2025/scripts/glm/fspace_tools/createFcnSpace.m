function S = createFcnSpace(fcnType, fcnParams, varargin)
%CREATEFCNSPACE convenience function for creating function spaces
%
% S = CREATEFCNSPACE(FTYPE, PARAMS) creates FunctionSpace object S using
% function specified by string FTYPE. FTYPE may be one of the following
% options:
% 'vm'      - Von-Mises PDF
% 'guas'    - 1-D Gaussian PDF
% 'gaus2'   - 2-D Gaussian PDF
%
% PARAMS is a 1-D cell array containing the parameters for each of the
% functions. Each element should contain the values for one function 
% parameter, specified as an F*N array where F is the number of function
% basis vectors and N is the dimensionality of the parameter. The required
% parameters for the different function types are:
%
% vm:
%       {1} kappa       (concentration)
%       {2} thetaHat    (center angle)
%
% gaus/gaus2
%       {1} mu          (center angle)
%       {2} sigma       (concentration)

inp = inputParser();
inp.addParameter('plottingCoords', {});
inp.addParameter('plottingLims', []);
inp.addParameter('plottingNGrid', 101);
inp.parse(varargin{:});
P = inp.Results;

% Define the common function and its parameters
[fcnCommon, nDimsIn, nDimsOut, prmNames, fsName, plotFcn] = basisFcnPreset(fcnType);

nVals = cellfun(@(x) size(x, 1), fcnParams);
nFcns = max(nVals);
if nFcns > 1
    inds = find(nVals==1);
    fcnParams(inds) = cellfun(@(x) {repmat(x, nFcns, 1)}, fcnParams(inds));
end

plottingCoords = parsePlottingCoords(fcnType, nDimsIn, fcnParams, P);

for d = 1:nFcns
    prmValsTmp = cellfun(@(x) {x(d, :)}, fcnParams);
    prmValsStr = prmValsTmp;
    for a = 1:numel(prmValsStr)
        val = prmValsStr{a};
        if isnumeric(val)
            if ~isempty(val)
                str = sprintf('%.3f, ', val);
                prmValsStr{a} = ['[' str(1:end-2)  ']' ];
            end
        end
    end
    args = [prmNames(:)'; prmValsStr(:)'];
    prmStr = sprintf('%s = %s, ', args{:});
    prmStr(end-1:end) = '';
    name = sprintf('%s (%s)', fsName, prmStr);
    fcn = @(X) fcnCommon(X, prmValsTmp{:});
    bf = BasisFunction(fcn, nDimsIn, nDimsOut, fcnType);
    bf.name = name;
    bFcns(d, 1) = bf;
end

S = FunctionSpace(bFcns);
S.name = fsName;
S.plotFcn = plotFcn;
S.defaultPlotX = plottingCoords;

end

function plottingCoords = parsePlottingCoords(fcnType, nDimsIn, fcnParams, P)

if isempty(P.plottingCoords)
    
    nGrid = P.plottingNGrid;
    lims = [];
    
    if isempty(P.plottingLims)
        switch lower(fcnType)
            case {'vm', 'gaus', 'gaus2', 'rcos'}
                mu = fcnParams{1};
                lims = [min(mu); max(mu)];
                lims = num2cell(lims, 1);
            case 'rcoslog'
                mu = fcnParams{1};
                s = fcnParams{3};
                muLin = exp(mu-s);
                lims = {[0, max(muLin)*1.5]};
            case 'conv'
                x = zeros(nGrid, 1);
                idx = ceil(nGrid/2);
                x(idx) = 1;
                plottingCoords = {x}; return;
            case 'conv-causal'
                x = zeros(nGrid, 1);
                x(1) = 1;
                plottingCoords = {x}; return;
        end
    else
        lims = P.plottingLims;
    end
    
    for d = 1:nDimsIn
        plottingCoords{d} = linspace(lims{d}(1), lims{d}(2), nGrid)';
    end
else
    plottingCoords = P.plottingCoords;
end

end