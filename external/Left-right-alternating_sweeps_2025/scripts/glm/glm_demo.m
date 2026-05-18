% This script demonstrates how to fit the GLM-based "shift model" to
% single-unit spiking data.

clear
S = SweepsSettings;
S.dataRoot_ = "~/Temp/sweeps/data"; % change this to the root folder where the data files are saved

%% Load the session dataset file

fn = fullfile(S.dataRoot, "navigation", "of", "25843_1.mat");
load(fn, "Dsession");

sessionType = "of";
arenaDiameter = 1.5;

%% Load the precomputed basis-function data for the GLM

% (Could add further steps illustrating how we calculate the gridded values)

% Option 1: load precomputed data
fn = fullfile(S.dataRoot, "basisfunc_data", sprintf("of_%.1fm.mat", arenaDiameter));
bfuncData = load(fn, "fspaces", "interpD"); % for 1.5 m open field

% % Option 2: compute from scratch
% bfuncData = createBasisFunctionsFromScratch(arenaDiameter);

%% Decompose the GLM input variables

[inputVariables, inputVariablesD] = getDecomposedInputVariables(Dsession, bfuncData);

%% Fit the model

% We need to select a single-unit spike train to use as the model's
% response variable.
%
% The IDs of a few example units are listed below; uncomment one of them to
% test fitting the model on it.

% unitId = "2_1043"; % pure grid cell
% unitId = "1_0090"; % HD cell
unitId = "1_0253"; % ID cell
[spikeCounts, U] = getExampleSpikeTrain(Dsession, unitId);

% For plotting purposes, define coordinates of box corners
bx = 0.75 * [-1, -1, 1, 1, -1];
by = 0.75 * [-1, 1, 1, -1, -1];
boxCornerCoords = [bx(:), by(:)];

% Set the L2 penalty strengths for the two parameter types
penaltyAlpha = 0.01;    % for 'alpha' parameters (shift curve)
penaltyBeta = 0.0001;   % for 'beta' parameters (GLM)

% Define the partitioning of the data for cross-validation (optional).
% If crossvalidation is enabled, the model will be fitted K times, each
% time holding out 1/Kth of the time points (the 'test set') for
% computing the log-likelihood value.
%
k = 0; % number of crossvalidation folds (set to zero to disable crossvalidation; otherwise k must be 2 or more)
crossvalidation = KFoldCV(Dsession.nt, k, 1);

shiftData = struct("pos", inputVariables.pos, "angle", inputVariables.id);

[fitData, scores] = runShiftModelComparisonCv( ...
    spikeCounts, ...
    inputVariablesD, ...
    bfuncData.interpD, ...
    bfuncData.fspaces, ...
    shiftData, crossvalidation, ...
    "floatClass", "double", ...
    "alphaPenalty", penaltyAlpha, ...
    "betaPenalty", penaltyBeta, ...
    "box", boxCornerCoords, ...
    "plot", true, ...
    "display", "iter" );

%% Reconstruct the fitted tuning and compare to the pre-loaded fit

% The above procedure fitted all three versions of the model("null", 
% "const" and "theta"). The output results structure "fitData" contains a 
% separate sub-structure containing the results of each model.

% Here we will retreive the fit data for the full model ("theta"). We then
% look in the "all" subfield, which contains the data from the complete
% fitting to all observations in the input data (i.e. not cross-validated).
fitDataTheta = fitData.theta.all;

% Reconstruct tuning curves from the fitted parameters (basis function weights)
[tuning, ggs, gvs] = reconstructShiftModelFit(fitDataTheta, bfuncData.fspaces, sessionType);

% Plot all of the fitted tuning curves
fig = plotTuning(tuning, gvs, boxCornerCoords);
fig.Name = "New fit";

% Now compare the newly fitted tuning with the original fit from the loaded
% data. N.B. the loaded dataset only contains fits of the full version of
% the model ("theta"). So for comparison we should check it against our new
% fit of the "theta" model too.

tuningOriginal = U.smdl;
fig = plotTuning(tuningOriginal, gvs, boxCornerCoords);
fig.Name = "Original fit";


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS

function [spikeCounts, U] = getExampleSpikeTrain(Dsession, unitId)
Us = Dsession.units.mec;
if nargin == 1
    % select a random unit from module 1
    moduleId = 1;
    gridUnitIds = Dsession.unitAcorrClus.grid(1).unitIds;
    UsGrid = findUnits(Us, gridUnitIds);
    rng(7);
    U = randsample(UsGrid, 1);
    fprintf("Selected random example grid cell '%s' from total %u units in module #%u\n", U.id, numel(UsGrid), moduleId);
else
    U = findUnits(Us, unitId);
end
spikeCounts = spikeIndsToCounts(U.spikeInds, Dsession.nt);
spikeCounts = double(spikeCounts);
end

function [vars, varsD] = getDecomposedInputVariables(Dsession, bfuncData)
% Create a struct containing all of the raw input variables. Each of these
% will subsequently be decomposed using a set of basis functions
vars = struct();
vars.hd = Dsession.hd;
vars.id = Dsession.id;
vars.pos = [Dsession.x, Dsession.y];
vars.theta = [Dsession.theta];
vars = structfun(@double, vars, "uniformOutput", false);

% Decompose each of the input variables with its designated set of basis
% functions
fprintf("Decomposing GLM input variables ... ");
varsD = decomposeGlmData(bfuncData.fspaces, vars, []);
fprintf("done.\n");
end

function bfuncData = createBasisFunctionsFromScratch(arenaDiameter)
% lims = calcPosBinningLimits(arenaDiameter);

% adapted from 'save_smdl_basis_data'
posExtendFactor = 1.6; % extend basis-function grid beyond arena walls
posLimsP = posExtendFactor * arenaDiameter/2*[-1, 1];
posRangeP = diff(posLimsP);

% Create the dimension-reduced basis-function grid
% Configure GLM basis function grids
glmP = glmParams();
glmP.pcaThresh.pos = 0.99; % (default is 99)
glmP.basisPosDecompGridStep = posRangeP / 500; % 0.25-cm spacing (in 1.5 m box)
% glmP.basisPosDecompGridStep = posRangeP / 1000; % 0.25-cm spacing (in 1.5 m box)
glmP.basisPosSpacing = posRangeP / 30; % 30 x 30 grid of BFs
glmP.basisPosSigma = posRangeP / 15; % should give 92 PCA BFs
glmP.basisPosBoxPadding = 0; % don't extend decomposition grid beyond BF grid
fspaces = createDefaultFunctionSpaces(posLimsP'.*[1,1], false, glmP);

% Now evaluate the BFs at a grid of points to create a lookup table for
% interpolating the BF values when we fit the GLM
grid = linspace(posLimsP(1), posLimsP(2), 100);
griddedValues = DRFDecomp(fspaces.pos, grid, grid);

% [fspaces, griddedValues] = createDefaultFunctionSpaces(lims);
bfuncData.fspaces = fspaces;
bfuncData.interpD = griddedValues;
end

function fig = plotTuning(tuning, gvs, boxCornerCoords)

lim = [0, 15];

fig = figure();
tiledlayout("flow");

for varName = ["theta", "hd", "id"]
    nexttile();
    title("GLM: " + varName);
    x = rad2deg(gvs.(varName));
    y = tuning.beta.(varName);
    y = exp(y);
    semiAngularXPlot(x, y);
    ylim(lim);
    ylabel("Rate multiplication factor");
end

nexttile();
x = gvs.pos;
y = gvs.pos;
z = tuning.beta.pos;
z = exp(z);
imagesc(x, y, z);
hold on
plot(boxCornerCoords(:, 1), boxCornerCoords(:, 2), 'w', 'lineWidth', 2);
axis image xy
clim(lim);
title("pos");

% If this model has a shift curve, plot it
nexttile();
if isfield(tuning, "alpha")
    a = tuning.alpha;
    if isstruct(a) && isfield(a, "theta")
        % Original format returned by reconstructShiftModelFit():
        % contributions from intercept and theta-modulated
        % component are represented separately
        doPlot = true;
        y = a.intercept + a.theta;
    elseif isnumeric(tuning.alpha)
        % Format used in saved data: intercept and theta-mod components
        % combined into a single shift curve
        doPlot = true;
        y = a;
    else
        doPlot= false;
    end

    if doPlot
        x = rad2deg(gvs.theta);
        semiAngularXPlot(x, y);
        title("Shift curve");
        xlabel("Theta phase");
        ylabel("Position shift (m)");
        ylim([-0.1, 0.2]);
    end
end

end

function semiAngularXPlot(x, y)
    plot(x, y);
    xlim([-180, 180]);
    xticks([-180, 0, 180]);
    yline(0);
    % xticklabels(["-180", "0", "\pi"]);
end