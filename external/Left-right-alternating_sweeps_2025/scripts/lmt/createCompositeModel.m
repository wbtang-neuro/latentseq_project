function MComp = createCompositeModel(D, units, P)
% 
% INPUTS
% D     - session dataset structure
% units - struct array with elements corresponding to single units. Must
%         contains fields "spikeCounts" and "id",

arguments
    D                   (1,1) struct
    units                     struct

    P.useHd             (1,1) logical         = 1
    P.useId             (1,1) logical         = 1
    P.usePos            (1,1) logical         = 1
    P.usePopRate        (1,1) logical         = 1
    P.useTheta          (1,1) logical         = 1

    P.useGpu            (1,1) logical         = 1

    P.nIter             (1,1) {mustBeNumeric} = 150
    P.XstoreIter        (1,1) logical = 0
    P.plotInterval      (1,1) {mustBeNumeric} = 10
end

Y      = full([units.spikeCounts]);
Y      = single(Y);
if P.useGpu
    Y = gpuArray(Y);
end

% Put the models in order of increasing complexity. May not make much
% difference, but in the first few steps it might help to remove the
% "easiest" variance first.

M = {};
if P.usePopRate,   M{end+1} = createPopRateGlm(Y);        end
if P.useTheta,     M{end+1} = createThetaLMT(Y, D);       end
if P.useHd,        M{end+1} = createDirLMT(Y, D, "hd");   end
if P.useId,        M{end+1} = createDirLMT(Y, D, "id");   end
if P.usePos,       M{end+1} = createPosLMT(Y, D);         end

models = struct();

for i = 1:numel(M)
    mdl = M{i};
    mdl.unitIds = [units.id];
    % mdl.dt = D.dt;
    mdl.F = single(mdl.F);
    if P.useGpu
        mdl.F = gpuArray(mdl.F);
    end
    if isa(mdl, "LmtModel")
        mdl.XstoreIter = P.XstoreIter;
        mdl.binTimes = D.t;     
        mdl.Xinit = single(mdl.Xinit);
        if P.useGpu
            mdl.Xinit = gpuArray(mdl.Xinit);
        end
        mdl.annealIterStart = mdl.XclampIter;

    elseif isa(mdl, "PoissonGLM")
        mdl.X = single(mdl.X);
        mdl.F = single(mdl.F);
        if P.useGpu
            mdl.X = gpuArray(mdl.X);
            mdl.F = gpuArray(mdl.F);
        end
    end
    models.(mdl.name) = mdl;
end


% Create pos plotting functions (requires access to ID model)
if isfield(models, "pos")
    fcns1 = createLmtPlotFcns(D, models, "standard");
    fcns2 = createLmtPlotFcns(D, models, "pos");
    models.pos.plotFcns = [fcns1, fcns2];
end

% Finally assemble the "composite" model from the separate components
MComp = CompositeModel(Y, M, [units.id]);
MComp.niter = P.nIter;
MComp.plotInterval = P.plotInterval;
MComp.initialize();

end

function mdl = createThetaLMT(Y, D)
% Create LMT with "fixed" input variable for theta phase

Xinit = D.theta;
mdl = LmtModel(Y, Xinit);

mdl.isCircular = 1;
mdl.name = "theta";

mdl.hparams.lamff  = 0;       % No L1 penalty
mdl.hparams.rhoff = 1000;
mdl.hparams.lenffR = 2; % penalize smoothness less

mdl.enableXStep = 0;
mdl.annealIterStart = Inf;
mdl.hparamRanges.sigma = [1, 1];       % don't anneal

models.theta = mdl;
mdl.plotFcns = createLmtPlotFcns(D, models, "circ");

end

function mdl = createDirLMT(Y, D, name)
% Create LMT for 1-D circular variable

% initialize the model with HD as the initial latent state
mdl = LmtModel(Y, D.hd);
mdl.isCircular = 1;
mdl.name = name;

mdl.hparams.rhoxx   = 100;
mdl.hparams.rhoff   = 0.1;
mdl.hparams.lenxx   = 0.1;
mdl.hparams.lamff   = 1;
mdl.hparams.lenffR  = 0.5;

mdl.hparamRanges.sigma = [3, 1, 100]; % smaller starting value (hold closer to init values)

mdl.XclampIter = 30;
mdl.enableXStep = name=="id";

models.(name) = mdl;
mdl.plotFcns = createLmtPlotFcns(D, models, "circ");

end

function mdl = createPosLMT(Y, D)
% Create LMT for 2-D position

mdl = LmtModel(Y, [D.x, D.y]);
mdl.name = "pos";

% Specify theta phase range to use when aligning the latent and true
% positions
if ~isempty(D.theta)
    phir = [2, 5]; % reset phase
    phi = mod(D.theta - phir(1), 2*pi) + phir(1);
    vphi = phi > phir(1) & phi < phir(2);
    assert(numel(vphi)==mdl.nt);
    mdl.XAlignInds = vphi;
end

% set hparams and initialize

% Clamp pos X for longer, so that other models can remove residual
% variance before it's allowed to leak into pos X.
mdl.XclampIter = 40;

mdl.hparams.rhoff = 0.1;
mdl.hparams.rhoxx = 10; % 1000
mdl.hparams.lamff = 1; % remove this for single-module stuff
mdl.hparams.lenffR = 6;

mdl.hparamRanges.sigma = [3, 1, 100];
mdl.XboundPercentile = 0.01;

end

function mdl = createPopRateGlm(Y)

% Create the GLM design matrix. There will be two columns (input
% variables): (1) the intercept term, and (2) the smoothed population
% log firing rate.
popRate = log(mean(Y, 2));
popRate = gsmooth(popRate, 3);
popRate = max(popRate, -5);
XGlm = [ones(size(Y, 1), 1), zscore(popRate)];

mdl = PoissonGLM(Y, XGlm);
mdl.name = "glm_poprate";
mdl.optimizeIntercept = 1;

mdl.plotFcns{1} = @(fig, mdl) glmPlots(fig, mdl, "glm");
mdl.plotFcns{2} = @(fig, mdl) basicPlots(fig, mdl);

end