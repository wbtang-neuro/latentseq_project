function P = glmParams()


% Time bin size used for all variables
P.glmBinWidth = 0.010;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COVARIATE BASIS EXPANSION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.basisDir = "\\forskning.it.ntnu.no\ntnu\mh-kin\moser\richarga\misc\to_abraham\misc\model_common_data";

% Basis functions
%
P.posRange = {[0 1.5], [0 1.5]};
P.basisPosSpacing = 0.1;            % tracking data units
P.basisPosSigma = 0.02;             %
P.basisPosBoxPadding = 0.3;         % scales with basis function width, not grid size
P.basisPosDecompGridStep = 0.0025; % finer than necessary, load results rather than recomputing
% P.basisPosDecompGridStep = 0.01; % if need to recompute

P.speedThresh = 0.05;              % m/s
P.basisSpeedN = 50;
P.basisSpeedFirstPeak = P.speedThresh;
P.basisSpeedLastPeak = 3;           % m/s
P.basisSpeedSigma = 0.05;
P.basisSpeedDecompMu = 0.15;
P.basisSpeedDecompSigma = 0.8;

P.basisAngN = 50;
P.basisAngKappa = 10;

P.basisPostSpikeN = 16;
P.basisPostSpikeLastPeak = 0.300;
P.basisPostSpikeB = 0.02;              % log->lin as B 0->inf
P.postSpikeBinWidth = 0.0001;

P.posInterpGridSpacing = 0.01;
P.posInterpGridPadding = 2.5;

P.pcaThresh = struct( ...
    'pos',   .99, ...
    'theta', .80, ...
    'hd',    .80, ...
    'sd',    .80);

[P.version, P.versionDate] = shiftModelVersion();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.decodePosSigmaT = 0.015;             % seconds
P.decodePosBoxPadding = 0;
P.decodePosGridSpacing = 0.02;

P.decodeHdSigmaT = 0.015;
P.decodeAngGridN = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MISCELLANEOUS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.floatClass = 'double';
P.posSmSpeed = 0.5;                 % seconds

P.nSpikesMin = 200; % unit must have at least this quantity in ALL CV folds
P.meanRateMin = 0.2;

% Crossvalidation
P.cvNFold = 10;
P.cvNDiv = 10;

P.posRateMapSmooth = 2;
P.posRateMapN = 31;
P.posRateMapNSpatialInfo = 10;

P.angTuningCurveSmooth = 1;
P.angTuningCurveN = 31;

P.sdMinThetaNSpikes = 10;
P.sdMinThetaMvl = 0.4;
P.sdMinHdMvl = 0.2;
P.sdNPcs = 5;

P.vnames = ["pos", "hd", "sd", "theta_mua"]; % postSpike isn't listed because it's special
P.vcirc = [0, 1, 1, 1];

P.modelNames = ["basic_no_shift", "hd_no_shift", "sd_no_shift", "all_no_shift", "all_sd_shift"];

P.modelVariables.basic_no_shift = ["theta_mua", "pos", "postSpike"];
P.modelVariables.hd_no_shift = ["theta_mua", "pos", "postSpike", "hd"];
P.modelVariables.sd_no_shift = ["theta_mua", "pos", "postSpike", "sd"];
P.modelVariables.all_no_shift = ["theta_mua", "pos", "postSpike", "hd", "sd"];
P.modelVariables.all_sd_shift = ["theta_mua", "pos", "postSpike", "hd", "sd"];

P.dt = P.glmBinWidth;

P.rootDataDir = "N:\richarga\misc\to_abraham\data\sfn2019_3";

end