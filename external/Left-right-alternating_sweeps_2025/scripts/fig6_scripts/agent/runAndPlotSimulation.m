function [mdl, optimSweepDirs, scores, scoresM, F, Fpen] = ...
    runAndPlotSimulation(ax, xy, dt, varargin)
inp = inputParser();

% These parameters are forwardded to createDefaultModel(). If left blank here,
% they'll be ignored and the defaults in createDefaultModel() will be used.
inp.addParameter("posGridRange", []);
inp.addParameter("posGridStep", []);
inp.addParameter("modules", []);
inp.addParameter("params", struct());

% These params are used within this function
inp.addParameter("sweepInds", []);
inp.addParameter("sweepDirs", []);
inp.addParameter("contours", "none");
inp.addParameter("randSeed", 0);
inp.addParameter("plotType", "basic_mono");
inp.parse(varargin{:});
P = inp.Results;

xp = xy(:,1);
yp = xy(:,2);

if isstruct(P.params)
    Pfwd = struct2pvpIgnoreEmpty(P, ["modules", "params", "posGridRange", "posGridStep"]);
    mdl = createDefaultCoverageModel(Pfwd{:});
else
    % % "params" contains a configured model instance: 
    % % take a copy and use directly
    mdl = P.params;
    mdl = mdl.copy();
    if ~isempty(P.posGridRange)
        mdl.createSquarePosGrid(P.posGridStep, P.posGridRange);
    end
end

mdl.axes = ax;
mdl.plotType = P.plotType;

if isempty(P.sweepInds)
    nSweeps = numel(xp);
else
    nSweeps = numel(P.sweepInds);
end

mdl.plotSweepHistoryLength = nSweeps;

if mdl.sweepProfileType == "gaussian"
    contourLevel = 0.8;
elseif mdl.sweepProfileType == "vm_invsqd"
    contourLevel = 0.1;
end

if P.contours == "final"
    mdl.contourLevels = contourLevel;
elseif P.contours == "all"
    mdl.contourLevels = contourLevel;
    mdl.plotAllContours = true;
end
% mdl.createSquarePosGrid(P.posGridStep, P.posGridRange);
mdl.initialize();
rng(P.randSeed); % make sure initial sweep angle is the same

scores = struct();
scoresM = struct();
[optimSweepDirs, scores.all, scoresM.all, F, Fpen] = mdl.run(xp, yp, dt, P.sweepInds, P.sweepDirs);

scores.best     = squeeze(min(scores.all,   [], 2));
scores.worst    = squeeze(max(scores.all,   [], 2));
scoresM.best    = squeeze(min(scoresM.all,  [], 2));
scoresM.worst   = squeeze(max(scoresM.all,  [], 2));

end