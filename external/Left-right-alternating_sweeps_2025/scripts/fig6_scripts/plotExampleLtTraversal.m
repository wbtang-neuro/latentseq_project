function [mdl, mdl0, overlapScoresM, altScores] = ...
    plotExampleLtTraversal(xy, P)

arguments
    xy
    P.sweepTraceDecayMode = "exponential"
    P.traceUpdateMode = "parallel"
    P.tau = 1 % exponential forgetting factor
    P.kappa = 5
    P.sweepTraceDiscreteN = 1; % discrete memory length
    P.timeInterval = 1/8 % this is only used when tau < 1
    P.randSeed = 3
    P.modules = "off"
    P.footprintType = "vm_invsqd"
    P.plotType = "basic_mono"
    P.contours = "none"
    P.arrowLength = 0
    P.lims = [-1.6, 1.9, -0.32, 0.32]
    P.monteCarloReps = 0
    P.parent = []
end

if isempty(P.parent)
    P.parent = figure();
end

if get(P.parent, 'type') == "figure"
    fig = P.parent;
    P.parent = tiledlayout("flow", "Parent", fig);
end

ax = nexttile(P.parent, [1, 4]);

mdl = createDefaultCoverageModel(footprintType=P.footprintType, modules=P.modules);

mdl.sweepTraceDecayMode = P.sweepTraceDecayMode;
mdl.sweepTraceDiscreteN = P.sweepTraceDiscreteN;
mdl.tau = P.tau;
mdl.kappa = P.kappa;
mdl.traceUpdateMode = P.traceUpdateMode;
mdl.addFinalSweepToTrace = true;
mdl.arrowLength = P.arrowLength;

mdl0 = mdl.copy(); % return uninitialized model

[mdl, ~, ~, scoresM] = ...
    runAndPlotSimulation(ax, xy, P.timeInterval, ...
    "params",           mdl, ...
    "randSeed",         P.randSeed, ...
    "contours",         P.contours, ...
    "plotType",         P.plotType);

if P.plotType == "basic_mono"
    ax.CLim = [0, 3];
end
axis(P.lims);

if P.monteCarloReps==0
    overlapScoresM = [];
    altScores = [];
    return;
end

delete(gcp("nocreate"));
parpool threads
[~, overlapScoresTot, overlapScoresM, altScores] = runMontecarloSimulation(xy, mdl0, P.monteCarloReps, [], true);
delete(gcp);

nexttile(P.parent);

% Total overlap scores
plotOverlapScores(overlapScoresTot.best, overlapScoresTot.worst);
title("Total overlap scores");

% Individual module score breakdown
if ~strcmpi(P.modules, "off")
    for m = 1:numel(overlapScoresM)
        nexttile();
        plotOverlapScores(overlapScoresM(m).best, overlapScoresM(m).worst);
        title("Overlap M" + m);
    end
end

% Alternation scores (with conf. intervals)
nexttile();
plotAltScores(altScores, "errorbar");
title("Alt. scores")

% Alternation scores (without conf. intervals)
nexttile();
plotAltScores(altScores, "line");
title("Alt. scores")

end

function h = plotAltScores(altScores, plotStyle)
cols = {[0.8, 0, 0], [0, 0.7, 0], [0, 0, 1]};
ndirs = size(altScores, 2);
for d = 1:ndirs
    y = squeeze(altScores(:, d, :));
    h(d) = plotOneScore(y, cols{d}, plotStyle);
end
legend(h, "Dir #" + (1:ndirs));
axis square
xlabel("Sweep #");
ylabel("Score");
ylim([0, 1]);
end

function plotOverlapScores(scoresBest, scoresWorst)

cols = get(0, 'defaultAxesColorOrder');
hbest = plotOneScore(scoresBest, cols(1, :), "errorbar");
hworst = plotOneScore(scoresWorst, cols(2, :), "errorbar");

axis square
xlabel("Sweep #");
ylabel("Overlap");
legend([hbest, hworst], ["Best choice", "Worst choice"]);
end

function [hline, herr] = plotOneScore(y, col, style)
x = (1:size(y, 1))';
ymu = mean(y, 2);
yci = prctile(y', [5, 95])';

if style=="errorbar"
    eneg = ymu-yci(:, 1);
    epos = yci(:, 2)-ymu;
    herr = errorbar(x, ymu, eneg, epos, "color", col, "lineStyle", "none", 'capSize', 2);
    hline = plot(x, ymu, "-o", "color", col, 'markerSize', 3);
elseif style=="fill"
    v = ~isnan(ymu);
    xv = x(v);
    ylo = yci(v,1);
    yhi = yci(v,2);
    herr = fill([xv; flip(xv)], [ylo; flip(yhi)], col, "faceAlpha", 0.5, "edgeColor", "none");
    hline = plot(x, ymu, "color", col, 'lineWidth', 1.5);
elseif style=="line"
    hline = plot(x, ymu, "-o", "color", col, 'markerSize', 3);
    herr = [];
else
    error("Invalid plotting style");
end

end
