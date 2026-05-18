function montecarloLtSimulationDispersion(P)
% Monte-carlo simulation of linear-track simulation, with different
% dispersion (kappa) params

arguments
    P.nKappaSteps = 20
    P.nReps = 100
end

S = SweepsSettings();

xy = createLinearTrajectory(20);

% Generate a new random seed for incrementing within the MC loop.
mdl0 = createDefaultCoverageModel();

kappaSteps = linspace(1, 10, P.nKappaSteps);

tic
randSeed = 0;
delete(gcp("nocreate"));

parpool threads
for s = 1:P.nKappaSteps
    fprintf("kappa step %u of %u... \n", s, P.nKappaSteps);
    mdl0.kappa = kappaSteps(s);
    [sweepDirs{s}, scores(s)] = runMontecarloSimulation(xy, mdl0, P.nReps, randSeed, true);
end
delete(gcp);

% Plot scores
figure
names = ["best", "worst", "diff", "alternation"];
cols = struct("best", [0, 0, 1], "worst", [0.85, 0, 0], "diff", [0, 0, 0], "alternation", [0, 0, 0]);
tiledlayout("flow");
nexttile();

clear h
for name = names
    y = cat(3, scores.(name));
    y = permute(y, [3, 1, 2]); % rearrange to [dispersion, time, reps]
    % average over timepoints from 10 onwards (stable period) and simulation reps
    ymu = mean(y(:, 10:end, :), [2,3], "omitnan");
    yci = prctile(y(:, 10:end, :), [5, 95], [2,3]);
    if name == "diff" || name == "alternation"
        nexttile();
    end
    h.(name) = plot(kappaSteps,ymu, "-o", "color", cols.(name), 'markerSize', 3);
    eneg = ymu-yci(:, 1);
    epos = yci(:, 2)-ymu;
    errorbar(kappaSteps, ymu, eneg, epos, "color", cols.(name), "lineStyle", "none", 'capSize', 2);
    if name == "alternation"
        yline(S.alternation_score_chance_level, 'Color', [.5,.5,.5]); % TODO: check this hard-coded value is correct
        yticks(0:.5:1)
        ylim([0,1])
        title("Alt. score");
    elseif name == "diff"
        title("Difference");
    else
        title("Overlap");
    end
    xlabel("\kappa", "interpreter", "tex");
    xlim([min(kappaSteps), max(kappaSteps)]);
    set(gca, 'FontSize',12)
end

legend([h.best, h.worst], ["best", "worst"]);

% Get typical directions (make one hist for each dispersion value)
edgs = linspace(-100,100, 80);
cens = edg2cen(edgs);
hall = zeros(numel(cens), P.nKappaSteps);
for s = 1:P.nKappaSteps
    sdirs = sweepDirs{s}(3:end, :); % only consider stable period
    sdirs = rad2deg(sdirs(:));
    hall(:,s) = histcounts(sdirs, edgs, 'Normalization','probability');
end

imagesc(nexttile, hall, 'XData', kappaSteps, 'YData', cens); 
xlabel("\kappa", "interpreter", "tex"), ylabel("Sweep dir. (deg)");
c = colorbar;
axis tight;
set(gca, 'FontSize',12)
ax = gca; 
ax.CLim = [0,.15];
c.Ticks = [0,.15];

