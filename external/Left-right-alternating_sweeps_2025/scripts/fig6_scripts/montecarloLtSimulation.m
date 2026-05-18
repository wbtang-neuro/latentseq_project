% Run LT simulation many times with random initial conditions to get 
% distribution of penalties and alternation
% 
% Or maybe repeat this for many kappa values
% change penalty to overlap

function montecarloLtSimulation(nReps)

if nargin < 1, nReps = 1000; end

S = SweepsSettings();

xy = createLinearTrajectory(20);
nt = size(xy, 1);

mdl0 = createDefaultCoverageModel();
delete(gcp("nocreate"));
parpool threads
[sweepDirs, scores] = runMontecarloSimulation(xy, mdl0, nReps, [], true);
delete(gcp);

% Plot scores
figure("tag", "fig_score_errbar");

names = ["best", "worst", "diff", "alternation"];
cols = struct("best", [0, 0, 1], "worst", [0.85, 0, 0], "diff", [0, 0, 0], "alternation", [0, 0, 0]);
tiledlayout(1,3);
nexttile();

clear h
for name = names
    y = scores.(name);
    ymu = mean(y, 2);
    yci = prctile(y', [5, 95])';
    eneg = ymu-yci(:, 1);
    epos = yci(:, 2)-ymu;

    if name == "diff" || name == "alternation"
        nexttile();
    end
    h.(name) = plot(ymu, "-o", "color", cols.(name), 'markerSize', 3);
    errorbar(1:nt, ymu, eneg, epos, ...
        "color", cols.(name), ...
        "lineStyle", "none", ...
        'capSize', 2, ...
        'tag', "errbar_"+name);

    if name == "alternation"
        yline(S.alternation_score_chance_level, 'Color', [.5,.5,.5]);
        yticks(0:.5:1)
        ylim([0,1])
        title("Alt. score");
    elseif name == "diff"
        ylim([0, 0.003]);
        title("Optimality");
    else
        ylim([0, 0.003]);
        title("Overlap");
    end
    xlabel("Sweep number");

    xlim([1, nt]);
    set(gca, 'FontSize',12)
end

legend([h.best, h.worst], ["best", "worst"]);

% Check when alternation becomes significantly higher than chance

figure();
randAltScore = S.alternation_score_chance_level;
edgs = linspace(0,1,30);
h = [];
timeInds = 2:5;
nplt = numel(timeInds);
cols = jet(nplt);

for i = 1:nplt
    tind = timeInds(i);
    col = cols(i, :);
    hist = histcounts(scores.alternation(tind, :), edgs, 'Normalization','probability');
    h(end+1) = area(edg2cen(edgs), hist, 'EdgeColor', col, 'FaceColor', col, 'faceAlpha', 0.3, 'lineWidth', 1);
    [res.pvals(i), res.hypot(i), stats] = signtest(scores.alternation(tind, :)', randAltScore, 'tail','right');
    res.signs(i) = stats.sign;
    res.zvals(i) = stats.zval;
end
xline(randAltScore);
xlabel("Alternation score");
ylabel("Probability");
title("Distribution of alt. scores per time point");
legend(h, "t = "+timeInds);

% Quantify alternation across all trials, at the final valid time point
figure();
plot(timeInds, res.pvals);
yline(.05);
ylabel("Alt. score pval");
xlabel("Time step");
title("Significance of alternation over time steps")

% Mean and sem at final valid timepoint
muscore = mean(scores.alternation(19, :));
stdscore = std(scores.alternation(19, :));
semscore = std(scores.alternation(19, :))./sqrt(nReps);
sprintf("Alternation Mean: %.2f, SD: %.4f, SEM: %.4f", muscore, stdscore, semscore)

end