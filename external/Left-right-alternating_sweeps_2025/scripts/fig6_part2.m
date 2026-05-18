% Fig 6 / Extended Data Fig. 12
% Part 2: comparisons between empirical data and artificial agent simulation

%% Load batch-run data
sweepsSetup
S.dataRoot_ = "~/Temp/sweeps/data"; % set this to your local data path

% Load the exported ID data source file
tmp = load(fullfile(S.dataRoot, "covmdl4", "fig6_id_data_v2.mat"));
batchMdlData = tmp.data;
nrecs = numel(batchMdlData);
batchMdlDataDir = S.filepath("covmdl4", "runs");

recsWithoutAlternatingId = [26035, 26018, 27764];
recHasAlternation = ~ismember([batchMdlData.animal_ID], recsWithoutAlternatingId);

% Sweep directions are most predictable when the animal moves at high
% speeds, so we focus on high-speed periods when running the model
% parameter search.
speedThresh = 0.2;

%% Load and plot parameter fits for each recording (ED Fig. 12e)

clear bestScores tauBest kappaBest batchDataAll

% initialize plots
cmap = seismic();
figure("colormap", cmap); layoutCorr = tiledlayout(4, 4); sgtitle("ID correlation");
figure("colormap", cmap); layoutAltScore = tiledlayout(4, 4); sgtitle("Alt. score");
layouts = {layoutCorr, layoutAltScore};

for f = 1:nrecs

    d0 = batchMdlData(f);

    assert(isscalar(d0));
    fn = fullfile(batchMdlDataDir, string(d0.animal_ID) + ".mat");
    fprintf("%s ...\n", fn);
    data = load(fn);
    batchDataAll(f) = data;

    tauSteps = data.steps.tau;
    kappaSteps = data.steps.kappa;

    vt = d0.speed >= speedThresh;
    id = d0.id;

    % the cell matrices containing angle data have dims [tau, dispersion]
    mdld = data.sweepAnglesDriven;
    mdldFree = data.sweepAnglesFree;
    mdldEgo = cellfun(@(a) {wrapToPi(a-d0.hd)}, mdld);
    idEgo = wrapToPi(d0.id-d0.hd);

    offsets = cellfun(@(a) {wrapToPi(a-id)}, mdld);

    % calculate correlation/error stats using only the data above speed
    % threshold
    corrs = cellfun(@(aego) circ_corrcc(aego(vt), idEgo(vt)), mdldEgo);
    meanErrs = cellfun(@(a) mean(abs(a(vt))), offsets);
    meanAltScoreFree = cellfun(@(a) meanAltScore(a, vt), mdldFree);

    z = {corrs, meanAltScoreFree};
    for n = 1:2

        % find best tau and dispersion for this statistic
        [bestScores(f, n), idxbest] = max(z{n}, [], "all");
        [iTauBest, iKappaBest] = ind2sub(size(meanErrs), idxbest);
        tauBest(f, n) = tauSteps(iTauBest);
        kappaBest(f, n) = kappaSteps(iKappaBest);

        % Plot the score matrix
        ax = nexttile(layouts{n});
        imagesc(ax, z{n}');
        axis(ax, "image");
        xticks(ax, []);
        yticks(ax, []);
        zmx = max(z{n}(:));
        plot(ax, iTauBest, iKappaBest, "bx", "markerSize", 10, "lineWidth", 2);
        xlabel(ax, "\tau", "interpreter", "tex");
        ylabel(ax, "\kappa", "interpreter", "tex")
        title(ax, d0.animal_ID);
        if n==1
            clim(ax, [0, 0.63]);
        else
            clim(ax, [-1, 1]);
        end

        colorbar

    end
    
    % Get the sweep directions for the self-driving agent, using the
    % tau and kappa values with the highest alternation score.
    batchMdlData(f).idmdlfree = mdldFree{iTauBest, iKappaBest};

    drawnow();

end

clear corrs meanErrs meanAltScoreFree

%% Calculate the average best kappa and tau values for the empirically driven agent

% find the average best tau and kappa values across all animals

% forgetting (tau)
optimTau = median(tauBest(:, 1));        % median best tau value
[~,iTau] = min(abs(tauSteps-optimTau));  % index of the best tau value
optimTau = tauSteps(iTau);

% dispersion (kappa)
optimKappa = median(kappaBest(:, 1));
[~,iKappa] = min(abs(kappaSteps-optimKappa));
optimKappa = kappaSteps(iKappa);

% get the chosen sweep directions from the selected model, for each recording
for r = 1:nrecs
    batchMdlData(r).idmdl = batchDataAll(r).sweepAnglesDriven{iTau, iKappa};
end

fprintf("Median values: tau %.3f, kappa %.3f\n", optimTau, optimKappa);

%% Print results for ED legend

idx = recHasAlternation;
optimtypes = ["Alternation", "ID corr"];
for i = 1:2
    disp(optimtypes(i)+":")
    muforgetting = mean(tauBest(idx, i));
    semforgetting = std(tauBest(idx, i))./sqrt(numel(tauBest(idx, i)));
    mudispersion = mean(kappaBest(idx, i));
    semdispersion = std(kappaBest(idx, i))./sqrt(numel(kappaBest(idx, i)));
    fprintf("Mean forgetting: %.3f + %.3f\n Mean dispersion; %.3f +%.3f\n", muforgetting, semforgetting, mudispersion, semdispersion);
end

%% Plot 2D histograms of ego directions (empirical ID vs. model sweep dir)

figure()
tl=tiledlayout(4, 4);
nrecs = numel(batchMdlData);

clear allStats

for r = 1:nrecs
    res = batchMdlData(r);
    vt = res.speed>.2;
    egoid_mdl = circ_dist(res.idmdl, res.hd);
    egoid = circ_dist(res.id, res.hd);
    egoid_mdl(~vt) = nan;
    egoid(~vt) = nan;
    stats = egoego_heatmap(ax=nexttile,x=egoid, y=egoid_mdl, computeStats=true, labels=[sprintf("Empirical\nEgo angle"),sprintf("Ego angle\nModel")]);
    [stats.coAlternation, stats.coAlternationPvals(r)] = computeCoAlternation([egoid, egoid_mdl]);
    offset = circ_dist(res.idmdl, res.id);
    stats.offsetMu(r) = circ_mean(offset(vt));
    stats.offsetStd(r) = circ_std(offset(vt));
    title(res.animal_ID);

    allStats(r) = stats;
end

%% Calculate summary statistics 

% (exclude datasets without clear alternating empirical ID)
stats = allStats(recHasAlternation);

sem = @(x) std(x)/sqrt(numel(x));

% Run stats tests
muFractionSameSide = mean([stats.psame])
semFractionSameSame = sem([stats.psame])
muCorr = mean([stats.corr])
semCorr = sem([stats.corr])
dirOffset = rad2deg(circ_mean([stats.offsetMu]'))
dirSem = rad2deg(circ_std([stats.offsetMu]')) / sqrt(sum(recHasAlternation))

%% Example ego direction polar histograms at different speeds (Fig. 6f)

irec = find([batchMdlData.animal_ID]==25843);
speedRanges = {[.05, .15], [.15, .30], [.30, .45], [.45, Inf]};
nSpeedRanges = numel(speedRanges);

figure();
layout = tiledlayout(2, nSpeedRanges+1);

data = batchMdlData(irec);
labels = {"Simulated", "Empirical"};
cols = {S.col_covmodel, S.col_id};
dirData = {data.idmdlfree, data.id};

for d = 1:2
    egoDirs = wrapToPi(dirData{d} - data.hd);
    nexttile
    text(0.5, 0.5, labels(d), 'sc', ...
        "HorizontalAlignment","center", ...
        "VerticalAlignment","middle", ...
        "FontSize", 14);
    axis off
    for r = 1:numel(speedRanges)
        speedRange = speedRanges{r};
        validSpeed = data.speed>speedRange(1) & data.speed<=speedRange(2);
        simplePolarAxes(layout);
        polarhistogram(egoDirs(validSpeed), S.gv.angular, "FaceColor", cols{d});
        title(sprintf("%.0f-%.0f cm/s", speedRange(1)*100, speedRange(2)*100))
    end
end

%% Alternation strength vs. running speed (Fig. 6g)

% Set up speed-binning that we'll use in some analyses
xEdgeData.speed = (0.05 : 0.05 : 0.75)';
xEdgeData.straightness = (0 : 0.05 : 1)';

xLabelData.speed = "Running speed (cm/s)";
xLabelData.straightness = "Path straightness";

dirTypes = ["empirical", "simulated"];
cols = struct("empirical", S.col_id, "simulated", S.col_covmodel);
minPoints = 50; % require at least this many observations (time points) per x-bin
alpha = 0.05;

clear res
for xVarname = ["speed", "straightness"]

    edges = xEdgeData.(xVarname);
    bins = edg2cen(edges);
    nbins = numel(bins);

    clear meanAltScores corrs pvals
    for dtyp = dirTypes
        meanAltScores.(dtyp) = nan(nrecs, nbins);
        corrs.(dtyp) = nan(nrecs, 1);
        pvals.(dtyp) = nan(nrecs, 1);
    end

    for r = 1:nrecs

        D = batchMdlData(r);

        if xVarname=="speed"
            xValues = D.speed;
        else
            window = [10, 10];
            xValues = 1 ./ movingStraightness([D.x, D.y], window);
        end

        dirValues.empirical = D.id;
        dirValues.simulated = D.idmdlfree;

        for dtyp = dirTypes
            altScoreValues.(dtyp) = alternationScore4(dirValues.(dtyp));
        end

        % calculate mean alternation scores in each x-data bin (for both empirical and model angles)
        for s = 1:nbins
            thrL = edges(s);
            thrU = edges(s+1);
            validTimeBin = xValues>=thrL & xValues<thrU;
            if sum(validTimeBin) > minPoints
                for dtyp = dirTypes
                    altScores = altScoreValues.(dtyp)(validTimeBin);
                    meanAltScores.(dtyp)(r, s) = mean(altScores, "omitnan");
                end
            end
        end

        for dtyp = dirTypes
            % Correlate the X bin values with the binned alternation score in the bin
            y = meanAltScores.(dtyp)(r, :)';
            validBin = ~isnan(y);
            x = bins(validBin);
            y = y(validBin);
            [corrs.(dtyp)(r), pvals.(dtyp)(r)] = corr(x, y);
        end

    end
    
    for dtyp = dirTypes
        fprintf("%s vs alt. score (%s):", xVarname, dtyp)
        r = corrs.(dtyp)(recHasAlternation);
        p = pvals.(dtyp)(recHasAlternation);
        calcMeanSemCorrelation(r, p, alpha);
    end

    % Make line plot showing relationship for each animal
    % nrecs = size(recs, 1);
    figure();
    for dtyp = dirTypes
        nexttile();
        y = meanAltScores.(dtyp)(recHasAlternation, :)';
        plot(bins, y, ".-");
        yline(S.alternation_score_chance_level, "k");
        xlabel(xLabelData.(xVarname));
        ylabel("Alternation score");
        ylim([0.3, 1.0]);
        title(dtyp);
    end
    sgtitle(xVarname);
    meanAltScores.bins = bins;
    res.(xVarname) = meanAltScores;
end

%% Functions

function meanScore = meanAltScore(dirs, vt)
scores = alternationScore4(dirs);
scores = scores(vt);
meanScore = mean(scores, "omitnan"); % first and last values are always NaN
end

function [nSignificant, rmean, rsem] = calcMeanSemCorrelation(corrs, pvals, alpha)
% Display mean and SEM correlation r-value for a set of correlations

% Count the number of significant correlations
isSignificant = pvals<alpha;
nSignificant = sum(isSignificant);
nTotal = numel(corrs);
fprintf("Significant in %u/%u animals\n", nSignificant, nTotal);

% Mean and SEM of r-values where significant
rmean = mean(corrs(isSignificant));
rsem = std(corrs(isSignificant))./sqrt(nSignificant);
fprintf("r=%.3f SEM=%.3f\n\n", rmean, rsem)
end