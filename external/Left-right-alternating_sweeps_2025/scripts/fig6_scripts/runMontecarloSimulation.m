function [sweepDirs, overlapScoresTot, overlapScoresM, altScores] = ...
    runMontecarloSimulation( posxy, mdl0, nreps, randSeed, verbose)

% This function runs many repeats of a self-driving agent simulation,
% to quantify the variability of the selected sweep directions and
% associated results. The initial sweep direction in every run is
% pseudorandom, which is main source of variability. We also apply
% dithering noise to the input position coordinates, to prevent
% the position binning causing quantization artifacts. 

% N.B. because we're using a parfor loop, we need to be careful to ensure
% that workers use different random seeds (otherwise we may get identical
% outcomes). We do this by designating a random seed for the first
% iteration, which is then incremented by the loop counter. This way, each
% iteration is guaranteed to have a distinct random seed.

if nargin < 4 || isempty(randSeed), randSeed = 0; end
if nargin < 5, verbose = false; end

parfor i = 1:nreps
    if verbose
        fprintf("Rep. %u of %u...\n", i, nreps);
    end

    % Set the random seed for this iteration
    rng(randSeed+i);
    mdl = mdl0.copy();
    mdl.initialize();

    % Apply dithering to position coords. This shifts the position
    % trajectory by a random offset up to half a position bin from the
    % original value.
    dgrid = diff(mdl.posGridX([1, 2]));
    xoff = (rand() - 0.5)*dgrid;
    yoff = (rand() - 0.5)*dgrid;
    x = posxy(:,1) + xoff;
    y = posxy(:,2) + yoff;

    % Run the simulation
    [sweepDirs(:, :, i), ~, overlapM(:,:,:,i)] = mdl.run(x, y);
end

nt = size(posxy, 1);
ndirs = size(sweepDirs, 2);
altScores = zeros(nt, ndirs, nreps);
% calculate alternation score for each separate direction
for i = 1:nreps
    for d = 1:ndirs
        altScores(:, d, i) = alternationScore4(sweepDirs(:, d, i));
    end
end

overlapTotal = squeeze(sum(overlapM, 3));
overlapScoresTot = calcOverlapScores(overlapTotal);
if ndirs==1 % only include alt. score if there is a sensible "total"
    overlapScoresTot.alternation = squeeze(altScores);
end

nModules = size(overlapM, 3);
for m = 1:nModules
    overlapOneModule = squeeze(overlapM(:,:,m,:));
    overlapScoresM(m) = calcOverlapScores(overlapOneModule);
end

% % extract summary scores for all runs
% % (dims are [time, iter], or [time, iter, direction] for "scores.all")
% scores = struct();
% scores.all = permute(overlap, [1, 3, 2]);
% scores.best = squeeze(min(overlap, [], 2));
% scores.worst = squeeze(max(overlap, [], 2));
% scores.diff = scores.worst - scores.best;
% scores.alternation = altScores;

end

function scores = calcOverlapScores(overlap)
scores = struct();
scores.all = permute(overlap, [1, 3, 2]);
scores.best = squeeze(min(overlap, [], 2));
scores.worst = squeeze(max(overlap, [], 2));
scores.diff = scores.worst - scores.best;
end