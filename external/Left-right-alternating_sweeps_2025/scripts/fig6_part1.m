% Fig 6 / Extended Data Fig. 12
% Part 1: artificial agent simulations

sweepsSetup;

%% Plot module footprints (ED fig. 12a)

plotModuleFootprints();

%% Illustrate scale-free footprint model (ED fig. 12b)

plotScaleFreeFootprintModel();

%% Illustrate self-driving and empirical simulation modes (Fig. 6a,d)

plotSimulationTypes();

%% Illustrate forgetting (ED fig. 12 e)

figure();

% Exponential forgetting
nexttile();
plotForgetting(sweepTraceDecayMode="exponential", tau=1);
title("Tau = 1 (no forgetting)");

nexttile();
plotForgetting(sweepTraceDecayMode="exponential", tau=0.95, clipPercentile=99.95);
title("Tau < 1 (exponential forgetting)");


%% Extended linear trajectory with monte-carlo repeats (Fig. 6b,c)

% Enable monte-carlo repeats by setting the 'monteCarloReps' parameter 
% to a large value, e.g. 1000.

pos = createLinearTrajectory(20);

[~, mdl0] = plotExampleLtTraversal(pos, ...
    modules="off", ...
    footprintType="vm_invsqd", ...
    arrowLength=0.2, ...
    plotType="basic", ...
    monteCarloReps=0);

%% Show the same linear-path simulation as an animation

pos = createLinearTrajectory(20);
plotSweepSimAnimationLT(pos, mdl0);

%% Multi-module version (ED Fig. 12i)

% Enable monte-carlo repeats by setting the 'monteCarloReps' parameter 
% to a large value, e.g. 1000.

pos = createLinearTrajectory(20);

moduleOptions = ["coordinated", "independent", "independent"];
updatingOptions = ["parallel", "parallel", "serial"];

for i = 1:3
    plotExampleLtTraversal(pos, ...
        modules=moduleOptions(i), ...
        footprintType="gaussian", ...
        arrowLength="peak", ...
        plotType="basic", ...
        contours="all", ...
        traceUpdateMode=updatingOptions(i), ...
        monteCarloReps=0);
    sgtitle(sprintf("Modules=%s, updating=%s", moduleOptions(i), updatingOptions(i)), interpreter="none");
end


%% Illustrate effect of changing kappa (ED Fig. 12d)

xy = createLinearTrajectory(3);

for k = [1, 3, 10]
    plotExampleLtTraversal(xy, ...
        kappa=k, ...
        arrowLength=0.2, ...
        randSeed=0, ...
        lims = [-0.5, 0.5, -0.3, 0.3]);
    sgtitle("kappa = " + k);
end

%% Linear movement with discrete forgetting (ED Fig. 12f)

xy = createLinearTrajectory(10);

mdl = plotExampleLtTraversal(xy, ...
    contours="none", ...
    sweepTraceDecayMode="discrete", ...
    sweepTraceDiscreteN=1, ...
    arrowLength=0.15);

%% Monte-carlo linear-path simulation (Fig. 6c, ED Fig. 12c)
% (GPU recommended)

montecarloLtSimulation(1000);

%% Monte-carlo linear-path simulation with different dispersion values (ED Fig. 12c)

% (GPU recommended)
montecarloLtSimulationDispersion(nKappaSteps=20, nReps=100);
