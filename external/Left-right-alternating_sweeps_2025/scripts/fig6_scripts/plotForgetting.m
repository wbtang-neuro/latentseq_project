function mdl = plotForgetting(P)

arguments
    P.sweepTraceDecayMode = "exponential"
    P.tau = 1; % exponential forgetting factor (range 0 - 1)
    P.discreteMemoryLength = 1;
    P.modules = []
    P.plotType = "basic_mono"
    P.clipPercentile = 99.8
    P.scale = 1 % scale factor for path coordinates and axis limits
end

% Create a bendy self-intersecting locomotor path
xy = simulateLocomotorPath(100);
xy = xy*P.scale;

% define evenly spaced points on the path where sweeps will occur
mdl = createDefaultCoverageModel();
if mdl.sweepProfileType=="gaussian"
    % with the gaussian footprint type, the weight is concentrated further
    % away from the sweep origin, so we need to adjust the sweep positions
    % to compensate for this.
    isweep = 5 : 15 : 90;
else
    isweep = 11 : 18 : 90;
end

posGrid = -1 : 0.005 : 1;

modelParams = struct( ...
    "sweepTraceDecayMode", P.sweepTraceDecayMode, ...
    "tau",      P.tau, ...
    "sweepTraceDiscreteN", P.discreteMemoryLength, ...
    "posGridX", posGrid, ...
    "posGridY", posGrid, ...
    "traceNormalization", "trace", ...
    "clipPercentile", P.clipPercentile, ...
    "addFinalSweepToTrace", false, ...
    "arrowLength", 0.15*P.scale);

mdl = runAndPlotSimulation(gca, xy, 1, ...
    "sweepInds",isweep, ...
    "params",   modelParams, ...
    "modules",  P.modules, ...
    "plotType", P.plotType, ...
    "contours", "final");

h = mdl.plotHandles;
set([h.contourAll{:}], 'faceColor', 'none', 'edgeAlpha', 0.5);
% delete([h.sweepPast.real, h.sweepPast.sim, h.sweep.real, h.sweep.sim]);
plot(xy(isweep, 1), xy(isweep, 2), 'r.', "markerSize", 20);
axis([-0.45, 0.4, -0.6, 0.35]*P.scale);
clim([0, 1]);

end


function xy = simulateLocomotorPath(nt)
% Simulate a locomotor path that loops back on itself
rad = 0.3;
w = 2*pi * linspace(0.18, 0.816, nt)';
y = rad*cos(w);
x = -0.5 * rad*sin(2*w);

% Reinterpolate the path so that the "speed" is constant.
%
% First, calculate the scalar distance between adjacent time steps
displacement = [0; hypot(diff(x), diff(y))];
cumDisplacement = cumsum(displacement);
pathLength = cumDisplacement(end);

% Finally, interpolate the positions of points with equal spacing along the
% path.
cdi = linspace(0, pathLength, nt)';
xy = interp1(cumDisplacement, [x, y], cdi);

end