function plotSimulationTypes()

% clear();
sweepsSetup;

npos = 3;

arrowLength = 0.3;
arrowLineArgs = {"lineWidth", 2};

x = linspace(-0.5, 0.5, npos)';
y = zeros(size(x));
empiricalSweepDirs = [0.5, -0.5, 0.4]';
colEmpirical = "r";

figure()
tiledlayout(2, 1);
names = ["Self-driving", "Empirically driven"];

for n = 1:2
    nexttile()
    title(names(n));
    if n==1
        % self-driven mode: we don't supply any empirical directions to the
        % model
        inputSweepDirs = [];
    else
        inputSweepDirs = empiricalSweepDirs;
    end

    % Run the simulation
    [mdl, optimSweepDirs] = runAndPlotSimulation(gca, [x, y], 1, ...
        "sweepDirs", inputSweepDirs, ...
        "params", struct("addFinalSweepToTrace", false, "ditherSweepDirections", false) );

    for i = 1:npos
        if i<npos
            alpha = 0.5;
        else
            alpha = 1;
        end
        if n==1 || i==npos
            agentArgs = {"agentType", "robot", "width", 0.2};
        else
            agentArgs = {"agentType", "rat", "width", 0.4};
        end
        plotAgent(x(i), y(i), "alpha", alpha, agentArgs{:});
    end

    plotArrow(x, y, optimSweepDirs, arrowLength, "color", S.col_covmodel, arrowLineArgs{:});

    if n==2
        plotArrow(x, y, empiricalSweepDirs, arrowLength, "color", colEmpirical, arrowLineArgs{:});
    end

    axis equal off
    xlim([-1.2, 1.2]);
    ylim([-0.6, 0.6]);
    clim([0, 1]);
end

end

function h = plotArrow(x0, y0, direction, arrowLength, varargin)
[u, v] = pol2cart(direction, arrowLength);
x = [x0, x0+u, nan(size(u))]';
y = [y0, y0+v, nan(size(v))]';
h = plot(x(:), y(:), varargin{:});
end