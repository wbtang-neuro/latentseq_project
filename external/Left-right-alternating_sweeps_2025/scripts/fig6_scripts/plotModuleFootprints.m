function plotModuleFootprints(P)
% Draw series of simulated modules EDa

arguments
    P.nModules = 3
    P.gridFieldWidth = 0.15         % field width as fraction of grid spacing
    P.footprintType = "gaussian"

    % Here we model sweep lengths as being a fixed ratio of the grid
    % spacing for each module. The empirical ratio is ~ 0.22, but for
    % illustrative purposes here it helps to use a larger ratio here, so
    % that the sweep offsets can be easily seen.
    P.sweepLenToSpacingRatio = 1/3; 
end

moduleScaleRatio = sqrt(2);     % geometric ratio between spacings of grid modules

moduleSpacings = moduleScaleRatio .^ (0:P.nModules-1);
moduleFieldRadiuses = P.gridFieldWidth*moduleSpacings;
moduleSweepLengths = P.sweepLenToSpacingRatio*moduleSpacings;


lim = max(moduleSpacings) * 0.75;
nPosGrid = 1000;
gy = linspace(-lim, lim, nPosGrid);
gx = linspace(-lim, lim, nPosGrid);
[xx,yy] = meshgrid(gx, gy);

figure();
tiledlayout(2, P.nModules, "TileIndexing", "columnmajor");

moduleWeights = zeros(P.nModules, 1);

% Create empty 2D spaces for accumulating grid fields
z3 = zeros(nPosGrid);
zW = zeros(nPosGrid);

for m = 1:P.nModules

    % plot grid
    nexttile();
    % gmod = rg.spatial.grid.GridModule;
    % gmod.spacing = moduleSpacings(m);
    % gmod.fieldWidth = moduleFieldRadiuses(m);
    % gcell = gmod.zeroPhaseGrid();
    % zgrid = gcell.pdf(gx, gy);
    zgrid = hexGridIntensity(gx, gy);
    zgrid = zgrid./max(zgrid(:));

    zgrid3 = zeros([size(zgrid), 3]);
    zgrid3(:, :,m) = zgrid;

    imagesc(gx, gy, zgrid3);
    axis image xy off
    title("module " + m)
    clim([0, 1]);

    % plot sweep footprint
    % calculate the module footprint (a 2D Gaussian function)
    dy = yy-moduleSweepLengths(m);
    dx = xx-0;
    dd = hypot(dx, dy);
    zmod = normpdf(dd, 0, moduleFieldRadiuses(m));
    zmod = zmod./max(zmod(:));
    w = 1./sum(zmod(:));
    zmodW = zmod*w;
    moduleWeights(m) = w;

    nexttile();
    zmod3 = zeros([size(zmod), 3]);
    zmod3(:, :, m) = zmod;
    imagesc(gx, gy, zmod3);
    axis image xy off

    plot(0, 0, 'r+');
    xline(0, 'w');
    yline(0, 'w');

    z3 = z3 + zmod3;
    zW = zW + zmodW;
end

% plot unnnormalized and normalized sums
figure("colormap", bone());
zdat = {z3, zW};
names = ["unnormalized", "normalized"];
tiledlayout(1, 2)
sgtitle("Summed module footprints");

dirgrid = linspace(-pi, pi, 101);

for n = 1:2
    nexttile();
    imagesc(gx, gy, zdat{n});
    plot(0, 0, 'r+');
    xline(0, 'w');
    yline(0, 'w');
    title(names(n));
    axis image xy off

    % Plot rings marking each constituent grid field
    for m = 1:P.nModules
        [u, v] = pol2cart(dirgrid, moduleFieldRadiuses(m));
        yf = u+moduleSweepLengths(m);
        xf = v;
        col = [0, 0, 0];
        col(m) = 1;
        plot(xf, yf, "color", col);
    end
    cl = prctile(zdat{n}(:), 99);
    clim([0, cl]);
end

% Short behavioral path with sweep traces
% (N.B. this is purely for illustration; the sweep angles are hard-coded,
% not chosen by the agent).
nt = 100;

scale = 0.7;
px = linspace(-scale, scale, nt)';
py = linspace(-scale, scale, nt)';
rng(3); % set random seed for reproducibility
py = py + scale*randn(size(py));
py = gsmooth(py, 10);

sweepInds = [20, 40, 60, 80];
sweepDirs = deg2rad([-20, 100, 0, 120]);

figure("colormap", bone());

for plotType =  ["basic", "basic_mono"]
    nexttile();
    mdl = createDefaultCoverageModel("footprintType", P.footprintType, modules="coordinated");
    mdl.addFinalSweepToTrace = true;
    mdl.sigma = P.gridFieldWidth / P.sweepLenToSpacingRatio;
    mdl.firstModuleSweepLength = moduleSweepLengths(1);

    mdl = runAndPlotSimulation(gca, [px, py], 1, ...
        params = mdl, ...
        sweepInds=sweepInds, ...
        sweepDirs=sweepDirs, ...
        posGridRange=[-2, 2]*scale, ...
        posGridStep=0.01, ...
        plotType=plotType );

    % delete lines indicating optimal sweep angle (not relevant here)
    h = mdl.plotHandles;
    
end

end