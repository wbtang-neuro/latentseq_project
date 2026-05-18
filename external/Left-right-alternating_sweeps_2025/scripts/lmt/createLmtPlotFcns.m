function fcns = createLmtPlotFcns(D, models, preset, nsmooth)
if nargin < 3, preset = "standard"; end
if nargin < 4 || isempty(nsmooth), nsmooth = 0; end

D.models = models;
D.nsmooth = nsmooth;
% D.smoothX = 0.03/D.dt;

preset = lower(preset);
fcns = {};
switch preset
    case "standard"
        fcns{end+1} = @(fig, mdl) lmt2DPosPlots(fig, mdl, "x_time", D);
%         fcns{2} = @(fig, mdl) lmt2DPosPlots(fig, mdl, "x_2d", D);
        fcns{end+1} = @(fig, mdl) basicPlots(fig, mdl);
        fcns{end+1} = @(fig,mdl)  unitPlot(fig, mdl, "pos",[], [], D);
        fcns{end+1} = @(fig,mdl)  unitPlot(fig, mdl, "tc", [], [], D);
    case "circ"
        fcns{end+1} = @(fig, mdl) lmt2DPosPlots(fig, mdl, "x_time", D);
        fcns{end+1} = @(fig, mdl) basicPlots(fig, mdl);

    case "pos"
        fcns{end+1} = @(fig, mdl) lmt2DPosPlots(fig, mdl, "x_2d", D);
        fcns{end+1} = @(fig, mdl) lmt2DPosPlots(fig, mdl, "f_peaks", D);

end


end