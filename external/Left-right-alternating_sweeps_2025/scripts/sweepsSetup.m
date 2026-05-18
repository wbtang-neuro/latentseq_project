% Sweeps analysis setup script
clearvars -except useGpu
clc
delete(gcp("nocreate"));

S = SweepsSettings();

% Set some default display options
set(0, "defaultFigureColormap", S.colormap());
set(0, "defaultFigureColor", "w");
set(0, "defaultFigureRenderer", "painters");
set(0, "defaultTextInterpreter", "none");
set(0, "defaultFigureWindowStyle", "docked");
set(0, "defaultAxesNextPlot", "add");
set(0, "defaultPolaraxesNextPlot", "add")
set(0, "defaultScatterMarkerFaceColor", "flat");
set(0, "defaultScatterMarkerEdgeColor", "none");
set(0, "defaultLegendAutoUpdate", "off");

fprintf("Using code root '%s'\n", S.codeRoot());
fprintf("Using data root '%s'\n", S.dataRoot());