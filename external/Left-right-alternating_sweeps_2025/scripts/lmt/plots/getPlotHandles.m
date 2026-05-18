function h = getPlotHandles(mdl, plotType, hSet)

persistent allHandles

if ~nargin
    % calling without arguments purges cache
    clear allHandles
    return;
end

if isempty(allHandles)
    allHandles = struct();
end

if ~isfield(allHandles, plotType)
    allHandles.(plotType) = [];
end

h = allHandles.(plotType);
mdlName = mdl.nameString();

if nargin < 3 || isempty(hSet)
    if ~isfield(h, mdlName)
        h.(mdlName) = [];
    end
    h = h.(mdlName);
else
    allHandles.(plotType).(mdlName) = hSet;
end

end