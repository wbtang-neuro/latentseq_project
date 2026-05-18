function formatCircAxes(ax, xy, label, varargin)
% Format axes and X/Y data series (line, bar, scatter)

if nargin < 3, label = []; end

S = SweepsSettings();

inp = inputParser();
inp.addParameter("nreps", 2);
inp.addParameter("outputUnits", "auto");
% inp.addParameter("label", []);
inp.addParameter("limits", []); % specify in OUTPUT units
inp.addParameter("interval", []);
inp.addParameter("string", []);
inp.addParameter("center", []);

inp.parse(varargin{:});

P = inp.Results;

if nargin < 2 || isempty(xy), xy = "both"; end

outUnits = lower(P.outputUnits);
xy = lower(xy);

if strcmpi(outUnits, "auto")
    outUnits = S.angleUnits;
end

[inUnits, axFmt, axNoFmt, dataObjs, propsFmt, propsNoFmt] = parseAxes(ax, xy);

if inUnits == "radians"
    toRadFcn = @(x) x;
elseif inUnits == "degrees"
     toRadFcn = @deg2rad;
end

if outUnits == "radians"
    lim = [0, 2*pi];
    cen = pi;
    interval = pi;
    if ~isempty(label)
        label = sprintf("%s (rad.)", label);
    end
elseif outUnits == "degrees"
    lim = [0, 360];
    cen = 180;
    interval = 180;
    if ~isempty(label)
        label = sprintf("%s (°)", label);
    end
end

if isempty(P.center), P.center = cen; end
lim = lim - mean(lim) + P.center;

if isempty(P.limits), P.limits = lim; end
if isempty(P.interval), P.interval = interval; end

% % Axis ticks
ticks = P.limits(1)*P.nreps : P.interval : P.limits(2)*P.nreps;

% % Tick label
if outUnits == "radians"
    num = ticks / pi;
    tlabs = num + "\pi";
    tlabs(num==0) = "0";
    tlabs(num==1) = "\pi";
end

% For each axis (i.e. for x, or y, or both)
for a = 1:numel(axFmt)
    
    % Format the axis itself
    axsub = axFmt(a);
    axsub.Limits = P.limits*P.nreps;

    if isprop(axsub, "TickValues")
        axsub.TickValues = ticks;
    else
        axsub.Ticks = ticks;
    end

    if ~isempty(label)
        axsub.Label.String = label;
    end

    if outUnits == "radians"
        axsub.TickLabels = tlabs;
    end

    % Now format the corresponding data values
    for prop = propsFmt
        inds = repdata(dataObjs, prop, P.nreps, P.limits, toRadFcn, outUnits);
    end

    for prop = propsNoFmt
        repdata(dataObjs, prop, P.nreps, "none", toRadFcn, outUnits, inds);
    end

end

end


function [inUnits, axFmt, axNoFmt, dataObjs, propFmt, propNoFmt] = parseAxes(h, xy)

% if h.Type == "axes"
    if xy == "both"
        axFmt = [h.XAxis, h.YAxis];
        propFmt = ["XData", "YData"];
    elseif xy == "x"
        axFmt = h.XAxis;
        propFmt = ["XData"];
    elseif xy == "y"
        axFmt = h.YAxis;
        propFmt = ["YData"];
    end
    axNoFmt = setdiff([h.XAxis, h.YAxis], axFmt);
    propNoFmt = setdiff(["XData", "YData", "CData", "ZData"], propFmt);
% elseif h.Type == "colorbar"
%     axFmt = h;
%     axNoFmt = [];
%     propFmt = [];
%     propNoFmt = [];
% end

minDegRange = 175; % assume that range of 180+ indicates degrees
axranges = arrayfun(@(a) diff(a.Limits), axFmt);

if all(axranges > minDegRange)
    inUnits = "degrees";
else
    inUnits = "radians";
end

dataObjs = [];
for type = ["line", "bar", "scatter"] % doesn't work with image or patch data
    dataObjs = [dataObjs; findobj(h, "type", type)];
end

% hline = findobj(h, "type", "line");
% hbar = findobj(h, "type", "bar");
% hscatter = findobj(h, "type", "scatter");
% himg = findobj(h, "type", "image");
% dataObjs = [hline; hbar; hscatter];

end

function irepAll = repdata(dataObjs, prop, nreps, lim0, parseAngleFcn, outUnits, irepAll0)

if nargin < 7, irepAll0 = cell(size(dataObjs)); end

isCircData = ~strcmpi(lim0, "none");

if isCircData
    cyc = diff(lim0);
    dLow = lim0(1) * nreps; % lower limit
end

irepAll = {};

for n = 1:numel(dataObjs)

    h = dataObjs(n);

    if ~isprop(h, prop), continue; end
    a = h.(prop)(:);
    if isempty(a) || isscalar(a), continue; end

    if isCircData

        % Convert data to output units
        a = parseAngleFcn(a); % convert a to radians
        if outUnits == "degrees"
            a = rad2deg(a);
            wrapThresh = 180;
        else
            wrapThresh = pi;
        end
        ash = mod(a-dLow, cyc) + dLow; % all vars in output units
        % shift = ash-a;
        a = ash;
        [a, ~, ~, indsIn] = nanPadCircWrap(a, wrapThresh);
        indsIn(indsIn==0) = nan;
        itmp = indsIn;
%         ijump = find(diff(a))
%         if all(diff(a)>0)
%             doSort = false;
%         else
%             warning("%s values are not monotonically in increasing after applying mod()", prop);
%             doSort = true;
%         end
        
    else
        itmp = (1:numel(a))';
    end

    arep = [];
    irep = []; % is it OK that this is cleared with every obj? the output will correspond to the last obj

    indsSpecified = ~isempty(irepAll0{n});

    if indsSpecified
        irep = irepAll0{n};
        arep = nan(size(irep));
        v = ~isnan(irep);
        arep(v) = a(irep(v));
    else
        for r = 1:nreps
            atmp = a;
            if isCircData
                atmp = atmp + (r-1)*cyc;
            end
            arep = [arep; nan; atmp];
            irep = [irep; nan; itmp];
            %         ar = [ar; a+a0; nan]; % why insert the NaN? it breaks the line...
        end
        % [arep, isort] = sort(arep, "ascend");
        % irep = irep(isort);
        irepAll{n} = irep;
    end

    % if indsSpecified
    %     arep = arep(inds0{n});
    %     irep = irep(inds0{n});
    % elseif doSort
            % If the circular wrapping broke the continuity of the data,
            % restore it by sorting the values ascending. The other properties
            % will then be sorted according to the same indices.
            % [arep, isort] = sort(arep, "ascend");
            % irep = irep(isort);
    % end
% 
    dataObjs(n).(prop) = arep;

end

end