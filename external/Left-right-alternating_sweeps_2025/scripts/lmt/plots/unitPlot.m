function unitPlot(fig, mdl, plotType, spikeTimes, iunit, plotData)

if nargin < 5 || isempty(iunit)
    iunit = 1:mdl.nunits;
end

persistent mdl0

if isempty(mdl0) || ~isfield(mdl0, "plotType") || (mdl0.(plotType) ~= mdl)
    % if current model is different from the last-used model, draw
    % everything from scratch (purge the handles)
    getPlotHandles(); 
    clf(fig);
end
    

tag = sprintf("unit_%s", plotType);
H = getPlotHandles(mdl, tag);

MAX_NPLOT = 12;
MAX_SPIKES_NPLOT = 1e4;
DITHER = 0.005; % sigma for display dither

scatterArgs = {"filled", "markerEdgeColor", "none", "markerFaceAlpha", 0.5};


if ~checkHandle(H, "button", 1)
    H.button = uicontrol("parent", fig, "style", "pushbutton", ...
        "units", "normalized", "position", [0.05, 0.8, 0.05, 0.05], ...
        "string", "REFRESH", "backgroundColor", [0, 0, 0]+0.3, ...
        "foregroundColor", [1, 1, 1], "fontSize", 16, ...
        "callback", @pushbuttonCallback);
end
holdUnits = lower(H.button.String) == "refresh";
if holdUnits && isfield(H, "lastUnitInds")
    % Keep units from the last draw
    iunit = H.lastUnitInds;
    if max(iunit) < mdl.nunits
        holdUnits = false;
    end
else
    % Randomly select a subset of units to plot
    if numel(iunit) > MAX_NPLOT
        irand = randperm(numel(iunit), MAX_NPLOT);
        iunit = sort(iunit(irand));
    end
    H.lastUnitInds = iunit;
end

nplot = numel(iunit);
aspect = 1;
splx = ceil(aspect * sqrt(nplot));
sply = ceil(nplot / splx);


nt = mdl.nt;
ndims = mdl.ndims;
dt = plotData.dt;
d = DITHER * (rand(nt, ndims, "like", mdl.Xinit) - 0.5);
X{1} = gather(mdl.Xinit + d);
X{2} = gather(mdl.XAligned + d);
    
lim = [min(X{2}); max(X{2})];
gv = mdl.Fgridv;
gsz = cellfun(@(s) numel(gather(s)), gv);

labelstring = "";

F = gather(mdl.logTuningCurvesF());
rng = gather(prctile(F(:), [1, 99]));
% flim = 1.2*[-1, 1]*abs(max(rng));
flim = mean(rng) + 0.6*[-1, 1]*diff(rng);
flimClip = rng;
% idat = mdl.getIterData();

for n = 1:nplot
    iu = iunit(n);
    spikes = gather(find(mdl.Y(:, iu) > 0));
    
    if ~checkHandle(H, "axes", n)
        ax = subplot(sply, splx, n, "parent", fig);
        ax.Title.FontSize = 5;
        ax.Title.HorizontalAlignment = "left";
        ax.Title.Units = "normalized";
        ax.Title.Position(1) = 0;
        ax.Title.Interpreter = "none";
        hold(ax, "on");
        H.axes(n) = ax;
    else
        ax = H.axes(n);
    end
    
    switch plotType
        case "pos"
            hasTrue = checkHandle(H, "true", n);
            hasLatent = checkHandle(H, "latent", n);
            if ~hasLatent || ~hasTrue
                hold(ax, "on");
                H.true(n) = plot3(ax, nan, nan, nan, "k.", "markerSize", 2);
                H.latent(n) = scatter3(ax, nan, nan, nan, 5, nan, scatterArgs{:});
                ax.Color = "none";
                ax.XAxis.Visible = "off";
                ax.YAxis.Visible = "off";
                ax.CLim = [-pi, pi];
                ax.Colormap = hsv();
                axis(ax, "equal");
                view(ax, [0, 89]);
            end
            ax.XLim = lim(:, 1);
            ax.YLim = lim(:, 2);
            
            if numel(spikes) > MAX_SPIKES_NPLOT
                spikes = randsample(spikes, MAX_SPIKES_NPLOT);
            end
            hlines = [H.true(n), H.latent(n)];
            for i=1:2
                hln = hlines(i);
                if ~isempty(X{i})
                    spos = X{i}(spikes, :);
                    stheta = plotData.theta(spikes);
                    hln.XData = spos(:, 1);
                    hln.YData = spos(:, 2);
                    hln.ZData = rand(size(spos, 1), 1);
                end
                H.latent(n).CData = stheta;
            end
            
        case "acorr"
            st = spikeTimes{iu};
            [c,b] = rg.signal.xCorrPointProcess(st, st, 0.02, 51);
            c(b==0) = 0;
            if checkHandle(H, "acg", n)
                H.acg(n).YData = c;
            else
                H.acg(n) = bar(ax, c, "barWidth", 1,  ...
                    "edgeColor", "none", "faceColor", 'k');
                axis(ax, "tight", "off");
            end
            
        case {"tc", "tuningcurve"}
            f = gather(F(:, iu));
            f = reshape(f, gsz);
            if ndims == 1
                if checkHandle(H, "tc_line", n)
                    H.tc_line(n).YData = f;
                else
                    H.tc_line(n) = line(ax, gv{1}, f, "color", 'k');
                end
                ax.YLim = flim;
            elseif ndims >= 2
                if ndims > 2
                    f = f(:, :, :); % squeeze all dims from 3+ into 3rd dim
                    f = squeeze(mean(f, 3));
                end
                if checkHandle(H, "tc_image", n)
                    H.tc_image(n).CData = f';
                else
                    H.tc_image(n) = imagesc(ax, gv{1:2}, f');
                    axis(ax, "xy", "image", "off");
                end
                ax.CLim = flimClip;
            end
            
            labelstring = sprintf("%.1f %.1f", flimClip);
            
        otherwise
            error("Invalid plot type '%s'", plotType);
    end
    
    mnrate = mean(mdl.Y(:, 1)/dt);
    nid = string(mdl.unitIds(iu));
    
    baselabel = sprintf("%s, %.1f Hz", nid, mnrate);
    if labelstring == ""
        str = baselabel;
    else
        str = sprintf("%s\n[ %s ]", baselabel, labelstring);
    end

    ax.Title.String = char(str);
    
end

getPlotHandles(mdl, tag, H);

mdl0.(plotType) = mdl;

end

function pushbuttonCallback(src, ~)
button = src;
currentMode = lower(button.String);
if currentMode == "refresh"
    newMode = "hold";
elseif currentMode == "hold"
    newMode = "refresh";
end  
button.String = upper(newMode);
end