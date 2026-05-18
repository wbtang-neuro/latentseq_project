function lmt2DPosPlots(varargin)

fig = varargin{1};
mdl = varargin{2};
plotType = varargin{3};
plotData = varargin{4};

H = getPlotHandles(mdl, plotType);

X = gather(mdl.X);
nsmooth = plotData.nsmooth;
if nsmooth
    if mdl.isCircular
        xfcn = @(x) gather(gsmoothcirc(x, nsmooth));
    else
        xfcn = @(x) gather(gsmooth(x, nsmooth));
    end
else
    xfcn = @(x) gather(x);
end

% XOld = xfcn(mdl.stateLastIter.X);

Xbnd = mdl.FgridBounds;

Xinit = xfcn(mdl.Xinit);
XA = xfcn(mdl.XAligned);
% XOldA = mdl.alignX(XOld);
XbndA = sort(gather(mdl.alignX(Xbnd')), 1);

if mdl.isCircular
    cen = 0;
else
    cen = mean(XA);
end
ndims = size(X, 2);
nt = mdl.nt;
dt = plotData.dt;

maxPoints = 1e5;
if nt > maxPoints
    isub = round(linspace(1, nt, maxPoints));
else
    isub = (1:nt)';
end

% scheme = figColScheme2;
scatterArgs = {"filled", "markerEdgeColor", "none", "markerFaceAlpha", 0.6};

switch lower(plotType)
    case "x_time"
        
        clf(fig);
        if isfield(plotData, "timeBins")
            ip = plotData.timeBins;
        else
            % plot the segment with highest theta amplitude
            plotlen = 10;
            nplt = 2*ceil(plotlen/2/dt)+1; % select odd number of points
            isegcenter = (nplt+1)/2;
            spd = movmean(plotData.speed, nplt);
            [~, imx] = max(spd);
            % ensure plotted segment is within data bounds
            i0 = max(0, imx-isegcenter);
            i0 = min(i0, nt-isegcenter+1);
            ip = i0 + (1:nplt)';
        end
        
        tplt = mdl.tgrid(ip);
        ndims = mdl.ndims;
        
        for di = 1:ndims
            ax = subplot(ndims,1,di, 'parent', fig);
            hold(ax, "on");
            cend = cen(di);
            makePlot = @(y, varargin) plotLocal(mdl, ax, tplt, y-cend, varargin{:});
            h = {};
            names = string([]);
            h{1} = makePlot(Xinit(ip,di), 'color', 'r', "lineWidth", 1);
            names(1) = "XInit";
            
            % Plot summed X
            h{end+1} = makePlot(XA(ip,di), 'color', [0, 0.7, 0]);
            names(end+1) = "X";

            % Previous X
            % h{end+1} = makePlot(XOldA(ip,di), ':', 'linewidth', 0.5, 'color', 0.7*[1, 1, 1]);
            % names(end+1) = "last X";
            
            yrng = XbndA(:, di) - cend;
            
            if mdl.isCircular
                ax.YLim = [0, 4*pi];
                ax.YTick = (0:4)*pi;
                ax.YTickLabel = string(0:4) + " \pi";
                ax.TickLabelInterpreter = "tex";
            else
                ax.YLim = yrng;
            end
            ax.XLim = tplt([1, end]);
            col = [1, 1, 1]*0.7;
            ax.XColor = col;
            ax.YColor = col;
            xlabel(ax, "Time / s", "color", col);
            ax.Color = "none";
            title(ax, sprintf("X dim %u", di));
        end
        
        h = [h{:}];
        legend(h, names);
        
    case "f_peaks"
        
        newplt = ~checkHandle(H, "f_peaks_path", mdl.ndims);
        R = gather(mdl.logTuningCurvesF());
        mx = max(R);
        % Fcen = gather(mdl.FcenterOfMass());
        [Fpk, ipk] = max(R);
        pkpos = gather(mdl.Fgrid(ipk, :));
        sigma = 0.01 * diff(mdl.FgridBounds, [], 2)';
        pkpos = pkpos + sigma .* randn(size(pkpos));
        pkpos = gather(pkpos);
        p = XA(isub, :) - cen;
        
        if newplt
            clf(fig)
            ax = subplot(2, 2, 1, "parent", fig);
            hold(ax, "on");
            H.f_peaks_path = plot(ax, nan, nan, 'k.', 'markerSize', 0.5);
            H.f_peaks_peaks = scatter(ax, nan, nan, 50, nan, "filled", "markerEdgeColor", "w", "lineWidth", 1);
            title(ax, "Peak positions");
            axis(ax, "equal", "off");
            dtip = dataTipTextRow("id", mdl.unitIds);
            H.f_peaks_peaks.DataTipTemplate.DataTipRows(end+1) = dtip;

            ax = subplot(2, 2, 2, "parent", fig);
            hold(ax, "on");
            H.f_peaks_histogram = imagesc(ax, mdl.Fgridv{:}, nan(mdl.nf));
            axis(ax, "image", "xy", "off");
            title(ax, "Grid occupancy histogram");
            
            
            % F average and variability
            ax = subplot(2, 2, 3, "parent", fig);
            H.f_peaks_meanimg = imagesc(ax, mdl.Fgridv{:}, nan(mdl.nf));
            axis(ax, "xy", "image", "off");
            title(ax, "Tuning density");
            
            ax = subplot(2, 2, 4, "parent", fig);
            H.f_peaks_stdimg = imagesc(ax, mdl.Fgridv{:}, nan(mdl.nf));
            axis(ax, "xy", "image", "off");
            title(ax, "Tuning std dev.");
        end
        
        % Update path
        h = H.f_peaks_path;
        h.XData = p(:, 1);
        h.YData = p(:, 2);
        
        % Update peak pos scatter
        h = H.f_peaks_peaks;
        h.XData = pkpos(:, 1) - cen(1);
        h.YData = pkpos(:, 2) - cen(2);
        h.CData = mx;
        ax = h.Parent;
        gg = gather(mdl.Fgrid);
        ax.XLim = [min(gg(:, 1)), max(gg(:, 1))];
        ax.YLim = [min(gg(:, 2)), max(gg(:, 2))];
        
        h = H.f_peaks_meanimg;
        z = mean(R, 2);
        h.CData = reshape(z, mdl.nf, mdl.nf);
        ax = h.Parent;
        lim = prctile(z(:), [1, 99]);
        ax.CLim = lim;
        title(ax, sprintf("Range = %.2f", diff(lim)));
        
        h = H.f_peaks_stdimg;
        z = std(R, [], 2);
        h.CData = reshape(z, mdl.nf, mdl.nf);
        ax = h.Parent;
        ax.CLim = [0, prctile(z(:), 99)];

        h = H.f_peaks_histogram;
        edges = cellfun(@(x) {cen2edg(gather(x))}, mdl.Fgridv);
        Xhist = histcounts2(gather(X(:, 1)), gather(X(:, 2)), edges{:});
        h.CData = Xhist;
        h.Parent.CLim = [0, gather(prctile(Xhist(:), 99))];
        
        
    case "x_2d"
        
        thetaPhase = plotData.theta;
        hd = plotData.hd;
        posOffset = vecnorm(mdl.XAligned - mdl.Xinit, 2, 2);
        
        try
            [sweepTheta, sweepRho] = getSweepDir(plotData.models.id, "id");
            h = wrapTo2Pi(sweepTheta) ./ (2*pi);
            v = sweepRho ./ prctile(sweepRho, 80);
            v = min(v, 1);
            colhsv = [h, ones(nt, 1), v];
            sweepThetaRgb = hsv2rgb(colhsv);
        catch e
            sweepThetaRgb = zeros(mdl.nt, 3);
            sweepThetaRgb(:, 1) = 1;
        end
        
        coldat = {thetaPhase, sweepThetaRgb, hd, posOffset};
        cmaps = {hsv, hsv, hsv, viridis};
        clims = {[-pi, pi], [], [-pi, pi], []};
        names = ["theta", "sweep", "hd", "shift magnitude"];
        
        H.axScatterLink = [];

        % grid coordinates box
        cen = gather(mean(mdl.X));
        gv = mdl.Fgridv;
        limx = gather(gv{1}([1, end]));
        limy = gather(gv{2}([1, end]));
        bx = limx([1, 1, 2, 2, 1]) - cen(1);
        by = limy([1, 2, 2, 1, 1]) - cen(2);
        
        for n = 1:numel(names)
            newplt = ~checkHandle(H, "scatter", n);
            p = XA(isub, :) - cen;
            p(:, 3) = 0;
            ctmp = gather(coldat{n});
            if ~isempty(ctmp)
                col = ctmp(isub, :);
            end
            
            clim = clims{n};

            if newplt
                ax = subplot(2, 2, n, "parent", fig);

                hold(ax, "on");
                H.path(n) = plot(ax, nan, nan, "k", "lineWidth", 0.3);

                H.axScatter(n) = ax;
                H.scatter(n) = scatter(ax, nan, nan, 2, scatterArgs{:});
                resizeobj(ax, 1.3);
                if ~isempty(clim)
                    ax.CLim = clim;
                end
                colormap(ax, cmaps{n});
                title(ax, names(n));
                if ndims==3
                    axis(ax, "off", "vis3d");
                    rotate3d(ax, 'on');
                    ax.Clipping = 'off';
                elseif ndims==2
                    axis(ax, "equal", "off")
                end
                ax.Clipping = "off";
                colorbar(ax);

                H.box(n) = plot(ax, nan, nan, "k");

            end
            
            h = H.scatter(n);
            h.XData = p(:, 1);
            h.YData = p(:, 2);
            h.ZData = p(:, 3);
            h.CData = col;

            h = H.path(n);
            h.XData = XA(:, 1) - cen(1);
            h.YData = XA(:, 2) - cen(2);

            H.box(n).XData = bx;
            H.box(n).YData = by;
            
            ax = H.axScatter(n);

            if isempty(clim)
                clim = gather(prctile(col(:), [1, 99]));
                if clim(2) > clim(1)
                    ax.CLim = clim;
                end
            end

        end
        
        if newplt
            H.axScatterLink = linkprop(H.axScatter, "CameraPosition");
        end
        
end

getPlotHandles(mdl, plotType, H);

end

function [theta, rho] = getSweepDir(mdl, method)

Xpos = double(gather(mdl.XAligned));
if method == "offset"
    dpos = Xpos - mdl.Xinit;
    dpos = double(gather(dpos));
elseif method == "filter"
    fs = 1/mdl.dt;
    [b, a] = butter(2, 2/(fs/2), 'high'); % HPF
    dpos = filtfilt(b, a, Xpos);
end
[theta, rho] = cart2pol(dpos(:, 1), dpos(:, 2));

end

function h = plotLocal(mdl, ax, x, y, varargin)

if mdl.isCircular
    x = [x; nan; x];
    y = wrapTo2Pi(y);
    y = [y; nan; y+2*pi];
    [x, y] = nanPadCircWrapXy(x, y);
end

h = plot(ax, x, y, varargin{:});

end