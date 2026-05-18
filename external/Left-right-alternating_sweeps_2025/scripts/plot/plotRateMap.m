function h = plotRateMap(rmaps, varargin)

S = SweepsSettings();
inp = inputParser();
inp.addParameter("colors",  S.col_example_cell());
inp.addParameter("showContours", false);
inp.addParameter("contourLevels", 85);
inp.addParameter("lineWidth", 1);
inp.addParameter("axes", []);
inp.addParameter("limits", 0.75*[-1, 1, -1, 1]);
inp.addParameter("showLabels", true);
inp.parse(varargin{:});
P = inp.Results;

ax = P.axes;
if isempty(ax), ax = gca; end

nrm = numel(rmaps);
rm0 = rmaps(1);
x = rm0.bins{1};
    
if rm0.ndim==1
    % Angular RM
    rho = [rmaps.z];
    rho = rho ./ max(rho);
    inds = [1:numel(x), 1];
    theta = x(inds);
    rho = rho(inds, :);

    if ax.Type == "axes"
        axNew = polaraxes();
        replaceAxes(ax, axNew)
        ax = axNew;
    end
    h.polar = polarplot(ax, theta, rho);
    if P.showLabels
        ax.ThetaTick = (1:4)*90;
        ax.ThetaTickLabel = string(rem(ax.ThetaTick, 360));
    end
    ax.RTick = [];
    for r = 1:nrm
        h.polar(r).LineWidth = P.lineWidth;
        h.polar(r).Color = P.colors(r, :);
    end
elseif rm0.ndim==2
    y = rm0.bins{2};
    % 2D-POS RM
    v = rm0.validBin';
    if nrm == 1
        z = rm0.z';
    else
        z = ones([size(rm0.z), 3]);
        exponents = [1.5, 1.5];
        for r = 1:nrm
            rm = rmaps(r);
            clipRange = prctile(rm.z(:), [20 98]);
            z = rg.color.subtractMappedColors( ...
                z, rm.z', P.colors(r, :), clipRange, exponents(r));
        end
    end
    z(~v) = nan;
    h.image = imagesc(ax, rm0.bins{:}, z, "alphaData", v);
    if P.showContours
        for r = 1:nrm
            rm = rmaps(r);
            z = rm.z';
            z(~rm.validBin') = nan;
            clevels = prctile(z(:), P.contourLevels);
            if isscalar(clevels)
                clevels = clevels.*[1, 1];
            end
            [~, h.contour(r)] = contour(ax, rm.bins{:}, z, clevels, ...
                'color', P.colors(r, :), 'lineWidth', P.lineWidth);
        end
    end
    ax.XAxis.Visible = "off";
    ax.YAxis.Visible = "off";
    ax.XLim = P.limits([1, 2]);
    ax.YLim = P.limits([3, 4]);
    [xx,yy] = meshgrid(x,y);
    vlim = restrictq(xx, ax.XLim) & restrictq(yy, ax.YLim);
    v = v&vlim;
    ax.CLim = [0, prctile(z(v), 99)];
    ax.DataAspectRatio = [1, 1, 1];
end

h.axes = ax; % axis obj may change if we convert normal 2d axes to polaraxes

end