function h = sweepIdTrajectoryPlotShaded(ax, t, x, y, trng, kwargs)

arguments
    ax
    t
    x
    y
    trng
    kwargs.col = "odd_even";
    kwargs.sweeps = []
    kwargs.id = []
    kwargs.box = "circle1.5"
    kwargs.plotPath = true;
    kwargs.plotScaleBar = false;
end
sweeps = kwargs.sweeps;
chk = kwargs.id;
col = kwargs.col;

if nargin < 5, trng = []; end

S = SweepsSettings();

% plot trajectory
vt = restrictq(t, trng);
if kwargs.plotPath
h.pos = plot(ax, x(vt), y(vt), "lineWidth", 1.5, "color", S.col_pos_true);
end

%% Plot sweeps
theta = linspace(0, 2*pi, 100);
rho = normpdf(theta, pi, .4);
rho = circshift(rho, 50);
clear swppatch
[swppatch(:, 1), swppatch(:, 2)] = pol2cart(theta, rho);
swppatch(end+1, :) = swppatch(1, :);
swppatch = swppatch*.1;
%%
% Apply time range
if ~isempty(sweeps)
    if isfield(sweeps, "pos")
       vsweep = [sweeps.tStart] > trng(1);% & [sweeps.tStop] < trng(2);
       sweeps = sweeps(vsweep);
    
        % Parse plotting colors
        nsweeps = numel(sweeps);
        if (isnumeric(col) && numel(col)==3) || strlength(col) == 1
            % numeric RGB values or single-character color code
            cols = repmat({col}, nsweeps, 1);
        elseif strcmpi(col, "odd_even")
            [cols(1:2:nsweeps)] = deal({S.col_cyc_odd});
            [cols(2:2:nsweeps)] = deal({S.col_cyc_even});
        end
          
        for s = 1:numel(sweeps)
            
            if isfield(sweeps, "pos")
                sweep = sweeps(s);
                p = sweep.pos(2:end, :);
                p = gsmooth(p, 1);
                dx = p(end, 1)-p(2,1);
                dy = p(end, 2)-p(2,2);
                sweepdir = atan2(dy, dx);
            elseif isfield(sweep.sweepdir)
                sweepdir = sweepdir(s);
            end
            thisSwppatch = rotate2d(swppatch, +sweepdir);
            thisSwppatch = thisSwppatch+[x(sweep.iStart), y(sweep.iStart)];
    
            h.sweeps(s) = patch(ax, thisSwppatch(:, 1), thisSwppatch(:, 2), cols{s}, "FaceAlpha", .5, "EdgeColor", "none");
     end
    else
        vcyc = restrictq(sweeps.tStart, trng);
        nsweeps = sum(vcyc);
        sweepdir = sweeps.sweepdir(vcyc);
        iStart = sweeps.iStart(vcyc);
        if (isnumeric(col) && numel(col)==3) || strlength(col) == 1
            % numeric RGB values or single-character color code
            cols = repmat({col}, nsweeps, 1);
        elseif strcmpi(col, "odd_even")
            [cols(1:2:nsweeps)] = deal({S.col_cyc_odd});
            [cols(2:2:nsweeps)] = deal({S.col_cyc_even});
        end
        for c = 1:nsweeps
            thisSwppatch = rotate2d(swppatch, +sweepdir(c));
            thisSwppatch = thisSwppatch+[x(iStart(c)), y(iStart(c))];
            h.sweeps(c) = patch(ax, thisSwppatch(:, 1), thisSwppatch(:, 2), cols{c}, "FaceAlpha", .5, "EdgeColor", "none");
        end
    end
end
%% Plot ID
if ~isempty(chk)
    vcyc = restrictq(chk.tStart, trng);
    ncyc = sum(vcyc);
    id = chk.id(vcyc);
    pkinds = chk.iStart(vcyc);
    [u, v] = pol2cart(id, 0.18);
    hold on
    clear hId
    for c = 1:ncyc
        xq = x(pkinds(c)) + [0, u(c)];
        yq = y(pkinds(c)) + [0, v(c)];
        hId(c) = plot(xq, yq, "lineWidth", 1.1, "Color",S.col_id);
    end
    
    set(hId(1:2:ncyc), "Color", S.col_cyc_odd);
    set(hId(2:2:ncyc), "Color", S.col_cyc_even);
end

% sweeps_plot_box(ax);
axis(ax, "equal", "off", "tight");

if kwargs.plotScaleBar
    x0 = ax.XLim(1);
y0 = ax.YLim(1);
ax.Clipping = "off";
ax.XLim = ax.XLim;
ax.YLim = ax.YLim;
h.scale_bar = plot(ax, x0+[0, 0.5], y0 - 0.1 + [0, 0], "k", "lineWidth", 2);
end


