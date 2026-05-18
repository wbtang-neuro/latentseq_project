function stop = plotOptim(beta, ~, ~, betaInds, plotData, plotType, varargin)
%PLOTOPTIM plotting function for shift-model tuning curves

inp = inputParser();
inp.addParameter('box', []);
inp.addParameter('grid', {});
inp.parse(varargin{:});
P = inp.Results;

stop = false;

z = accumulate(plotData, beta(betaInds));
ax = gca;

col = [0, 0, 0] + 0.5;

titleStr = plotType;
plotType = string(plotType);
patts = ["alpha_", "beta_"];
strs = split(plotType, "_");
if startsWith(plotType, patts)
    plotType = strs(2);
    paramType = strs(1);
else
    paramType = "beta";
    plotType = strs(1);
end

if paramType == "beta"
    z = exp(z);
    valstr = "exp(\beta)";
    % ylims = [0, 2];
    lims = [0, 15];
elseif paramType == "alpha"
    valstr = "shift (m)";
    lims = [-0.2 0.2];
end

switch plotType
    
    case 'pos'
        h = findobj(ax, 'type', 'image');
        if isempty(h)
            if isempty(P.grid)
                imagesc(z);
            else
                imagesc(P.grid{2}, P.grid{1}, z');
            end
            axis(ax, 'equal', 'xy', 'off');
            hcb = colorbar(ax);
            title(hcb, valstr);
            if ~isempty(P.box)
                x = P.box(:, 1);
                y = P.box(:, 2);
                line(ax, x, y, 'color', col);
            end
            clim(lims)
        else
            h.CData = z';
        end
        
        % % Set clim according to data inside box only
        % if ~isempty(P.box) && ~isempty(P.grid)
        %     [xx,yy] = meshgrid(P.grid{:});
        %     pad = 0.2;
        %     x = P.box(:, 1) + pad*[1 1 -1 -1 1]';
        %     y = P.box(:, 2) + pad*[1 -1 -1 1 1]';
        %     v = inpolygon(xx(:), yy(:), x, y);
        %     zv = z(v);
        %     ax.CLim = prctile(zv, [0.5, 99.5]);
        % end
        
    case {'hd', 'id', 'theta'}
        h = findobj(ax, 'type', 'line');
        x = (1:numel(z))'/numel(z)*2*pi;
        x = x-pi; % added 2024-11-20 for consistency with [-pi, pi] range elsewhere
        x = [x; x+2*pi];
        y = [z; z];
        if isempty(h)
            axis(ax, 'equal', 'square');
            % Set limits and plot axes
            ax.XLim = [-pi, 3*pi];
            ax.XTick = pi*(-1:3);
            ax.XTickLabel = ["-\pi", "0", "\pi", "2\pi", "3\pi"];
            line(ax, x, y, 'color', col);
            xgrid = pi*[0, 1, 2];
            for n = 1:numel(xgrid)
                x = xgrid(n);
                xline(ax, x, 'color', [0, 0, 0]+0.6, 'handleVisibility', 'off');
            end
            yline(ax, 0, 'color', [0, 0, 0]+0.6, 'handleVisibility', 'off');
            ylabel(ax, valstr, "interpreter", "tex");
            ylim(ax, lims);
        else
            h.XData = x;
            h.YData = y;
        end
        checkYLim(ax, z);
        
    case {'postSpike', 'speed'}
        h = findobj(ax, 'type', 'line');
        if isempty(h)
            if isempty(P.grid)
                x = 1:size(z, 1);
            else
                x = P.grid{1};
            end
            if strcmpi(plotType, 'speed')
                line(ax, x, z, 'color', col, 'lineWidth', 2);
            elseif strcmpi(plotType, 'postSpike')
                bar(ax, x, z, 'FaceColor', col, 'EdgeColor', 'none', 'barWidth', 1);
            end
            ax.XLim = x([1, end]);
            line(ax, x([1, end]), [1, 1], 'color', [0 0 0]+0.6, 'lineStyle', '--', ...
                'handleVisibility', 'off');
            ylabel(ax, valstr, "interpreter", "tex");
            ylim(ax, lims);
        else
            h.YData = z;
        end
        checkYLim(ax, z);
    otherwise
        error("Unknown plot type '%s'", plotType);
end

if isempty(h)
    title(ax, titleStr, "interpreter", "none");
end

end

function [z, zAll] = accumulate(basis, weights)
sz = size(basis);
z = zeros([sz(1:end-1), 1]);
szz = size(z);
basis = reshape(basis, [prod(sz(1:end-1)) sz(end)]);
zAll = basis .* weights(:)';
zAll = reshape(zAll, sz);
z = sum(zAll, numel(sz));
end

function checkYLim(ax, y)
yl = ax.YLim;
yr = yl(2)- yl(1);
d = yr*0.5;
if min(y) < yl(1)
    yl(1) = yl(1) - d;
end
if max(y) > yl(2);
    yl(2) = yl(2) + d;
end
ax.YLim = yl;
end