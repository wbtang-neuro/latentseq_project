function h = diagline(varargin)

args = varargin;
range = [];
ax = [];

if nargin
    % Axes will always be first arg, if specified
    if isgraphics(args{1}, 'axes')
        ax = args{1};
        args(1) = [];
    end

    % Range is next arg after axes
    if isnumeric(args{1})
        range = args{1};
        args(1) = [];
    end
end

if isempty(ax)
    ax = gca;
end

if isempty(range)
    xl = ax.XLim;
    yl = ax.YLim;
    lims = [xl; yl];
   range = [min(lims(:, 1)), max(lims(:, 2))];
end

ax.XLim = range;
ax.YLim = range;
h = plot(ax, range, range, args{:});

end