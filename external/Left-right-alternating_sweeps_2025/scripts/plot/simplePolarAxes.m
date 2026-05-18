function ax = simplePolarAxes(parent)
if nargin < 1 || isempty(parent), parent = gcf; end

if isa(parent, "matlab.graphics.layout.TiledChartLayout")
    axsib = parent.Children;
    if isempty(axsib)
        tile = 1;
    else
        tile = max(arrayfun(@(a) a.Layout.Tile, axsib)) + 1;
    end
else
    tile = [];
end

ax = polaraxes(parent);
ax.RTick = [];
ax.ThetaTick = 90:90:360;
ax.ThetaTickLabel = [];
ax.ThetaZeroLocation = "top";

if ~isempty(tile)
    ax.Layout.Tile = tile;
end

end