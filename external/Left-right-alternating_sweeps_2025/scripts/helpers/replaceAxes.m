function replaceAxes(ax, axNew)
% Replace the axes object AX within its parent view (either a figure or a
% tiledlayout), with the new axes object AXNEW
parent = ax.Parent;
if parent.Type == "tiledlayout"
    axPos = ax.Layout.Tile;
else
    axPos = ax.Position;
end
delete(ax);
axNew.Parent = parent;
if parent.Type == "tiledlayout"
    axNew.Layout.Tile = axPos;
else
    axNew.Position = axPos;
end
end