function Pwrapped = wrapPosToGridTile(P, gridAxes, tileShape)
% Wraps 2D position coordinates onto a unit rhombus defined by two adjecent
% grid axes
% Inputs:
% P [nt x 2]
% grid axes [x, y]
% tileShape either "rhombus" or "hexagon"
% 
% Use 2 axes for wrapping to rhombus, or 3 axes for wrapping to hexagon
% P is transformed into rhombus/hex space by using the basis vectors
% Pwrapped is transformed back to 2D space

if nargin < 3 || isempty(tileShape), tileShape = "hexagon"; end

np = size(P, 1);

if strcmpi(tileShape, "rhombus")
    assert(size(gridAxes, 1)==2, "Wrapping to rhombus tile requires two grid axes.");
    [~, Pwrapped] = wrapOneAxesPair(P, gridAxes, tileShape);
    % Pn = P/gridAxes;
    % Pnw = mod(Pn, 1);
    % Pw = Pnw*gridAxes;
elseif strcmpi(tileShape, "hexagon")

    % To find the hexagonal phase tile, we iterate through the three 
    % possible rhombuses, each of which is defined by a pairwise 
    % combination of two of the three grid axes.
    %
    % We know that the axes are supplied in counter-clockwise order, so we
    % can use a fixed set of combinations of axes, with each combination
    % being 60° apart:
    %   Rhombus 1: axes 1, 2
    %   Rhombus 2: axes 2, 3
    %   Rhombus 3: axes 3, -1 (inverted)
    axPairs = {[1, 2], [2, 3], [3, -1]};
    assert(size(gridAxes, 1)==3, "Wrapping to hexagonal tile requires three grid axes.");

    for a = 1:3
        iax = abs(axPairs{a});
        axSign = sign(axPairs{a});
        gridAxPair = gridAxes(iax, :) .* axSign(:);
        [~, Pwrapped(:, :, a), rad(:, a)] = wrapOneAxesPair(P, gridAxPair, tileShape);
    end

    % Now check which of the three possible wrappings is nearest to [0, 0].
    % Choosing the nearest point corresponds to wrapping on the hexagonal 
    % tile.
    [~, ibestAxPair] = min(rad, [], 2);
    [ii, jj] = ndgrid(1:np, 1:2);
    kk = repmat(ibestAxPair, 1, 2);
    inds = sub2ind(size(Pwrapped), ii(:), jj(:), kk(:));
    Pwrapped = reshape(Pwrapped(inds), [np, 2]);
end

end

function [Pn, Pw, r] = wrapOneAxesPair(P, gridAxes, tileShape)
    %Pn normalized position on one 
    Pn = P/gridAxes;
    if strcmpi(tileShape, "rhombus")
        % zero phase is corner of rhombus
        Pn = mod(Pn, 1);
    elseif strcmpi(tileShape, "hexagon")
        % zero phase is center of hexagon
%         Pn = mod(Pn, 1);
        Pn = mod(Pn+0.5, 1) - 0.5;
    end
    Pw = Pn*gridAxes;
    r = hypot(Pw(:, 1), Pw(:, 2));
end